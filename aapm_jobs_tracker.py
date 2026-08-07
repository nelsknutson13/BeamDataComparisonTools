#!/usr/bin/env python3
from __future__ import annotations
"""
AAPM Jobs Tracker
Scrapes careers.aapm.org/jobs (live + Wayback Machine) and builds a local
SQLite database for medical physics workforce demand analysis.

Dependencies (install once):
    conda run pip install selenium webdriver-manager
"""

import json
import re
import sqlite3
import subprocess
import threading
import time
import tkinter as tk
import urllib.request
from tkinter import ttk, scrolledtext, messagebox
from datetime import datetime, timedelta
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("TkAgg")
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    import pandas as pd
    _HAS_CHARTS = True
except ImportError:
    _HAS_CHARTS = False

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.edge.options import Options as EdgeOptions
from selenium.webdriver.common.by import By
from selenium.common.exceptions import NoSuchElementException, WebDriverException, TimeoutException
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# ── Constants ─────────────────────────────────────────────────────────────────
HERE       = Path(__file__).parent
DB_PATH    = HERE / "aapm_jobs.db"
DEBUG_HTML = HERE / "aapm_debug_page.html"
FETCH_DAYS = 14
LIVE_BASE  = "https://careers.aapm.org/jobs"
WB_BASE    = "https://web.archive.org/web/{ts}/https://careers.aapm.org/jobs/"

# HTTP-200 Wayback snapshots (collected 2026-08-05)
WAYBACK_200 = [
    "20150627010654", "20161025135327", "20161103041750", "20170113105932",
    "20170612112529", "20180210221859", "20180731113913", "20190718145820",
    "20190822141459", "20190920180003", "20191022065932", "20191113040048",
    "20191212205210", "20200322071940", "20200713060732", "20201001114451",
    "20201114180716", "20210119042942", "20210413122912", "20210729110250",
    "20210918185204", "20220119034216", "20220302131826", "20220619104514",
    "20220701155302", "20221129065423", "20230330195642", "20230501221007",
    "20230601043538", "20230701114247", "20230807055424", "20230926164404",
    "20231201162638", "20240102035615", "20240208050802", "20240301044926",
    "20240501155132", "20240814104916", "20241107085239", "20250125005916",
    "20250213182247",
    "20250311034536", "20250401104031", "20250501122054", "20250601090013",
    "20250701102509", "20250801101134", "20250901095306", "20251001091846",
    "20251103214352", "20251201094311", "20260101084109", "20260201081110",
    "20260613151026",
]

_RANK_PATTERNS = [
    ("Resident",            r"\bresident\b"),
    ("Fellow",              r"\bfellow\b"),
    ("Assistant Professor", r"\bassistant\s+professor\b"),
    ("Associate Professor", r"\bassociate\s+professor\b"),
    ("Professor",           r"\bprofessor\b"),
    ("Chief Physicist",     r"\bchief\s+(medical\s+)?physicist\b"),
    ("Director",            r"\bdirector\b"),
    ("Senior Physicist",    r"\bsenior\s+(medical\s+)?physicist\b"),
    ("Lead Physicist",      r"\blead\s+(medical\s+)?physicist\b"),
    ("Staff Physicist",     r"\bstaff\s+(medical\s+)?physicist\b"),
]

# ── Database ──────────────────────────────────────────────────────────────────
def _conn() -> sqlite3.Connection:
    c = sqlite3.connect(DB_PATH)
    c.row_factory = sqlite3.Row
    return c


def init_db():
    with _conn() as conn:
        conn.executescript("""
            CREATE TABLE IF NOT EXISTS jobs (
                job_id        TEXT PRIMARY KEY,
                title         TEXT,
                employer      TEXT,
                city          TEXT,
                state         TEXT,
                date_posted   TEXT,
                date_closes   TEXT,
                job_type      TEXT,
                position_type TEXT,
                rank_level    TEXT,
                salary_min    REAL,
                salary_max    REAL,
                salary_text   TEXT,
                url           TEXT,
                first_seen    TEXT,
                last_seen     TEXT,
                times_seen    INTEGER DEFAULT 1
            );
            CREATE TABLE IF NOT EXISTS fetch_log (
                id          INTEGER PRIMARY KEY AUTOINCREMENT,
                fetched_at  TEXT NOT NULL,
                source      TEXT NOT NULL,
                jobs_found  INTEGER DEFAULT 0,
                new_jobs    INTEGER DEFAULT 0
            );
            CREATE TABLE IF NOT EXISTS settings (
                key   TEXT PRIMARY KEY,
                value TEXT
            );
            CREATE TABLE IF NOT EXISTS match_stats (
                year                  INTEGER PRIMARY KEY,
                registered_applicants INTEGER,
                positions_offered     INTEGER,
                matches_filled        INTEGER,
                unmatched_applicants  INTEGER
            );
        """)
        # Seed approximate values read from the published MedPhys Match chart.
        # User should verify/correct from the official annual PDFs.
        _MATCH_SEED = [
            (2017, 290, 110, 105, 105),
            (2018, 270, 125, 120,  85),
            (2019, 270, 135, 130,  75),
            (2020, 265, 140, 130,  75),
            (2021, 280, 140, 125,  80),
            (2022, 265, 140, 130,  80),
            (2023, 310, 150, 140,  90),
            (2024, 320, 165, 155,  95),
            (2025, 350, 180, 170, 105),
            (2026, 430, 195, 195, 145),
        ]
        with _conn() as conn:
            conn.executemany(
                "INSERT OR IGNORE INTO match_stats "
                "(year,registered_applicants,positions_offered,matches_filled,unmatched_applicants)"
                " VALUES (?,?,?,?,?)",
                _MATCH_SEED
            )


def upsert_jobs(jobs: list[dict], source: str) -> tuple[int, int]:
    seen_at = datetime.now().isoformat()
    new = 0
    with _conn() as conn:
        for job in jobs:
            if not job.get("job_id"):
                continue
            exists = conn.execute(
                "SELECT job_id FROM jobs WHERE job_id=?", (job["job_id"],)
            ).fetchone()
            if exists:
                conn.execute(
                    """UPDATE jobs SET last_seen=?, times_seen=times_seen+1,
                       title=COALESCE(?,title), employer=COALESCE(?,employer),
                       salary_text=COALESCE(?,salary_text),
                       salary_min=COALESCE(?,salary_min),
                       salary_max=COALESCE(?,salary_max),
                       position_type=COALESCE(?,position_type),
                       rank_level=COALESCE(?,rank_level),
                       job_type=COALESCE(?,job_type),
                       date_closes=COALESCE(?,date_closes)
                       WHERE job_id=?""",
                    (seen_at, job.get("title"), job.get("employer"),
                     job.get("salary_text"), job.get("salary_min"),
                     job.get("salary_max"), job.get("position_type"),
                     job.get("rank_level"), job.get("job_type"),
                     job.get("date_closes"), job["job_id"])
                )
            else:
                conn.execute(
                    """INSERT INTO jobs
                       (job_id,title,employer,city,state,date_posted,date_closes,
                        job_type,position_type,rank_level,salary_min,salary_max,
                        salary_text,url,first_seen,last_seen)
                       VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                    (job["job_id"], job.get("title"), job.get("employer"),
                     job.get("city"), job.get("state"), job.get("date_posted"),
                     job.get("date_closes"), job.get("job_type"),
                     job.get("position_type"), job.get("rank_level"),
                     job.get("salary_min"), job.get("salary_max"),
                     job.get("salary_text"), job.get("url"), seen_at, seen_at)
                )
                new += 1
        conn.execute(
            "INSERT INTO fetch_log (fetched_at,source,jobs_found,new_jobs) VALUES (?,?,?,?)",
            (seen_at, source, len(jobs), new)
        )
        if source == "live":
            conn.execute(
                "INSERT OR REPLACE INTO settings (key,value) VALUES ('last_live_fetch',?)",
                (seen_at,)
            )
    return len(jobs), new


def get_last_fetch() -> datetime | None:
    with _conn() as conn:
        row = conn.execute(
            "SELECT value FROM settings WHERE key='last_live_fetch'"
        ).fetchone()
    return datetime.fromisoformat(row["value"]) if row else None


def get_wayback_done() -> set[str]:
    with _conn() as conn:
        rows = conn.execute(
            "SELECT source FROM fetch_log WHERE source != 'live'"
        ).fetchall()
    return {r["source"] for r in rows}


def get_stats() -> dict:
    with _conn() as conn:
        total   = conn.execute("SELECT COUNT(*) FROM jobs").fetchone()[0]
        first   = conn.execute("SELECT MIN(first_seen) FROM jobs").fetchone()[0]
        wb_done = conn.execute(
            "SELECT COUNT(*) FROM fetch_log WHERE source != 'live'"
        ).fetchone()[0]
    return {"total": total, "first_seen": first, "wb_done": wb_done}


# ── Scraper helpers ───────────────────────────────────────────────────────────
def _infer_rank(title: str) -> str | None:
    tl = title.lower()
    for rank, pattern in _RANK_PATTERNS:
        if re.search(pattern, tl):
            return rank
    return None


def _parse_salary(text: str) -> tuple[float | None, float | None]:
    nums = [float(n.replace(",", "")) for n in re.findall(r"\d[\d,]*(?:\.\d+)?", text or "")]
    if len(nums) >= 2:
        return nums[0], nums[1]
    if len(nums) == 1:
        return nums[0], None
    return None, None




def _make_driver():
    ua = (
        "user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
    )
    common_args = [
        "--headless", "--no-sandbox", "--disable-dev-shm-usage",
        "--disable-gpu", "--window-size=1920,1080",
        "--disable-blink-features=AutomationControlled", ua,
    ]
    # Try Chrome first, fall back to Edge
    try:
        opts = Options()
        for a in common_args:
            opts.add_argument(a)
        opts.add_experimental_option("excludeSwitches", ["enable-automation"])
        return webdriver.Chrome(options=opts)
    except WebDriverException:
        opts = EdgeOptions()
        for a in common_args:
            opts.add_argument(a)
        return webdriver.Edge(options=opts)



def _table_fields(card) -> dict:
    """Extract labeled table rows (Salary, Position Type, dates, etc.)."""
    result = {}
    for row in card.find_elements(By.CSS_SELECTOR, "table tr"):
        cells = row.find_elements(By.TAG_NAME, "td")
        if len(cells) < 2:
            continue
        label = cells[0].text.strip().lower()
        value = cells[1].text.strip()
        if not value:
            continue
        if "salary" in label:
            result["salary_text"] = value
            result["salary_min"], result["salary_max"] = _parse_salary(value)
        elif "position type" in label or "employment" in label:
            result["job_type"] = value
        elif "clos" in label or "deadline" in label:
            result["date_closes"] = value
        elif "post" in label or "open" in label:
            result["date_posted"] = value
        elif "rank" in label or "level" in label:
            result["position_type"] = value
    return result


def _scrape_cards_new(driver) -> list[dict]:
    """Current AAPM careers site (JS-rendered, post-redesign)."""
    cards = driver.find_elements(
        By.CSS_SELECTOR, "ul.js-search-regular-slot li.search-item[data-id]"
    )
    if not cards:
        return []

    jobs = []
    for card in cards:
        job_id = card.get_attribute("data-id")
        if not job_id:
            continue

        job: dict = {"job_id": job_id}

        try:
            link = card.find_element(By.CSS_SELECTOR, "a.js-item-job-url")
            job["title"] = link.text.strip()
            job["url"]   = link.get_attribute("href") or ""
            if job["title"]:
                job["rank_level"] = _infer_rank(job["title"])
        except NoSuchElementException:
            continue

        loc_divs = card.find_elements(By.CSS_SELECTOR, "div.location div.small")
        if len(loc_divs) >= 1:
            job["employer"] = loc_divs[0].text.strip()
        if len(loc_divs) >= 2:
            loc = loc_divs[1].text.strip()
            parts = [p.strip() for p in loc.split(",")]
            job["city"]  = parts[0]
            job["state"] = parts[1] if len(parts) > 1 else None

        job.update(_table_fields(card))
        jobs.append(job)

    return jobs


def _scrape_cards_old(driver) -> list[dict]:
    """Pre-redesign AAPM careers site (server-side rendered, used in Wayback ~2015–2022)."""
    cards = driver.find_elements(
        By.CSS_SELECTOR, "div.bti-ui-job-result-wrapper"
    )
    if not cards:
        return []

    jobs = []
    for card in cards:
        try:
            link = card.find_element(By.CSS_SELECTOR, "a.bti-job-detail-link")
        except NoSuchElementException:
            continue

        job_id = link.get_attribute("job")
        if not job_id:
            continue

        job: dict = {"job_id": job_id}
        job["title"] = link.text.strip()
        job["url"]   = link.get_attribute("href") or ""
        if job["title"]:
            job["rank_level"] = _infer_rank(job["title"])

        # Location: "City, State" in div.bti-ui-job-result-detail-location
        try:
            loc_el = card.find_element(By.CSS_SELECTOR, "div.bti-ui-job-result-detail-location")
            loc = loc_el.text.strip()
            parts = [p.strip() for p in loc.split(",")]
            job["city"]  = parts[0] if parts else None
            job["state"] = parts[1] if len(parts) > 1 else None
        except NoSuchElementException:
            pass

        # Row 2 typically contains employer/organization name
        try:
            row2 = card.find_element(By.CSS_SELECTOR, "div.bti-ui-job-detail-row2")
            text = row2.text.strip()
            if text:
                job["employer"] = text
        except NoSuchElementException:
            pass

        # Right column may contain salary, job type, or dates as labeled text
        try:
            right = card.find_element(By.CSS_SELECTOR, "div.bti-ui-job-result-right-col")
            for line in right.text.splitlines():
                line = line.strip()
                if not line:
                    continue
                ll = line.lower()
                if "salary" in ll:
                    job.setdefault("salary_text", line)
                    mn, mx = _parse_salary(line)
                    job.setdefault("salary_min", mn)
                    job.setdefault("salary_max", mx)
                elif "clos" in ll or "deadline" in ll:
                    job.setdefault("date_closes", line)
                elif "post" in ll:
                    job.setdefault("date_posted", line)
                elif any(w in ll for w in ("full", "part", "contract", "temporary")):
                    job.setdefault("job_type", line)
        except NoSuchElementException:
            pass

        jobs.append(job)

    return jobs


def _scrape_cards_mid(driver) -> list[dict]:
    """Intermediate BTI format (~2021–2024): card-body bti-job-search-row with data-job attribute."""
    cards = driver.find_elements(By.CSS_SELECTOR, "div.bti-job-search-row[data-job]")
    if not cards:
        return []
    jobs = []
    for card in cards:
        job_id = card.get_attribute("data-job")
        if not job_id:
            continue
        job: dict = {"job_id": job_id}
        try:
            link = card.find_element(By.CSS_SELECTOR, "div.card-title a")
            job["title"] = link.text.strip()
            job["url"]   = link.get_attribute("href") or ""
            if job["title"]:
                job["rank_level"] = _infer_rank(job["title"])
        except NoSuchElementException:
            continue
        try:
            job["employer"] = card.find_element(By.CSS_SELECTOR, "div.card-subtitle").text.strip()
        except NoSuchElementException:
            pass
        try:
            loc = card.find_element(By.CSS_SELECTOR, "p.card-text").text.strip()
            parts = [p.strip() for p in loc.split(",")]
            job["city"]  = parts[0] if parts else None
            job["state"] = parts[1] if len(parts) > 1 else None
        except NoSuchElementException:
            pass
        jobs.append(job)
    return jobs


def _scrape_cards(driver) -> list[dict]:
    """Dispatcher: detect site era and route to the correct extractor."""
    if driver.find_elements(By.CSS_SELECTOR, "ul.js-search-regular-slot"):
        return _scrape_cards_new(driver)
    if driver.find_elements(By.CSS_SELECTOR, "div.bti-job-search-row[data-job]"):
        return _scrape_cards_mid(driver)
    if driver.find_elements(By.CSS_SELECTOR, "div.bti-ui-job-result-wrapper"):
        return _scrape_cards_old(driver)
    return []


def _enrich_detail(driver, job: dict) -> dict:
    """Visit detail page to pick up any fields missing from the list view."""
    url = job.get("url")
    if not url:
        return job
    try:
        driver.get(url)
        time.sleep(2)
    except WebDriverException:
        return job
    # Re-use the same table-row parser on the detail page
    job.update({k: v for k, v in _table_fields(driver).items() if not job.get(k)})
    return job


def scrape(base_url: str, log_fn, scrape_details: bool = False,
           save_debug: bool = False) -> list[dict]:
    all_jobs: list[dict] = []
    seen_ids: set[str] = set()   # deduplicate within this scraping session
    driver = _make_driver()
    try:
        is_wayback = "web.archive.org" in base_url
        page_num = 1
        while True:
            if is_wayback:
                # Wayback only archived the bare URL — don't append live-site params
                url = base_url if page_num == 1 else f"{base_url}?page={page_num}"
            else:
                sep = "&" if "?" in base_url else "?"
                url = f"{base_url}{sep}locale=en&page={page_num}&sort=date"
            log_fn(f"  Page {page_num} → {url}")
            try:
                driver.get(url)
                if is_wayback:
                    _card_sel = ("li.search-item[data-id], "
                                 "div.bti-job-search-row[data-job], "
                                 "div.bti-ui-job-result-wrapper")
                    # Wait up to 90s for cards — Wayback pages can be slow to
                    # load the toolbar, execute JS, and complete the API call
                    try:
                        WebDriverWait(driver, 90).until(
                            EC.presence_of_element_located((By.CSS_SELECTOR, _card_sel))
                        )
                    except TimeoutException:
                        # Cards didn't auto-load — try clicking "Find a job" to trigger search
                        try:
                            btn = driver.find_element(
                                By.XPATH,
                                "//button[contains(text(),'Find a job') or contains(text(),'Search')] | "
                                "//a[contains(text(),'Find a job')] | "
                                "//input[@type='submit']"
                            )
                            btn.click()
                            WebDriverWait(driver, 45).until(
                                EC.presence_of_element_located((By.CSS_SELECTOR, _card_sel))
                            )
                        except Exception:
                            pass
                    driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
                    time.sleep(2)
                else:
                    time.sleep(3)
            except WebDriverException as e:
                log_fn(f"  Navigation error: {e}")
                break

            if page_num == 1 and save_debug:
                DEBUG_HTML.write_text(driver.page_source, encoding="utf-8")
                log_fn(f"  Debug HTML saved → {DEBUG_HTML.name}")

            page_jobs = _scrape_cards(driver)
            # Drop any IDs already seen this session (catches Wayback looping back to page 1)
            fresh = [j for j in page_jobs if j.get("job_id") not in seen_ids]
            seen_ids.update(j["job_id"] for j in page_jobs if j.get("job_id"))
            log_fn(f"  {len(page_jobs)} listings on page {page_num}")
            if not fresh:
                break

            if scrape_details:
                for i, job in enumerate(fresh):
                    log_fn(f"    Detail {i+1}/{len(fresh)}: {job.get('title','?')[:55]}")
                    fresh[i] = _enrich_detail(driver, job)
                    time.sleep(1)

            all_jobs.extend(fresh)

            # New site: "Load More" button signals another page
            has_more = bool(driver.find_elements(By.CSS_SELECTOR, "button.js-load-more"))
            if not has_more:
                has_more = bool(driver.find_elements(By.CSS_SELECTOR, "a.bti-job-search-load-more"))
            # Old site: if we got fresh results try the next page;
            # the in-session dedup above will catch it if Wayback loops back
            if not has_more and driver.find_elements(
                By.CSS_SELECTOR, "div.bti-ui-job-result-wrapper"
            ):
                has_more = True
            if has_more and page_num < 20:
                page_num += 1
            else:
                break
    finally:
        driver.quit()

    return all_jobs


def fetch_cdx_jobs(timestamp: str, log_fn=None, window_days: int = 20) -> list[dict]:
    """Query Wayback CDX API for individual job pages archived near a timestamp."""
    ts_dt    = datetime.strptime(timestamp[:8], "%Y%m%d")
    from_str = (ts_dt - timedelta(days=window_days)).strftime("%Y%m%d")
    to_str   = (ts_dt + timedelta(days=window_days)).strftime("%Y%m%d")

    # Job IDs in the new era (2020s+) start with "2" (20XXXXXX, 21XXXXXX, etc.)
    # Old era started with "3". The broad wildcard jobs/* returns nothing;
    # prefixing with the first digit works. As of early 2025, IDs are ~20.8M
    # and may roll into 21M+ during 2025, so we use jobs/2* to catch all of them.
    cdx_url = (
        "https://web.archive.org/cdx/search/cdx"
        f"?url=careers.aapm.org/jobs/2*&output=json"
        f"&from={from_str}&to={to_str}"
        f"&fl=original&collapse=original&limit=5000"
    )
    if log_fn:
        log_fn(f"  CDX window: {from_str[:4]}-{from_str[4:6]}-{from_str[6:]} → "
               f"{to_str[:4]}-{to_str[4:6]}-{to_str[6:]}")
        log_fn(f"  CDX URL: {cdx_url}")

    data = None
    for attempt in range(3):
        try:
            req = urllib.request.Request(cdx_url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                data = json.loads(resp.read())
            break
        except Exception as e:
            if log_fn:
                log_fn(f"  CDX error (attempt {attempt+1}/3): {e}")
            if attempt < 2:
                time.sleep(5)

    if data is None:
        return []

    if log_fn:
        log_fn(f"  CDX raw rows: {len(data)}")
    if len(data) <= 1:
        return []

    jobs = []
    seen_ids: set[str] = set()
    for row in data[1:]:
        url = row[0]
        m = re.search(r"/jobs/(\d+)", url)
        if not m:
            continue
        job_id = m.group(1)
        if job_id in seen_ids:
            continue
        seen_ids.add(job_id)

        slug  = re.search(r"/jobs/\d+/([^/?#]+)", url)
        title = slug.group(1).replace("-", " ").title() if slug else None

        jobs.append({
            "job_id":     job_id,
            "url":        f"https://careers.aapm.org/jobs/{job_id}/",
            "title":      title,
            "rank_level": _infer_rank(title) if title else None,
        })

    return jobs


def _fetch_url_with_retry(url: str, timeout: int = 60, retries: int = 3,
                          log_fn=None) -> bytes | None:
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as e:
            if log_fn:
                log_fn(f"    fetch error (attempt {attempt+1}/{retries}): {e}")
            if attempt < retries - 1:
                time.sleep(5)
    return None


def _parse_sitemap_xml(content: bytes, wb_timestamp: str,
                       log_fn=None) -> list[str]:
    """Parse sitemap XML (or sitemap index) and return job-page URLs."""
    import xml.etree.ElementTree as ET
    NS = "http://www.sitemaps.org/schemas/sitemap/0.9"

    # Wayback sometimes injects an HTML toolbar before the XML — strip it
    for marker in (b"<?xml", b"<urlset", b"<sitemapindex"):
        idx = content.find(marker)
        if idx > 0:
            content = content[idx:]
            break

    try:
        root = ET.fromstring(content)
    except ET.ParseError as e:
        if log_fn:
            log_fn(f"    XML parse error: {e}")
        return []

    tag = root.tag.replace(f"{{{NS}}}", "")
    job_urls: list[str] = []

    if tag == "sitemapindex":
        # Follow each child sitemap – try the Wayback-archived version first
        for sm in root.findall(f"{{{NS}}}sitemap"):
            loc = sm.findtext(f"{{{NS}}}loc") or ""
            # Strip any Wayback prefix that Wayback injected into the XML
            live_loc = re.sub(r"https?://web\.archive\.org/web/\d+/", "", loc)
            wb_sub = f"https://web.archive.org/web/{wb_timestamp}/{live_loc}"
            if log_fn:
                log_fn(f"    sub-sitemap → {live_loc}")
            sub_content = _fetch_url_with_retry(wb_sub, log_fn=log_fn)
            if sub_content:
                job_urls.extend(_parse_sitemap_xml(sub_content, wb_timestamp, log_fn))
    else:
        for url_el in root.findall(f"{{{NS}}}url"):
            loc = url_el.findtext(f"{{{NS}}}loc") or ""
            live_loc = re.sub(r"https?://web\.archive\.org/web/\d+/", "", loc)
            if re.search(r"/jobs/\d+", live_loc):
                job_urls.append(live_loc)

    return job_urls


def fetch_sitemap_jobs(wb_timestamp: str, log_fn=None) -> list[dict]:
    """Fetch an archived careers.aapm.org sitemap and return job dicts."""
    sitemap_url = (
        f"https://web.archive.org/web/{wb_timestamp}"
        "/https://careers.aapm.org/sitemap.xml"
    )
    if log_fn:
        log_fn(f"  Sitemap: {sitemap_url}")

    content = _fetch_url_with_retry(sitemap_url, log_fn=log_fn)
    if not content:
        return []

    job_page_urls = _parse_sitemap_xml(content, wb_timestamp, log_fn)
    if log_fn:
        log_fn(f"  Sitemap job URLs found: {len(job_page_urls)}")

    jobs: list[dict] = []
    seen_ids: set[str] = set()
    for url in job_page_urls:
        m = re.search(r"/jobs/(\d+)", url)
        if not m:
            continue
        job_id = m.group(1)
        if job_id in seen_ids:
            continue
        seen_ids.add(job_id)
        slug  = re.search(r"/jobs/\d+/([^/?#]+)", url)
        title = slug.group(1).replace("-", " ").title() if slug else None
        jobs.append({
            "job_id":     job_id,
            "url":        url,
            "title":      title,
            "rank_level": _infer_rank(title) if title else None,
        })
    return jobs


# Known Wayback sitemap archives: (sitemap_timestamp, nearest_snapshot_ymd)
# Pulled from CDX: careers.aapm.org/sitemap.xml
_SITEMAP_ARCHIVES: list[tuple[str, str]] = [
    ("20250125005611", "20250125"),
    ("20250311050619", "20250311"),
    ("20250401103513", "20250401"),
    ("20250501110625", "20250501"),
]


def nearest_sitemap_ts(snapshot_ts: str) -> str | None:
    """Return the closest archived sitemap timestamp for a given snapshot."""
    snap_dt = datetime.strptime(snapshot_ts[:8], "%Y%m%d")
    best_ts, best_delta = None, timedelta(days=999)
    for sm_ts, _ in _SITEMAP_ARCHIVES:
        sm_dt = datetime.strptime(sm_ts[:8], "%Y%m%d")
        delta = abs(snap_dt - sm_dt)
        if delta < best_delta and delta <= timedelta(days=45):
            best_ts, best_delta = sm_ts, delta
    return best_ts


# ── GUI ───────────────────────────────────────────────────────────────────────
class App:
    def __init__(self, root: tk.Tk):
        self.root = root
        root.title("AAPM Jobs Tracker")
        root.resizable(True, True)
        self._busy = False
        self._build()
        self._refresh_stats()
        self._refresh_charts()
        self._refresh_trending()
        self._refresh_match_stats()
        self._check_autofetch()

    def _build(self):
        r = self.root
        r.columnconfigure(0, weight=1)
        r.rowconfigure(2, weight=1)

        # Stats
        sf = ttk.LabelFrame(r, text="Database", padding=8)
        sf.grid(row=0, column=0, sticky="ew", padx=10, pady=(10, 4))
        sf.columnconfigure((0, 1, 2, 3), weight=1)

        self.lbl_total = ttk.Label(sf, text="—", font=("Segoe UI", 14, "bold"))
        self.lbl_first = ttk.Label(sf, text="—")
        self.lbl_wb    = ttk.Label(sf, text="—")
        self.lbl_next  = ttk.Label(sf, text="—")

        for col, (val, lbl) in enumerate([
            (self.lbl_total, "Unique jobs"),
            (self.lbl_first, "Earliest data"),
            (self.lbl_wb,    "Wayback snapshots"),
            (self.lbl_next,  "Next auto-fetch"),
        ]):
            ttk.Label(sf, text=lbl, foreground="gray").grid(row=0, column=col, sticky="w")
            val.grid(row=1, column=col, sticky="w", padx=(0, 20))

        # Controls
        cf = ttk.Frame(r, padding=(10, 4))
        cf.grid(row=1, column=0, sticky="ew")

        self.btn_live = ttk.Button(cf, text="Fetch Now (Live)", command=self._fetch_live)
        self.btn_live.pack(side="left", padx=(0, 6))

        self.btn_wb = ttk.Button(cf, text="Fetch Historic (Wayback)", command=self._fetch_wayback)
        self.btn_wb.pack(side="left", padx=(0, 6))

        self.var_details = tk.BooleanVar(value=False)
        ttk.Checkbutton(cf, text="Scrape detail pages (slower, more data)",
                        variable=self.var_details).pack(side="left", padx=(10, 0))

        ttk.Button(cf, text="Retry Zero-Result Snapshots",
                   command=self._retry_zero_snapshots).pack(side="left", padx=(10, 0))

        ttk.Button(cf, text="Fill CDX Gaps",
                   command=self._fetch_cdx_gaps).pack(side="left", padx=(6, 0))

        ttk.Button(cf, text="Re-run Low Counts",
                   command=self._rerun_low_count_snapshots).pack(side="left", padx=(6, 0))

        ttk.Button(cf, text="Test 2024 Snapshot",
                   command=self._test_wayback_snapshot).pack(side="left", padx=(6, 0))

        ttk.Button(cf, text="Open Debug HTML",
                   command=self._open_debug).pack(side="right")

        # Notebook: Log + Charts
        nb = ttk.Notebook(r)
        nb.grid(row=2, column=0, sticky="nsew", padx=10, pady=(4, 4))

        # Log tab
        log_tab = ttk.Frame(nb, padding=4)
        nb.add(log_tab, text="Log")
        log_tab.rowconfigure(0, weight=1)
        log_tab.columnconfigure(0, weight=1)
        self.log = scrolledtext.ScrolledText(
            log_tab, height=22, wrap="word", font=("Consolas", 9), state="disabled"
        )
        self.log.grid(row=0, column=0, sticky="nsew")

        # Charts tab
        charts_tab = ttk.Frame(nb, padding=4)
        nb.add(charts_tab, text="Charts")
        charts_tab.rowconfigure(0, weight=1)
        charts_tab.columnconfigure(0, weight=1)

        if _HAS_CHARTS:
            self._fig = Figure(figsize=(10, 6), dpi=96, tight_layout=True)
            canvas = FigureCanvasTkAgg(self._fig, master=charts_tab)
            canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
            self._canvas = canvas
        else:
            ttk.Label(charts_tab,
                      text="Charts unavailable — install matplotlib and pandas",
                      foreground="gray").grid(row=0, column=0)

        # Trending tab
        trend_tab = ttk.Frame(nb, padding=4)
        nb.add(trend_tab, text="Trending")
        trend_tab.rowconfigure(0, weight=1)
        trend_tab.columnconfigure(0, weight=1)

        if _HAS_CHARTS:
            self._trend_fig = Figure(figsize=(10, 6), dpi=96, tight_layout=True)
            trend_canvas = FigureCanvasTkAgg(self._trend_fig, master=trend_tab)
            trend_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
            self._trend_canvas = trend_canvas
        else:
            ttk.Label(trend_tab,
                      text="Charts unavailable — install matplotlib and pandas",
                      foreground="gray").grid(row=0, column=0)

        # Match Stats tab
        ms_tab = ttk.Frame(nb, padding=4)
        nb.add(ms_tab, text="Match Stats")
        ms_tab.rowconfigure(0, weight=2)
        ms_tab.rowconfigure(1, weight=1)
        ms_tab.columnconfigure(0, weight=1)

        if _HAS_CHARTS:
            self._ms_fig = Figure(figsize=(10, 4), dpi=96, tight_layout=True)
            ms_canvas = FigureCanvasTkAgg(self._ms_fig, master=ms_tab)
            ms_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
            self._ms_canvas = ms_canvas
        else:
            ttk.Label(ms_tab, text="Charts unavailable",
                      foreground="gray").grid(row=0, column=0)

        ms_lower = ttk.Frame(ms_tab)
        ms_lower.grid(row=1, column=0, sticky="nsew", pady=(4, 0))
        ms_lower.rowconfigure(1, weight=1)
        ms_lower.columnconfigure(0, weight=1)

        # Controls row
        ms_ctrl = ttk.Frame(ms_lower)
        ms_ctrl.grid(row=0, column=0, sticky="ew", pady=(0, 4))
        ttk.Label(ms_ctrl, text="Year:").pack(side="left")
        self._ms_year = ttk.Entry(ms_ctrl, width=6); self._ms_year.pack(side="left", padx=4)
        ttk.Label(ms_ctrl, text="Registered:").pack(side="left")
        self._ms_reg  = ttk.Entry(ms_ctrl, width=6); self._ms_reg.pack(side="left", padx=4)
        ttk.Label(ms_ctrl, text="Offered:").pack(side="left")
        self._ms_off  = ttk.Entry(ms_ctrl, width=6); self._ms_off.pack(side="left", padx=4)
        ttk.Label(ms_ctrl, text="Matched:").pack(side="left")
        self._ms_mat  = ttk.Entry(ms_ctrl, width=6); self._ms_mat.pack(side="left", padx=4)
        ttk.Label(ms_ctrl, text="Unmatched:").pack(side="left")
        self._ms_unm  = ttk.Entry(ms_ctrl, width=6); self._ms_unm.pack(side="left", padx=4)
        ttk.Button(ms_ctrl, text="Add / Update",
                   command=self._ms_upsert).pack(side="left", padx=(8, 4))
        ttk.Button(ms_ctrl, text="Delete Selected",
                   command=self._ms_delete).pack(side="left")

        # Treeview
        ms_cols = ("year","registered_applicants","positions_offered",
                   "matches_filled","unmatched_applicants")
        self._ms_tree = ttk.Treeview(ms_lower, columns=ms_cols, show="headings", height=6)
        for col, hdr, w in zip(ms_cols,
                                ("Year","Registered","Offered","Matched","Unmatched"),
                                (60, 100, 90, 90, 110)):
            self._ms_tree.heading(col, text=hdr)
            self._ms_tree.column(col, width=w, anchor="center")
        self._ms_tree.grid(row=1, column=0, sticky="nsew")
        self._ms_tree.bind("<<TreeviewSelect>>", self._ms_on_select)
        sb = ttk.Scrollbar(ms_lower, orient="vertical", command=self._ms_tree.yview)
        sb.grid(row=1, column=1, sticky="ns")
        self._ms_tree.configure(yscrollcommand=sb.set)

        # Database tab
        db_tab = ttk.Frame(nb, padding=4)
        nb.add(db_tab, text="Database")
        db_tab.rowconfigure(1, weight=1)
        db_tab.columnconfigure(0, weight=1)

        db_ctrl = ttk.Frame(db_tab)
        db_ctrl.grid(row=0, column=0, sticky="ew", pady=(0, 4))

        self.var_table = tk.StringVar(value="fetch_log")
        for tbl in ("fetch_log", "jobs", "settings"):
            ttk.Radiobutton(db_ctrl, text=tbl, variable=self.var_table,
                            value=tbl, command=self._load_table).pack(side="left", padx=(0, 6))
        ttk.Button(db_ctrl, text="Refresh", command=self._load_table).pack(side="left", padx=(12, 0))
        ttk.Button(db_ctrl, text="Delete Selected", command=self._delete_selected).pack(side="left", padx=(6, 0))

        tree_frame = ttk.Frame(db_tab)
        tree_frame.grid(row=1, column=0, sticky="nsew")
        tree_frame.rowconfigure(0, weight=1)
        tree_frame.columnconfigure(0, weight=1)
        self._tree = ttk.Treeview(tree_frame, selectmode="extended", show="headings")
        vsb = ttk.Scrollbar(tree_frame, orient="vertical",   command=self._tree.yview)
        hsb = ttk.Scrollbar(tree_frame, orient="horizontal", command=self._tree.xview)
        self._tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        self._tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        self._tree_cols: list[str] = []

        sql_frame = ttk.LabelFrame(db_tab, text="Run SQL", padding=4)
        sql_frame.grid(row=2, column=0, sticky="ew", pady=(4, 0))
        sql_frame.columnconfigure(0, weight=1)
        self.sql_entry = tk.Text(sql_frame, height=2, font=("Consolas", 9))
        self.sql_entry.grid(row=0, column=0, sticky="ew")
        ttk.Button(sql_frame, text="Run", command=self._run_sql).grid(row=0, column=1, padx=(4, 0), sticky="ns")

        self.progress = ttk.Progressbar(r, mode="indeterminate")
        self.progress.grid(row=3, column=0, sticky="ew", padx=10, pady=(0, 8))

    def _log(self, msg: str):
        def _append():
            self.log.configure(state="normal")
            self.log.insert("end", f"{datetime.now().strftime('%H:%M:%S')}  {msg}\n")
            self.log.see("end")
            self.log.configure(state="disabled")
        self.root.after(0, _append)

    def _refresh_stats(self):
        s = get_stats()
        self.lbl_total.configure(text=str(s["total"]))
        self.lbl_first.configure(text=(s["first_seen"] or "—")[:10])
        self.lbl_wb.configure(text=f"{s['wb_done']} / {len(WAYBACK_200)}")
        last = get_last_fetch()
        if last:
            nxt = last + timedelta(days=FETCH_DAYS)
            self.lbl_next.configure(text=nxt.strftime("%Y-%m-%d"))
        else:
            self.lbl_next.configure(text="Now")

    def _refresh_charts(self):
        if not _HAS_CHARTS:
            return
        with _conn() as conn:
            df = pd.read_sql("SELECT * FROM jobs", conn)
        if df.empty:
            return

        C  = "#2B5DA8"   # primary blue
        C2 = "#6B9DD4"   # lighter blue
        BG = "#F5F6F8"
        GRID = "#D8DCE8"

        self._fig.clear()
        self._fig.patch.set_facecolor(BG)
        axes = self._fig.subplots(2, 2)

        def _style(ax, title):
            ax.set_facecolor(BG)
            ax.set_title(title, fontsize=9, fontweight="bold", pad=6)
            ax.tick_params(labelsize=8)
            ax.spines[["top", "right"]].set_visible(False)
            ax.spines[["left", "bottom"]].set_color(GRID)
            ax.xaxis.grid(True, color=GRID, linewidth=0.6)
            ax.set_axisbelow(True)

        # ── Rank breakdown ────────────────────────────────────────────────
        ax = axes[0, 0]
        ranks = df["rank_level"].fillna("Unspecified")
        cnt = ranks.value_counts().sort_values()
        bars = ax.barh(cnt.index, cnt.values, color=C, height=0.6)
        for bar, v in zip(bars, cnt.values):
            ax.text(v + 0.3, bar.get_y() + bar.get_height() / 2,
                    str(v), va="center", fontsize=7, color="#1A2030")
        _style(ax, "Rank / Position")
        ax.set_xlabel("Listings", fontsize=8)

        # ── Top states ───────────────────────────────────────────────────
        ax = axes[0, 1]
        states = df["state"].dropna()
        cnt = states.value_counts().head(12).sort_values()
        bars = ax.barh(cnt.index, cnt.values, color=C2, height=0.6)
        for bar, v in zip(bars, cnt.values):
            ax.text(v + 0.3, bar.get_y() + bar.get_height() / 2,
                    str(v), va="center", fontsize=7, color="#1A2030")
        _style(ax, "Top States")
        ax.set_xlabel("Listings", fontsize=8)

        # ── Job type ─────────────────────────────────────────────────────
        ax = axes[1, 0]
        jt = df["job_type"].fillna("Not listed")
        cnt = jt.value_counts().sort_values()
        bars = ax.barh(cnt.index, cnt.values, color=C, height=0.6)
        for bar, v in zip(bars, cnt.values):
            ax.text(v + 0.3, bar.get_y() + bar.get_height() / 2,
                    str(v), va="center", fontsize=7, color="#1A2030")
        _style(ax, "Job Type")
        ax.set_xlabel("Listings", fontsize=8)

        # ── Salary coverage + distribution ────────────────────────────────
        ax = axes[1, 1]
        has_sal = df["salary_min"].notna()
        n_sal = has_sal.sum()
        pct = 100 * n_sal / max(len(df), 1)
        if n_sal >= 3:
            sal_df = df[has_sal].copy()
            sal_df["rank_label"] = sal_df["rank_level"].fillna("Unspecified")
            grp = sal_df.groupby("rank_label")["salary_min"]
            labels = [k for k, _ in grp]
            data   = [v.values / 1000 for _, v in grp]
            ax.boxplot(data, labels=labels, vert=False, patch_artist=True,
                            medianprops=dict(color="#1A2030", linewidth=1.5),
                            boxprops=dict(facecolor=C2, linewidth=0.8),
                            whiskerprops=dict(linewidth=0.8),
                            capprops=dict(linewidth=0.8),
                            flierprops=dict(marker="o", markersize=3,
                                            markerfacecolor=C, linewidth=0))
            _style(ax, f"Salary by Rank  ({n_sal}/{len(df)} listed, {pct:.0f}%)")
            ax.set_xlabel("Salary Min ($k)", fontsize=8)
        else:
            ax.text(0.5, 0.5,
                    f"Salary listed on {n_sal} of {len(df)} jobs ({pct:.0f}%)\n"
                    "(not enough data for distribution)",
                    ha="center", va="center", fontsize=9, color="#4A5568",
                    transform=ax.transAxes)
            _style(ax, "Salary Coverage")
            ax.set_xticks([])
            ax.set_yticks([])

        self._canvas.draw()

    def _refresh_trending(self):
        if not _HAS_CHARTS:
            return
        with _conn() as conn:
            log_df = pd.read_sql(
                "SELECT fetched_at, source, jobs_found, new_jobs FROM fetch_log ORDER BY fetched_at",
                conn,
            )
        if log_df.empty:
            return

        def _plot_date(row):
            if row["source"] == "live":
                return pd.to_datetime(row["fetched_at"])
            return pd.to_datetime(str(row["source"])[:8], format="%Y%m%d")

        log_df["date"] = log_df.apply(_plot_date, axis=1)
        log_df = log_df.sort_values("date").reset_index(drop=True)

        C    = "#2B5DA8"
        C2   = "#E07B39"
        BG   = "#F5F6F8"
        GRID = "#D8DCE8"

        self._trend_fig.clear()
        self._trend_fig.patch.set_facecolor(BG)
        ax = self._trend_fig.add_subplot(111)

        import matplotlib.dates as mdates

        ax.plot(log_df["date"], log_df["jobs_found"], color=C, linewidth=2,
                marker="o", markersize=5, zorder=3, label="Active listings")
        ax.fill_between(log_df["date"], log_df["jobs_found"], alpha=0.10, color=C)

        ax.plot(log_df["date"], log_df["new_jobs"], color=C2, linewidth=1.5,
                marker="o", markersize=4, zorder=3, label="New listings")
        ax.fill_between(log_df["date"], log_df["new_jobs"], alpha=0.10, color=C2)

        ax.set_facecolor(BG)
        ax.set_title("Job Postings Over Time", fontsize=11, fontweight="bold", pad=8)
        ax.set_ylabel("Listings", fontsize=9)
        ax.tick_params(labelsize=8, axis="x", rotation=30)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color(GRID)
        ax.yaxis.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.legend(fontsize=8, framealpha=0.7)

        self._trend_canvas.draw()

    def _set_busy(self, busy: bool):
        self._busy = busy
        state = "disabled" if busy else "normal"
        self.btn_live.configure(state=state)
        self.btn_wb.configure(state=state)
        if busy:
            self.progress.start(12)
        else:
            self.progress.stop()
            self._refresh_stats()
            self._refresh_charts()
            self._refresh_trending()

    def _run_in_thread(self, fn):
        def target():
            try:
                fn()
            except Exception as e:
                self._log(f"ERROR: {e}")
            finally:
                self.root.after(0, lambda: self._set_busy(False))
        self._set_busy(True)
        threading.Thread(target=target, daemon=True).start()

    def _fetch_live(self):
        if self._busy:
            return
        details = self.var_details.get()
        self._log("=== Live fetch started ===")

        def task():
            jobs = scrape(LIVE_BASE, self._log, scrape_details=details, save_debug=True)
            if jobs:
                found, new = upsert_jobs(jobs, "live")
                self._log(f"Done — {found} listings found, {new} new to DB.")
            else:
                self._log("No jobs extracted. Open Debug HTML to inspect page structure.")

        self._run_in_thread(task)

    def _fetch_wayback(self):
        if self._busy:
            return
        done      = get_wayback_done()
        remaining = [ts for ts in WAYBACK_200 if ts not in done]
        if not remaining:
            messagebox.showinfo("Wayback", "All historic snapshots already fetched.")
            return
        if not messagebox.askyesno(
            "Fetch Historic",
            f"{len(remaining)} snapshots remaining ({WAYBACK_200[0][:8]} – present).\n"
            "This may take several hours. Proceed?"
        ):
            return
        self._log(f"=== Wayback fetch: {len(remaining)} snapshots ===")

        def task():
            for i, ts in enumerate(remaining):
                dt_str = f"{ts[:4]}-{ts[4:6]}-{ts[6:8]}"
                self._log(f"[{i+1}/{len(remaining)}] {dt_str}")
                url  = WB_BASE.format(ts=ts)
                jobs = scrape(url, self._log, scrape_details=False, save_debug=False)
                found, new = upsert_jobs(jobs, ts)
                if found:
                    self._log(f"  → {found} listings, {new} new.")
                else:
                    self._log("  → No listings extracted (logged as done).")
                time.sleep(4)

        self._run_in_thread(task)

    _TABLE_PK = {"jobs": "job_id", "fetch_log": "id", "settings": "key"}

    def _load_table(self):
        table = self.var_table.get()
        with _conn() as conn:
            cur = conn.execute(f"SELECT * FROM {table} LIMIT 1000")
            cols = [d[0] for d in cur.description]
            rows = cur.fetchall()
        self._tree_cols = cols
        self._tree.delete(*self._tree.get_children())
        self._tree["columns"] = cols
        for col in cols:
            w = 180 if col in ("title", "url", "employer", "fetched_at", "first_seen", "last_seen") else 80
            self._tree.heading(col, text=col)
            self._tree.column(col, width=w, minwidth=40, stretch=False)
        for row in rows:
            self._tree.insert("", "end", values=list(row))

    def _delete_selected(self):
        selected = self._tree.selection()
        if not selected:
            messagebox.showinfo("Delete", "No rows selected.")
            return
        table = self.var_table.get()
        pk_col = self._TABLE_PK.get(table)
        if not pk_col or pk_col not in self._tree_cols:
            messagebox.showwarning("Delete", "Cannot determine primary key for this view.")
            return
        pk_idx = self._tree_cols.index(pk_col)
        pk_vals = [self._tree.item(iid)["values"][pk_idx] for iid in selected]
        if not messagebox.askyesno("Delete", f"Delete {len(pk_vals)} row(s) from {table}?"):
            return
        with _conn() as conn:
            for pk in pk_vals:
                conn.execute(f"DELETE FROM {table} WHERE {pk_col}=?", (pk,))
        for iid in selected:
            self._tree.delete(iid)
        self._refresh_stats()
        self._log(f"Deleted {len(pk_vals)} row(s) from {table}.")

    def _run_sql(self):
        sql = self.sql_entry.get("1.0", "end").strip()
        if not sql:
            return
        try:
            with _conn() as conn:
                cur = conn.execute(sql)
                if cur.description:
                    cols = [d[0] for d in cur.description]
                    rows = cur.fetchall()
                    self._tree_cols = cols
                    self._tree.delete(*self._tree.get_children())
                    self._tree["columns"] = cols
                    for col in cols:
                        self._tree.heading(col, text=col)
                        self._tree.column(col, width=120, minwidth=40, stretch=False)
                    for row in rows:
                        self._tree.insert("", "end", values=list(row))
                    self._log(f"SQL OK — {len(rows)} rows.")
                else:
                    self._log(f"SQL OK — {cur.rowcount} row(s) affected.")
                    self._refresh_stats()
        except Exception as e:
            self._log(f"SQL error: {e}")

    def _fetch_cdx_gaps(self):
        if self._busy:
            return
        with _conn() as conn:
            zero_ts = [r["source"] for r in conn.execute(
                "SELECT source FROM fetch_log WHERE source != 'live' AND jobs_found = 0"
            ).fetchall()]
        if not zero_ts:
            messagebox.showinfo("CDX Fill", "No zero-result snapshots to fill.")
            return
        if not messagebox.askyesno(
            "Fill CDX Gaps",
            f"Query Wayback CDX API for {len(zero_ts)} zero-result snapshot(s).\n"
            f"Uses a ±20 day window around each date to find archived job pages.\nProceed?"
        ):
            return
        self._log(f"=== CDX gap fill: {len(zero_ts)} snapshots ===")

        def task():
            for i, ts in enumerate(sorted(zero_ts)):
                dt_str = f"{ts[:4]}-{ts[4:6]}-{ts[6:8]}"
                self._log(f"[{i+1}/{len(zero_ts)}] CDX {dt_str}")

                jobs = fetch_cdx_jobs(ts, self._log, window_days=20)
                if not jobs:
                    self._log("  → No results in ±20d window, expanding to ±60d...")
                    jobs = fetch_cdx_jobs(ts, self._log, window_days=60)

                # Always supplement with sitemap archive when one is nearby
                sm_ts = nearest_sitemap_ts(ts)
                if sm_ts:
                    self._log(f"  → Sitemap archive {sm_ts[:8]}...")
                    sm_jobs = fetch_sitemap_jobs(sm_ts, self._log)
                    if sm_jobs:
                        existing_ids = {j["job_id"] for j in jobs}
                        extra = [j for j in sm_jobs if j["job_id"] not in existing_ids]
                        jobs = jobs + extra
                        self._log(f"  → Sitemap added {len(extra)} more (CDX+sitemap total: {len(jobs)})")

                if jobs:
                    with _conn() as conn:
                        conn.execute("DELETE FROM fetch_log WHERE source=?", (ts,))
                    found, new = upsert_jobs(jobs, ts)
                    self._log(f"  → {found} jobs found, {new} new to DB.")
                else:
                    self._log("  → No data found.")
                time.sleep(2)

        self._run_in_thread(task)

    def _test_wayback_snapshot(self):
        if self._busy:
            return
        url = "https://web.archive.org/web/20250213182247/https://careers.aapm.org/jobs/"
        self._log(f"Testing: {url}")

        def task():
            jobs = scrape(url, self._log, scrape_details=False, save_debug=True)
            self._log(f"Test result: {len(jobs)} jobs found.")

        self._run_in_thread(task)

    def _retry_zero_snapshots(self):
        with _conn() as conn:
            zero_count = conn.execute(
                "SELECT COUNT(*) FROM fetch_log WHERE source != 'live' AND jobs_found = 0"
            ).fetchone()[0]
            all_count = conn.execute(
                "SELECT COUNT(*) FROM fetch_log WHERE source != 'live'"
            ).fetchone()[0]

        if all_count == 0:
            messagebox.showinfo("Retry", "No Wayback snapshots have been logged yet.")
            return

        choice = messagebox.askyesnocancel(
            "Reset Wayback Log",
            f"Zero-result snapshots: {zero_count}\n"
            f"All Wayback snapshots logged: {all_count}\n\n"
            f"Yes  → Clear zero-result only ({zero_count})\n"
            f"No   → Clear ALL Wayback entries ({all_count})\n"
            f"Cancel → Do nothing"
        )
        if choice is None:
            return
        with _conn() as conn:
            if choice:
                conn.execute("DELETE FROM fetch_log WHERE source != 'live' AND jobs_found = 0")
                self._log(f"Cleared {zero_count} zero-result Wayback snapshot(s).")
            else:
                conn.execute("DELETE FROM fetch_log WHERE source != 'live'")
                self._log(f"Cleared all {all_count} Wayback log entries — full re-fetch on next run.")
        self._refresh_stats()

    def _rerun_low_count_snapshots(self):
        """Reset non-zero but suspiciously-low Wayback entries so Fill CDX Gaps re-processes them."""
        THRESHOLD = 50
        with _conn() as conn:
            rows = conn.execute(
                "SELECT source, jobs_found FROM fetch_log "
                "WHERE source != 'live' AND jobs_found > 0 AND jobs_found < ?"
                " AND source >= '20250101'",
                (THRESHOLD,)
            ).fetchall()

        if not rows:
            messagebox.showinfo("Re-run Low Counts",
                                f"No Wayback snapshots with 1–{THRESHOLD-1} jobs found.")
            return

        detail = "\n".join(
            f"  {r[0][:4]}-{r[0][4:6]}-{r[0][6:8]}  ({r[1]} jobs)"
            for r in sorted(rows, key=lambda r: r[0])
        )
        if not messagebox.askyesno(
            "Re-run Low Counts",
            f"Reset {len(rows)} low-count Wayback snapshot(s) to zero so\n"
            f"Fill CDX Gaps re-processes them:\n\n{detail}\n\n"
            f"Proceed?"
        ):
            return

        with _conn() as conn:
            for ts, _ in rows:
                conn.execute("UPDATE fetch_log SET jobs_found = 0 WHERE source = ?", (ts,))
        self._log(f"Reset {len(rows)} low-count snapshot(s) to zero — run Fill CDX Gaps next.")
        self._refresh_stats()

    # ── Match Stats ───────────────────────────────────────────────────────────

    def _refresh_match_stats(self):
        with _conn() as conn:
            rows = conn.execute(
                "SELECT year,registered_applicants,positions_offered,"
                "matches_filled,unmatched_applicants FROM match_stats ORDER BY year"
            ).fetchall()

        self._ms_tree.delete(*self._ms_tree.get_children())
        for row in rows:
            self._ms_tree.insert("", "end", values=list(row))

        if not _HAS_CHARTS or not rows:
            return

        years = [r[0] for r in rows]
        reg   = [r[1] for r in rows]
        off   = [r[2] for r in rows]
        mat   = [r[3] for r in rows]
        unm   = [r[4] for r in rows]

        self._ms_fig.clear()
        ax = self._ms_fig.add_subplot(111)
        ax.plot(years, reg, "o-", color="#d62728", label="Registered Applicants")
        ax.plot(years, mat, "o-", color="#1f77b4", label="Matches / Positions Filled")
        ax.plot(years, unm, "o-", color="#ff7f0e", label="Unmatched Applicants")
        ax.plot(years, off, "o-", color="#2ca02c", label="Positions Offered")
        ax.set_title("MedPhys Match Statistics")
        ax.set_xlabel("Year")
        ax.set_ylabel("Count")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        self._ms_canvas.draw()

    def _ms_on_select(self, _=None):
        sel = self._ms_tree.selection()
        if not sel:
            return
        vals = self._ms_tree.item(sel[0])["values"]
        for entry, val in zip(
            (self._ms_year, self._ms_reg, self._ms_off, self._ms_mat, self._ms_unm), vals
        ):
            entry.delete(0, "end")
            entry.insert(0, str(val))

    def _ms_upsert(self):
        try:
            year = int(self._ms_year.get())
            reg  = int(self._ms_reg.get())
            off  = int(self._ms_off.get())
            mat  = int(self._ms_mat.get())
            unm  = int(self._ms_unm.get())
        except ValueError:
            messagebox.showerror("Match Stats", "All fields must be integers.")
            return
        with _conn() as conn:
            conn.execute(
                "INSERT OR REPLACE INTO match_stats "
                "(year,registered_applicants,positions_offered,matches_filled,unmatched_applicants)"
                " VALUES (?,?,?,?,?)",
                (year, reg, off, mat, unm)
            )
        self._refresh_match_stats()

    def _ms_delete(self):
        sel = self._ms_tree.selection()
        if not sel:
            return
        year = self._ms_tree.item(sel[0])["values"][0]
        if not messagebox.askyesno("Delete", f"Delete match stats for {year}?"):
            return
        with _conn() as conn:
            conn.execute("DELETE FROM match_stats WHERE year=?", (year,))
        self._refresh_match_stats()

    def _check_autofetch(self):
        last = get_last_fetch()
        if last and datetime.now() - last >= timedelta(days=FETCH_DAYS):
            self._log(f"Auto-fetch triggered (last: {last.strftime('%Y-%m-%d')})")
            self.root.after(2000, self._fetch_live)

    def _open_debug(self):
        if DEBUG_HTML.exists():
            subprocess.Popen(["start", "", str(DEBUG_HTML)], shell=True)
        else:
            messagebox.showinfo("Debug HTML", "No debug HTML yet — run a Live fetch first.")


# ── Entry point ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
    init_db()
    root = tk.Tk()
    App(root)
    root.mainloop()
