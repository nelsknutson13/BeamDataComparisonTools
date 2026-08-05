#!/usr/bin/env python3
from __future__ import annotations
"""
AAPM Jobs Tracker
Scrapes careers.aapm.org/jobs (live + Wayback Machine) and builds a local
SQLite database for medical physics workforce demand analysis.

Dependencies (install once):
    conda run pip install selenium webdriver-manager
"""

import re
import sqlite3
import subprocess
import threading
import time
import tkinter as tk
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
from selenium.common.exceptions import NoSuchElementException, WebDriverException

# ── Constants ─────────────────────────────────────────────────────────────────
HERE       = Path(__file__).parent
DB_PATH    = HERE / "aapm_jobs.db"
DEBUG_HTML = HERE / "aapm_debug_page.html"
FETCH_DAYS = 14
LIVE_BASE  = "https://careers.aapm.org/jobs"
WB_BASE    = "https://web.archive.org/web/{ts}/https://careers.aapm.org/jobs"

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
    "20240501155132", "20240814104916", "20241107085239", "20250213182247",
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
        """)


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


def _scrape_cards(driver) -> list[dict]:
    # Use only the regular (non-sponsored) slot — sponsored listings appear on every page
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

        # Title + URL
        try:
            link = card.find_element(By.CSS_SELECTOR, "a.js-item-job-url")
            job["title"] = link.text.strip()
            job["url"]   = link.get_attribute("href") or ""
            if job["title"]:
                job["rank_level"] = _infer_rank(job["title"])
        except NoSuchElementException:
            continue

        # Employer (first div.small) + Location (second div.small) in div.location
        loc_divs = card.find_elements(By.CSS_SELECTOR, "div.location div.small")
        if len(loc_divs) >= 1:
            job["employer"] = loc_divs[0].text.strip()
        if len(loc_divs) >= 2:
            loc = loc_divs[1].text.strip()
            parts = [p.strip() for p in loc.split(",")]
            job["city"]  = parts[0]
            job["state"] = parts[1] if len(parts) > 1 else None

        # Labeled table rows: salary, position type, dates
        job.update(_table_fields(card))

        jobs.append(job)

    return jobs


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
    driver = _make_driver()
    try:
        page_num = 1
        while True:
            sep = "&" if "?" in base_url else "?"
            url = f"{base_url}{sep}locale=en&page={page_num}&sort=date"
            log_fn(f"  Page {page_num} → {url}")
            try:
                driver.get(url)
                time.sleep(3)
            except WebDriverException as e:
                log_fn(f"  Navigation error: {e}")
                break

            if page_num == 1 and save_debug:
                DEBUG_HTML.write_text(driver.page_source, encoding="utf-8")
                log_fn(f"  Debug HTML saved → {DEBUG_HTML.name}")

            page_jobs = _scrape_cards(driver)
            log_fn(f"  {len(page_jobs)} listings on page {page_num}")
            if not page_jobs:
                break

            if scrape_details:
                for i, job in enumerate(page_jobs):
                    log_fn(f"    Detail {i+1}/{len(page_jobs)}: {job.get('title','?')[:55]}")
                    page_jobs[i] = _enrich_detail(driver, job)
                    time.sleep(1)

            all_jobs.extend(page_jobs)

            # Site uses a "Load More" button; presence means another page exists
            if driver.find_elements(By.CSS_SELECTOR, "button.js-load-more"):
                page_num += 1
            else:
                break
    finally:
        driver.quit()

    return all_jobs


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
                if jobs:
                    found, new = upsert_jobs(jobs, ts)
                    self._log(f"  → {found} listings, {new} new.")
                else:
                    self._log("  → No listings extracted.")
                time.sleep(4)

        self._run_in_thread(task)

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
