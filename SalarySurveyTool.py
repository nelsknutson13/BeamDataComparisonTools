import tkinter as tk
from tkinter import ttk, messagebox
import datetime

# ---------------------------------------------------------------------------
# Constants — match AAPM survey structure exactly
# ---------------------------------------------------------------------------
EXP_BINS    = ["0-2", "3-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30+"]
COLUMNS     = ["Average", "20th", "Median", "80th"]
DEGREE_CATS = ["PhD - With Cert", "Masters - With Cert", "PhD - No Cert", "Masters - No Cert"]
ROLES       = ["Physicist", "Site Lead", "Service Chief", "Director"]
ADMIN_PCT   = {"Physicist": 0, "Site Lead": 3, "Service Chief": 7, "Director": 10}

# 2024 pre-loaded data — Primary Income (full dollars)
# Source: AAPM Professional Survey Report 2024
BUILTIN_2024 = {
    "PhD - With Cert": {
        "3-4":   {"Average": 198200, "20th": 180000, "Median": 200000, "80th": 215900},
        "5-9":   {"Average": 217100, "20th": 189200, "Median": 215000, "80th": 244800},
        "10-14": {"Average": 248200, "20th": 210000, "Median": 240000, "80th": 290000},
        "15-19": {"Average": 270400, "20th": 233000, "Median": 260000, "80th": 319000},
        "20-24": {"Average": 275400, "20th": 231600, "Median": 269800, "80th": 330000},
        "25-29": {"Average": 285900, "20th": 220000, "Median": 284000, "80th": 350000},
        "30+":   {"Average": 278700, "20th": 200000, "Median": 283500, "80th": 359700},
    },
    "Masters - With Cert": {
        "5-9":   {"Average": 210700, "20th": 183900, "Median": 210000, "80th": 238800},
        "10-14": {"Average": 230400, "20th": 200000, "Median": 228000, "80th": 255200},
        "15-19": {"Average": 254700, "20th": 213300, "Median": 250000, "80th": 291000},
        "20-24": {"Average": 264300, "20th": 240000, "Median": 261700, "80th": 293300},
        "25-29": {"Average": 263900, "20th": 225000, "Median": 258500, "80th": 300000},
        "30+":   {"Average": 254100, "20th": 220000, "Median": 258600, "80th": 300000},
    },
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def parse_dollar(s):
    return float(str(s).replace("$", "").replace(",", "").strip())

def fmt_dollar(v):
    return f"${v:,.0f}"

def exp_to_bin(years):
    if years <= 2:  return "0-2"
    if years <= 4:  return "3-4"
    if years <= 9:  return "5-9"
    if years <= 14: return "10-14"
    if years <= 19: return "15-19"
    if years <= 24: return "20-24"
    if years <= 29: return "25-29"
    return "30+"

def estimate_percentile(salary, row):
    """Return estimated percentile string given salary and survey row dict."""
    pts = []
    for lbl, pct in [("20th", 20), ("Median", 50), ("80th", 80)]:
        v = row.get(lbl)
        if v is not None:
            try:
                pts.append((pct, float(v)))
            except (TypeError, ValueError):
                pass
    if not pts:
        return "n/a"
    pts.sort(key=lambda x: x[1])
    if salary < pts[0][1]:
        return f"< {pts[0][0]}th"
    if salary > pts[-1][1]:
        return f"> {pts[-1][0]}th"
    for i in range(len(pts) - 1):
        p1, v1 = pts[i]
        p2, v2 = pts[i + 1]
        if v1 <= salary <= v2:
            frac = (salary - v1) / (v2 - v1) if v2 != v1 else 0.5
            est = int(round(p1 + frac * (p2 - p1)))
            suf = "st" if est % 10 == 1 and est != 11 else (
                  "nd" if est % 10 == 2 and est != 12 else (
                  "rd" if est % 10 == 3 and est != 13 else "th"))
            return f"~{est}{suf}"
    return "n/a"

def bar_line(salary, row, label, width=38):
    """Single ASCII bar line: label: |──┼──●──|  ~Nth"""
    try:
        v20  = float(row["20th"])
        vmed = float(row["Median"])
        v80  = float(row["80th"])
    except (KeyError, TypeError, ValueError):
        return ""
    span = v80 - v20
    if span <= 0:
        return ""
    bar = ["─"] * width
    med_p = max(0, min(width - 1, int((vmed - v20) / span * (width - 1))))
    bar[med_p] = "┼"
    pct_str = estimate_percentile(salary, row)
    if salary < v20:
        return f"  {label:<14} ◄ |{''.join(bar)}|  {pct_str}"
    if salary > v80:
        return f"  {label:<14}   |{''.join(bar)}| ►  {pct_str}"
    sal_p = max(0, min(width - 1, int((salary - v20) / span * (width - 1))))
    bar[sal_p] = "◎" if sal_p == med_p else "●"
    return f"  {label:<14}   |{''.join(bar)}|  {pct_str}"


# ---------------------------------------------------------------------------
# Application
# ---------------------------------------------------------------------------
class SalarySurveyTool:
    def __init__(self, root):
        self.root = root
        self.root.title("AAPM Salary Survey Tool")
        self.root.geometry("1060x800")
        self.root.resizable(True, True)

        # survey_data[year][degree][exp_bin][col] = tk.StringVar
        self.survey_data = {}
        self.years = []

        self._add_year("2024", BUILTIN_2024)
        self._build_ui()

    # ── Data management ────────────────────────────────────────────────────────
    def _add_year(self, year, prefill=None):
        if year in self.survey_data:
            return
        self.years.append(year)
        self.years.sort()
        self.survey_data[year] = {}
        for deg in DEGREE_CATS:
            self.survey_data[year][deg] = {}
            for exp in EXP_BINS:
                self.survey_data[year][deg][exp] = {}
                for col in COLUMNS:
                    raw = ""
                    if prefill and deg in prefill:
                        v = prefill[deg].get(exp, {}).get(col, "")
                        if v:
                            raw = f"{v:,.0f}" if isinstance(v, (int, float)) else str(v)
                    self.survey_data[year][deg][exp][col] = tk.StringVar(value=raw)

    def _remove_year(self, year):
        if year in self.survey_data:
            del self.survey_data[year]
            self.years.remove(year)

    # ── UI ─────────────────────────────────────────────────────────────────────
    def _build_ui(self):
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except Exception:
            pass
        style.configure("H.TLabel", font=("Helvetica", 10, "bold"))

        nb = ttk.Notebook(self.root)
        nb.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

        self.tab_calc = ttk.Frame(nb)
        self.tab_data = ttk.Frame(nb)
        nb.add(self.tab_calc, text="  Salary Calculator  ")
        nb.add(self.tab_data, text="  Survey Data  ")

        self._build_calc_tab()
        self._build_data_tab()

    # ── Survey Data tab ────────────────────────────────────────────────────────
    def _build_data_tab(self):
        ttk.Label(self.tab_data,
            text="Primary Income percentiles ($/yr) by years experience.  "
                 "2024 PhD data pre-loaded.  Add prior years to compare.",
            foreground="gray").pack(pady=(8, 2), padx=10, anchor=tk.W)

        ctrl = ttk.Frame(self.tab_data)
        ctrl.pack(fill=tk.X, padx=10, pady=6)

        ttk.Label(ctrl, text="Year:").pack(side=tk.LEFT)
        self.d_year = tk.StringVar(value="2024")
        self.year_cb = ttk.Combobox(ctrl, textvariable=self.d_year,
                                    values=self.years, state="readonly", width=8)
        self.year_cb.pack(side=tk.LEFT, padx=4)
        self.year_cb.bind("<<ComboboxSelected>>", lambda _: self._refresh_table())

        self.new_yr = tk.StringVar()
        ttk.Entry(ctrl, textvariable=self.new_yr, width=7).pack(side=tk.LEFT, padx=(16, 2))
        ttk.Button(ctrl, text="Add Year",    command=self._add_year_ui).pack(side=tk.LEFT, padx=2)
        ttk.Button(ctrl, text="Remove Year", command=self._remove_year_ui).pack(side=tk.LEFT, padx=2)

        ttk.Separator(ctrl, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=12, pady=2)

        ttk.Label(ctrl, text="Degree:").pack(side=tk.LEFT)
        self.d_deg = tk.StringVar(value=DEGREE_CATS[0])
        deg_cb = ttk.Combobox(ctrl, textvariable=self.d_deg,
                               values=DEGREE_CATS, state="readonly", width=24)
        deg_cb.pack(side=tk.LEFT, padx=4)
        deg_cb.bind("<<ComboboxSelected>>", lambda _: self._refresh_table())

        # Scrollable table
        cont = ttk.Frame(self.tab_data)
        cont.pack(fill=tk.BOTH, expand=True, padx=10, pady=6)
        canvas = tk.Canvas(cont, highlightthickness=0)
        sb = ttk.Scrollbar(cont, orient=tk.VERTICAL, command=canvas.yview)
        canvas.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.tbl = ttk.Frame(canvas)
        canvas.create_window((0, 0), window=self.tbl, anchor="nw")
        self.tbl.bind("<Configure>", lambda _: canvas.configure(scrollregion=canvas.bbox("all")))
        self._refresh_table()

    def _add_year_ui(self):
        yr = self.new_yr.get().strip()
        if not (yr.isdigit() and len(yr) == 4):
            messagebox.showerror("Invalid", "Enter a 4-digit year."); return
        if yr in self.survey_data:
            messagebox.showinfo("Exists", f"{yr} already loaded."); return
        self._add_year(yr)
        self.year_cb.configure(values=self.years)
        self.d_year.set(yr)
        self.new_yr.set("")
        self._refresh_table()

    def _remove_year_ui(self):
        yr = self.d_year.get()
        if len(self.years) <= 1:
            messagebox.showwarning("Cannot Remove", "Need at least one year."); return
        if not messagebox.askyesno("Confirm", f"Remove all data for {yr}?"): return
        self._remove_year(yr)
        self.year_cb.configure(values=self.years)
        self.d_year.set(self.years[-1])
        self._refresh_table()

    def _refresh_table(self):
        for w in self.tbl.winfo_children():
            w.destroy()
        yr  = self.d_year.get()
        deg = self.d_deg.get()
        if yr not in self.survey_data:
            return
        ttk.Label(self.tbl, text="Exp (yrs)", style="H.TLabel", width=12
                  ).grid(row=0, column=0, padx=4, pady=3, sticky=tk.W)
        for c, col in enumerate(COLUMNS, 1):
            ttk.Label(self.tbl, text=col, style="H.TLabel", width=16
                      ).grid(row=0, column=c, padx=4, pady=3)
        for r, exp in enumerate(EXP_BINS, 1):
            ttk.Label(self.tbl, text=exp, width=12
                      ).grid(row=r, column=0, padx=4, pady=2, sticky=tk.W)
            for c, col in enumerate(COLUMNS, 1):
                ttk.Entry(self.tbl, textvariable=self.survey_data[yr][deg][exp][col], width=16
                          ).grid(row=r, column=c, padx=4, pady=2)

    # ── Calculator tab ─────────────────────────────────────────────────────────
    def _build_calc_tab(self):
        # Parameters
        plf = ttk.LabelFrame(self.tab_calc, text="Parameters")
        plf.pack(fill=tk.X, padx=10, pady=8)

        r1 = ttk.Frame(plf)
        r1.pack(fill=tk.X, padx=10, pady=4)
        ttk.Label(r1, text="Degree / Certification:", width=28).pack(side=tk.LEFT)
        self.c_deg = tk.StringVar(value=DEGREE_CATS[0])
        ttk.Combobox(r1, textvariable=self.c_deg, values=DEGREE_CATS,
                     state="readonly", width=26).pack(side=tk.LEFT, padx=6)

        r2 = ttk.Frame(plf)
        r2.pack(fill=tk.X, padx=10, pady=4)
        ttk.Label(r2, text="Years of Experience:", width=28).pack(side=tk.LEFT)
        self.c_exp = tk.StringVar(value="10")
        ttk.Spinbox(r2, from_=0, to=50, textvariable=self.c_exp, width=8
                    ).pack(side=tk.LEFT, padx=6)
        self.bin_lbl = ttk.Label(r2, text="", foreground="gray")
        self.bin_lbl.pack(side=tk.LEFT, padx=6)
        self.c_exp.trace_add("write", self._upd_bin)
        self._upd_bin()

        r3 = ttk.Frame(plf)
        r3.pack(fill=tk.X, padx=10, pady=4)
        ttk.Label(r3, text="Role:", width=28).pack(side=tk.LEFT)
        self.c_role = tk.StringVar(value=ROLES[0])
        role_cb = ttk.Combobox(r3, textvariable=self.c_role, values=ROLES,
                                state="readonly", width=22)
        role_cb.pack(side=tk.LEFT, padx=6)
        role_cb.bind("<<ComboboxSelected>>", self._on_role)

        # Admin pay
        self.adm_lf = ttk.LabelFrame(self.tab_calc, text="Administrative Pay Component")
        self.adm_lf.pack(fill=tk.X, padx=10, pady=2)

        ar1 = ttk.Frame(self.adm_lf)
        ar1.pack(fill=tk.X, padx=10, pady=4)
        ttk.Label(ar1, text="Admin Pay Type:", width=28).pack(side=tk.LEFT)
        self.adm_type = tk.StringVar(value="percentage")
        ttk.Radiobutton(ar1, text="% of base", variable=self.adm_type,
                        value="percentage", command=self._upd_adm_lbl).pack(side=tk.LEFT, padx=6)
        ttk.Radiobutton(ar1, text="Flat $",    variable=self.adm_type,
                        value="flat",       command=self._upd_adm_lbl).pack(side=tk.LEFT, padx=6)

        ar2 = ttk.Frame(self.adm_lf)
        ar2.pack(fill=tk.X, padx=10, pady=4)
        self.adm_lbl = ttk.Label(ar2, text="Admin Pay (%):", width=28)
        self.adm_lbl.pack(side=tk.LEFT)
        self.adm_val = tk.StringVar(value="0")
        self.adm_ent = ttk.Entry(ar2, textvariable=self.adm_val, width=12)
        self.adm_ent.pack(side=tk.LEFT, padx=6)
        ttk.Label(ar2,
                  text="Staff: 0%   Site Lead: 3%   Service Chief: 7%   Director: 10%",
                  foreground="gray").pack(side=tk.LEFT, padx=10)

        # Survey projection
        plf2 = ttk.LabelFrame(self.tab_calc, text="Survey Data Projection")
        plf2.pack(fill=tk.X, padx=10, pady=2)
        pr1 = ttk.Frame(plf2)
        pr1.pack(fill=tk.X, padx=10, pady=4)
        self.proj_on = tk.BooleanVar(value=True)
        ttk.Checkbutton(pr1, text="Project survey data forward to current year",
                        variable=self.proj_on).pack(side=tk.LEFT)
        pr2 = ttk.Frame(plf2)
        pr2.pack(fill=tk.X, padx=10, pady=4)
        ttk.Label(pr2, text="Annual growth rate (%):", width=28).pack(side=tk.LEFT)
        self.proj_rate = tk.StringVar(value="5.0")
        ttk.Entry(pr2, textvariable=self.proj_rate, width=8).pack(side=tk.LEFT, padx=6)
        ttk.Button(pr2, text="Calculate from Survey Data",
                   command=self._calc_growth_rate).pack(side=tk.LEFT, padx=10)
        self.proj_rate_lbl = ttk.Label(pr2, text="", foreground="gray")
        self.proj_rate_lbl.pack(side=tk.LEFT, padx=4)

        # Salary inputs
        slf = ttk.LabelFrame(self.tab_calc, text="Salary")
        slf.pack(fill=tk.X, padx=10, pady=6)

        sr1 = ttk.Frame(slf)
        sr1.pack(fill=tk.X, padx=10, pady=4)
        ttk.Label(sr1, text="Current Salary ($/yr):", width=28).pack(side=tk.LEFT)
        self.cur_sal = tk.StringVar()
        ttk.Entry(sr1, textvariable=self.cur_sal, width=16).pack(side=tk.LEFT, padx=6)

        sr2 = ttk.Frame(slf)
        sr2.pack(fill=tk.X, padx=10, pady=4)
        ttk.Label(sr2, text="Proposed Increase:", width=28).pack(side=tk.LEFT)
        self.inc_val = tk.StringVar()
        ttk.Entry(sr2, textvariable=self.inc_val, width=16).pack(side=tk.LEFT, padx=6)
        self.inc_type = tk.StringVar(value="dollar")
        ttk.Radiobutton(sr2, text="$ amount", variable=self.inc_type, value="dollar"
                        ).pack(side=tk.LEFT, padx=4)
        ttk.Radiobutton(sr2, text="% raise",  variable=self.inc_type, value="percent"
                        ).pack(side=tk.LEFT, padx=4)

        ttk.Button(self.tab_calc, text="  Calculate  ",
                   command=self._calc, padding=(8, 5)).pack(pady=8)

        rlf = ttk.LabelFrame(self.tab_calc, text="Results")
        rlf.pack(fill=tk.BOTH, expand=True, padx=10, pady=4)
        self.res = tk.Text(rlf, font=("Courier", 10), state=tk.DISABLED, background="#f8f8f8")
        rsby = ttk.Scrollbar(rlf, orient=tk.VERTICAL, command=self.res.yview)
        self.res.configure(yscrollcommand=rsby.set)
        rsby.pack(side=tk.RIGHT, fill=tk.Y)
        self.res.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        self._on_role(None)

    def _upd_bin(self, *_):
        try:
            b = exp_to_bin(int(self.c_exp.get()))
            self.bin_lbl.configure(text=f"→ survey bin: {b} yrs")
        except (ValueError, AttributeError):
            pass

    def _on_role(self, _):
        role = self.c_role.get()
        is_adm = role != "Physicist"
        state = tk.NORMAL if is_adm else tk.DISABLED
        self.adm_ent.configure(state=state)
        for child in self.adm_lf.winfo_children():
            for w in child.winfo_children():
                try:
                    w.configure(state=state)
                except Exception:
                    pass
        self.adm_val.set(str(ADMIN_PCT.get(role, 0)))

    def _upd_adm_lbl(self):
        lbl = "Admin Pay (%):" if self.adm_type.get() == "percentage" else "Admin Pay ($):"
        self.adm_lbl.configure(text=lbl)

    def _calc_growth_rate(self):
        """Compute CAGR from the oldest to newest loaded survey year using the
        Median salary for the currently selected degree and experience bin."""
        deg     = self.c_deg.get()
        exp_bin = exp_to_bin(int(self.c_exp.get()) if self.c_exp.get().isdigit() else 0)

        # Collect (year_int, median) pairs across loaded years
        points = []
        for yr in self.years:
            v = self.survey_data[yr][deg][exp_bin]["Median"].get().strip()
            if v:
                try:
                    points.append((int(yr), parse_dollar(v)))
                except ValueError:
                    pass

        if len(points) < 2:
            messagebox.showwarning(
                "Not Enough Data",
                f"Need Median values in at least 2 survey years for\n"
                f"  Degree: {deg}\n  Exp bin: {exp_bin}\n\n"
                f"Enter data in the Survey Data tab first.")
            return

        points.sort()
        yr_old, med_old = points[0]
        yr_new, med_new = points[-1]
        n = yr_new - yr_old

        cagr = (med_new / med_old) ** (1.0 / n) - 1
        self.proj_rate.set(f"{cagr * 100:.2f}")
        self.proj_rate_lbl.configure(
            text=f"CAGR {yr_old}→{yr_new}  "
                 f"({fmt_dollar(med_old)} → {fmt_dollar(med_new)},  "
                 f"n={n} yrs)")

    # ── Calculate ──────────────────────────────────────────────────────────────
    def _calc(self):
        try:
            years = int(self.c_exp.get())
        except ValueError:
            messagebox.showerror("Input Error", "Years of experience must be a whole number.")
            return

        deg     = self.c_deg.get()
        role    = self.c_role.get()
        exp_bin = exp_to_bin(years)

        try:
            cur = parse_dollar(self.cur_sal.get())
        except (ValueError, AttributeError):
            messagebox.showwarning("Missing Input", "Please enter your current salary.")
            return

        # Admin pay
        is_adm = role != "Physicist"
        adm_pay = 0.0
        adm_note = ""
        if is_adm:
            try:
                av = float(self.adm_val.get() or 0)
            except ValueError:
                av = 0.0
            if self.adm_type.get() == "percentage":
                adm_pay = cur * av / 100.0
                adm_note = f"{av:.1f}%"
            else:
                adm_pay = av
                adm_note = "flat"

        cur_total = cur + adm_pay

        # Proposed
        prop_sal = None
        prop_total = None
        prop_note = ""
        inc_str = self.inc_val.get().strip()
        if inc_str:
            try:
                iv = parse_dollar(inc_str)
                if self.inc_type.get() == "percent":
                    inc_amt = cur * iv / 100.0
                    prop_sal = cur + inc_amt
                    prop_note = f"+{iv:.1f}% = {fmt_dollar(inc_amt)}"
                else:
                    prop_sal = cur + iv
                    prop_note = f"+{fmt_dollar(iv)}"
                prop_total = prop_sal + adm_pay
            except ValueError:
                messagebox.showerror("Input Error", f"Invalid increase: '{inc_str}'")
                return

        # Pull survey rows across all years
        current_year = datetime.date.today().year
        do_proj = self.proj_on.get()
        try:
            growth_rate = float(self.proj_rate.get()) / 100.0
        except ValueError:
            growth_rate = 0.0

        year_rows = {}
        for yr in sorted(self.years):
            row = {}
            for col in COLUMNS:
                v = self.survey_data[yr][deg][exp_bin][col].get().strip()
                if v:
                    try:
                        row[col] = parse_dollar(v)
                    except ValueError:
                        pass
            if row:
                # Apply forward projection if enabled
                if do_proj and growth_rate > 0:
                    try:
                        years_elapsed = current_year - int(yr)
                    except ValueError:
                        years_elapsed = 0
                    if years_elapsed > 0:
                        factor = (1 + growth_rate) ** years_elapsed
                        row = {col: v * factor for col, v in row.items()}
                year_rows[yr] = row

        lines = self._render(deg, exp_bin, role, years,
                             cur, adm_pay, adm_note, cur_total,
                             prop_sal, prop_note, prop_total,
                             year_rows, do_proj, growth_rate, current_year)

        self.res.configure(state=tk.NORMAL)
        self.res.delete("1.0", tk.END)
        self.res.insert(tk.END, "\n".join(lines))
        self.res.configure(state=tk.DISABLED)

    def _render(self, deg, exp_bin, role, years_exp,
                cur, adm_pay, adm_note, cur_total,
                prop_sal, prop_note, prop_total, year_rows,
                do_proj=False, growth_rate=0.0, current_year=None):
        W = 62
        has_adm = adm_pay > 0

        proj_note = ""
        if do_proj and growth_rate > 0 and current_year:
            proj_note = f"  Projection  : {growth_rate*100:.1f}%/yr compounded to {current_year}"

        L = ["=" * W,
             "  AAPM SALARY SURVEY ANALYSIS",
             "=" * W,
             f"  Degree   : {deg}",
             f"  Exp      : {years_exp} yrs  (survey bin: {exp_bin})",
             f"  Role     : {role}",]
        if proj_note:
            L.append(proj_note)
        L.append("-" * W)
        L += ["",
              "  YOUR SALARY",
              f"  {'Current salary':<32}: {fmt_dollar(cur):>14}"]

        if has_adm:
            L.append(f"  {'Admin pay (' + adm_note + ')':<32}: {fmt_dollar(adm_pay):>14}")
            L.append("  " + "·" * 48)
            L.append(f"  {'Current total compensation':<32}: {fmt_dollar(cur_total):>14}")

        if prop_sal is not None:
            L.append("")
            L.append(f"  {'Proposed base  (' + prop_note + ')':<32}: {fmt_dollar(prop_sal):>14}")
            if has_adm:
                L.append(f"  {'Proposed total compensation':<32}: {fmt_dollar(prop_total):>14}")
            delta = prop_sal - cur
            pct   = delta / cur * 100 if cur else 0
            L.append(f"  {'Increase':<32}: {fmt_dollar(delta):>14}  ({pct:+.1f}%)")

        if not year_rows:
            L += ["", "  No survey data found for this degree/bin.",
                  "  Go to Survey Data tab to enter values.", "=" * W]
            return L

        # Multi-year survey table
        yr_list = sorted(year_rows.keys())
        col_w = 13
        L += ["",
              "=" * W,
              f"  SURVEY DATA  ({exp_bin} yrs exp,  {deg})",
              "=" * W]

        hdr = f"  {'':16}" + "".join(f"{yr:>{col_w}}" for yr in yr_list)
        L.append(hdr)
        L.append("  " + "-" * (16 + len(yr_list) * col_w))

        for col in COLUMNS:
            row_s = f"  {col:<16}"
            for yr in yr_list:
                v = year_rows[yr].get(col)
                row_s += f"{(fmt_dollar(v) if v else '---'):>{col_w}}"
            L.append(row_s)

        L.append("  " + "-" * (16 + len(yr_list) * col_w))

        # Percentile positions for current
        comp_lbl = "Total comp" if has_adm else "Current"
        comp_val  = cur_total
        row_c = f"  {comp_lbl:<16}"
        for yr in yr_list:
            row_c += f"{estimate_percentile(comp_val, year_rows[yr]):>{col_w}}"
        L.append(row_c)

        # Percentile positions for proposed
        if prop_total is not None:
            prop_lbl = "Prop total" if has_adm else "Proposed"
            prop_val  = prop_total
            row_p = f"  {prop_lbl:<16}"
            for yr in yr_list:
                row_p += f"{estimate_percentile(prop_val, year_rows[yr]):>{col_w}}"
            L.append(row_p)

        # Visual bar using most recent year
        latest = yr_list[-1]
        lr = year_rows[latest]
        if all(k in lr for k in ("20th", "Median", "80th")):
            L += ["",
                  f"  Position on {latest} scale:",
                  f"  20th={fmt_dollar(lr['20th'])}   "
                  f"Median={fmt_dollar(lr['Median'])}   "
                  f"80th={fmt_dollar(lr['80th'])}",
                  "  " + "─" * 54,
                  bar_line(comp_val, lr, comp_lbl)]
            if prop_total is not None:
                p_lbl = "Prop total" if has_adm else "Proposed"
                L.append(bar_line(prop_total, lr, p_lbl))

        L.append("=" * W)
        return L


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    root = tk.Tk()
    SalarySurveyTool(root)
    root.mainloop()
