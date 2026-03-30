#TRS 398 calculator. Written by Nels Knutson.  I assume no responsibility for its accuracy use
#Use at own risk.  

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import math
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from datetime import datetime
# ===== Chamber data =====
# Mode A: kQ(TPR20/10) = a * exp[b * (TPR20,10 - 0.57)]
# L [cm]: cavity length along beam axis for FFF volume-averaging correction.
CHAMBERS = {
    # NE chambers
    "NE 2561/2611A":       {"a": 1.07699, "b": -0.08732, "L": 3.07},
    "NE 2571":             {"a": 1.08918, "b": -0.09222, "L": 2.41},
    # PTW Farmer-type
    "PTW 30010":           {"a": 1.12594, "b": -0.10740, "L": 2.30},
    "PTW 30011":           {"a": 1.10850, "b": -0.10107, "L": 2.30},
    "PTW 30012":           {"a": 1.12442, "b": -0.10415, "L": 2.30},
    "PTW 30013":           {"a": 1.18273, "b": -0.13256, "L": 2.30},
    # PTW Semiflex / PinPoint
    "PTW 31010 Semiflex":  {"a": 1.23755, "b": -0.15295, "L": 1.65},
    "PTW 31013 Semiflex":  {"a": 1.19297, "b": -0.13366, "L": 1.65},
    "PTW 31016 PinPoint3D":{"a": 1.11650, "b": -0.10841, "L": 0.29},
    "PTW 31021 Semiflex3D":{"a": 1.29612, "b": -0.16514, "L": 1.65},
    "PTW 31022 PinPoint3D":{"a": 1.14435, "b": -0.11130, "L": 0.29},
    # Capintec
    "Capintec PR-06C":     {"a": 1.06833, "b": -0.08262, "L": 2.27},
    # Exradin
    "Exradin A1SL":        {"a": 1.21633, "b": -0.13351, "L": 0.57},
    "Exradin A12":         {"a": 1.09783, "b": -0.09544, "L": 2.58},
    "Exradin A12S":        {"a": 1.11499, "b": -0.10057, "L": 2.00},
    "Exradin A18":         {"a": 1.10487, "b": -0.09670, "L": 2.25},
    "Exradin A19":         {"a": 1.12024, "b": -0.10493, "L": 3.17},
    "Exradin A26":         {"a": 1.09587, "b": -0.09383, "L": 0.24},
    "Exradin A28":         {"a": 1.12453, "b": -0.10278, "L": 0.25},
    # IBA
    "IBA FC65-G":          {"a": 1.09752, "b": -0.09642, "L": 2.31},
    "IBA FC65-P":          {"a": 1.12374, "b": -0.10784, "L": 2.31},
    "IBA FC23-C":          {"a": 1.09189, "b": -0.09346, "L": 0.88},
    "IBA CC13":            {"a": 1.11441, "b": -0.10260, "L": 0.58},
    "IBA CC25":            {"a": 1.08981, "b": -0.09254, "L": 1.00},
    # Sun Nuclear
    "Sun Nuclear SNC600c": {"a": 1.06800, "b": -0.08485, "L": 2.27},
    "Sun Nuclear SNC125c": {"a": 1.09700, "b": -0.09749, "L": 0.705},
}
# ks coefficients for pulsed beams (TRS-398 Table 10 / 4.VII)
KS_PULSED = {
    2.0: (2.337, -3.636, 2.299),
    2.5: (1.474, -1.587, 1.114),
    3.0: (1.198, -0.875, 0.677),
    3.5: (1.080, -0.542, 0.463),
    4.0: (1.022, -0.363, 0.341),
    5.0: (0.975, -0.188, 0.214),
}

def get_kq(chamber_label, tpr):
    data = CHAMBERS.get(chamber_label)
    if not data:
        raise ValueError(f"k_Q not coded for '{chamber_label}'.")
    a = data["a"]
    b = data["b"]
    num = 1.0 + math.exp((a - 0.57) / b)
    den = 1.0 + math.exp((a - tpr) / b)
    return num / den


def get_kvol(label, tpr, sdd_cm, fff):
    """FFF volume-averaging correction (k_vol)."""
    if not fff:
        return 1.0
    data = CHAMBERS.get(label)
    if not data or data.get("L") is None:
        raise ValueError(f"k_vol not available for '{label}'.")
    if sdd_cm <= 0:
        raise ValueError("SDD must be > 0 cm for k_vol.")
    L = data["L"]
    Q = tpr
    # TRS-398 Rev.1 style: depends on Q, SDD, and L (for FFF beams).
    return 1.0 + (0.0062 * Q - 0.0036) * (100.0 / sdd_cm) ** 2 * (L ** 2)

def mean(vals):
    return sum(vals) / len(vals) if vals else 0.0

def parse_row(entries, label):
    vals = []
    for e in entries:
        txt = e.get().strip()
        if not txt:
            raise ValueError(f"Missing value in {label}.")
        vals.append(float(txt))
    return vals

def compute():
    beam = entry_beam.get().strip()
    chamber = chamber_combo.get().strip()
    user = entry_user.get().strip()
    machine = entry_machine.get().strip()
    chamber_sn  = entry_ch_sn.get().strip()
    chamber_cal = entry_ch_cal.get().strip()
    em_model    = entry_em_model.get().strip()
    em_sn       = entry_em_sn.get().strip()
    em_cal      = entry_em_cal.get().strip()
    if not chamber:
        messagebox.showerror("Input error", "Select a chamber.")
        return

    try:
        # Calibration factor ND,w in cGy/nC
        NDw = float(entry_NDw.get().strip())
        # Beam quality & kQ
        tpr = float(entry_TPR.get().strip())
        kQ = get_kq(chamber, tpr)

        # k_TP
        T = float(entry_T.get().strip())
        p_raw = float(entry_p.get().strip())
        T0 = float(entry_T0.get().strip())
        P0 = float(entry_P0.get().strip())
        # P0 and p are in the same units (hPa or mmHg) — ratio is dimensionless
        kTP = ((273.15 + T) / (273.15 + T0)) * (P0 / p_raw)


        # Readings (all in nC):
        # three readings each: M(+ , V_high), M(- , V_high), M(+ , V_low)
        plus_Vh  = parse_row(M_plus_Vh_entries,  "M(+ , V_high)")
        minus_Vh = parse_row(M_minus_Vh_entries, "M(- , V_high)")
        plus_Vl  = parse_row(M_plus_Vl_entries,  "M(+ , V_low)")

        M_plus  = mean(plus_Vh)   # routine polarity at |V_high|
        M_minus = mean(minus_Vh)  # opposite polarity at |V_high|
        M_low   = mean(plus_Vl)   # routine polarity at |V_low|
        # NEW: string versions of the three readings for the summary
        plus_Vh_str  = ", ".join(f"{v:.6f}" for v in plus_Vh)
        minus_Vh_str = ", ".join(f"{v:.6f}" for v in minus_Vh)
        plus_Vl_str  = ", ".join(f"{v:.6f}" for v in plus_Vl)
        
        
        if M_plus <= 0 or M_low <= 0:
            raise ValueError("M(+ , V_high) and M(+ , V_low) must be > 0.")

        # Reference reading
        M_ref = M_plus

        # k_pol (routine = +)
        M_mean_pm = (abs(M_plus) + abs(M_minus)) / 2.0
        if M_mean_pm == 0:
            raise ValueError("Invalid M(+/-) for k_pol.")
        kpol = (abs(M_plus)+abs(M_minus) )/(2* M_plus)

        # Recombination k_s (pulsed, two-voltage, TRS-398 Table 10)
        Vh_in = float(entry_Vh.get().strip())
        Vl_in = float(entry_Vl.get().strip())
        V1 = abs(Vh_in)
        V2 = abs(Vl_in)
        if V1 <= 0 or V2 <= 0 or V1 <= V2:
            raise ValueError("Require |V_high| > |V_low| > 0.")

        ratio_V = V1 / V2
        R = M_plus / M_low  # M1/M2

        # choose nearest tabulated V1/V2 for pulsed beams
        table_ratio = min(KS_PULSED.keys(), key=lambda x: abs(x - ratio_V))
        a0, a1, a2 = KS_PULSED[table_ratio]

        k_s = a0 + a1 * R + a2 * (R ** 2)


        # k_vol — use whatever is in the entry (user-typed or from Calc k_vol button)
        s = entry_kvol.get().strip()
        kvol = float(s) if s else 1.0

        # k_elec
        kele = float(entry_kelec.get().strip())

        # Corrected reading at z_ref [nC]
        Mcorr = M_ref * kTP * kpol * k_s * kvol * kele

        # Dose at reference depth [cGy]
        DwQ = Mcorr * NDw * kQ

        # Dose at dmax from clinical PDD(z_ref) [%]
        Ddmax = None
        pdd_txt = entry_PDD.get().strip()
        if pdd_txt:
            PDD = float(pdd_txt)
            if PDD <= 0 or PDD > 100:
                raise ValueError("Clinical PDD must be in (0,100].")
            Ddmax = DwQ / (PDD / 100.0)

                # --- Post-adjust (optional, simple) ---
        post_vals = [float(e.get().strip()) for e in post_M_entries if e.get().strip()]
        post_Mcorr = None
        post_DwQ = None
        post_Ddmax = None
        post_line = ""

        if post_vals:
            post_Mref = mean(post_vals)
            ratio = post_Mref / M_ref if M_ref != 0 else float("nan")

            # reuse same correction factors
            post_Mcorr = post_Mref * kTP * kpol * k_s * kvol * kele
            post_DwQ = post_Mcorr * NDw * kQ

            if pdd_txt:
                PDD = float(pdd_txt)
                if PDD > 0:
                    post_Ddmax = post_DwQ / (PDD / 100.0)

            post_line = (
                f"Post-adjust M(mean) = {post_Mref:.6f} nC  "
                f"(ratio vs orig M_ref = {ratio:.4f})"
            )


        # Update UI labels
        lbl_kQ.config(text=f"{kQ:.5f}")
        lbl_kTP.config(text=f"{kTP:.5f}")
        lbl_kpol.config(text=f"{kpol:.5f}")
        lbl_ks.config(text=f"{k_s:.5f}")
        lbl_kvol_val.config(text=f"{kvol:.5f}")
        lbl_Mcorr.config(text=f"{Mcorr:.6f} nC")
        lbl_DwQ.config(text=f"{DwQ:.3f} cGy")
        lbl_Ddmax.config(text=f"{Ddmax:.3f} cGy" if Ddmax is not None else "—")

                # Summary text (cleaned)
        now_str = datetime.now().strftime("%Y-%m-%d %H:%M")

        lines = [
            f"Date/Time: {now_str}",
            f"Beam: {beam}",
            f"Chamber: {chamber}",
            f"Chamber SN#: {chamber_sn}",
            f"Chamber Cal Date: {chamber_cal}",
            f"Electrometer Model: {em_model}",
            f"Electrometer SN#: {em_sn}",
            f"Electrometer Cal Date: {em_cal}",
            f"N_D,w(60Co) = {NDw:.3f} cGy/nC",
            f"TPR20,10 = {tpr}  -> k_Q = {kQ:.5f}  (a={CHAMBERS[chamber]['a']}, b={CHAMBERS[chamber]['b']})",
            f"T = {T:.1f} °C, p = {p_raw:.1f} {pressure_unit_var.get()}, T0 = {T0:.1f} °C, P0 = {P0:.2f} {pressure_unit_var.get()}",

            # Three raw readings + mean for V_high (+)
            f"M(Vref , V_high) readings [nC] = {plus_Vh_str}",
            f"M_ref = mean[M(+ , V_high)] = {M_ref:.6f} nC",
            f"M(Vref , V_high) mean = {M_plus:.6f} nC",

            # Three raw readings + mean for V_high (−)
            f"M(-1*Vref , -1*V_high) readings [nC] = {minus_Vh_str}",
            f"M(-1*Vref , V_high) mean = {M_minus:.6f} nC",

            f"|V_high| = {V1:.0f} V, |V_low| = {V2:.0f} V",

            # Three raw readings + mean for V_low
            f"M(Vref/2 , V_low) readings [nC] = {plus_Vl_str}",
            f"M(+ , V_low) mean = {M_low:.6f} nC",

            f"k_TP = {kTP:.5f}",
            f"k_pol = {kpol:.5f}",
            f"k_s = {k_s:.5f}",
            f"k_vol = {kvol:.5f}",
            f"k_elec = {kele:.5f}",
            f"Total correction (Σk) = {kTP * kpol * k_s * kvol * kele:.5f}",
            f"M_corr = {Mcorr:.6f} nC",
            f"D_w,Q(z_ref) = {DwQ:.3f} cGy",
        ]



        if Ddmax is not None:
            lines.append(f"D_w,Q(d_max) = {Ddmax:.3f} cGy")

        if post_line:
            lines.append(post_line)
        if post_DwQ is not None:
            lines.append(f"Post D_w,Q(z_ref) = {post_DwQ:.3f} cGy")
        if post_Ddmax is not None:
            lines.append(f"Post D_w,Q(d_max) = {post_Ddmax:.3f} cGy")


        txt_summary.config(state=tk.NORMAL)
        txt_summary.delete("1.0", tk.END)
        txt_summary.insert(tk.END, "\n".join(lines))
        txt_summary.config(state=tk.DISABLED)
        
        # stash for PDF export
        global last_report
        last_report = {
            "user": user,
            "machine": machine,
            "beam": beam,
            "chamber": chamber,
            "lines": lines,
        }
        
    except ValueError as e:
        messagebox.showerror("Input error", str(e))

#PDF Writer
def save_pdf_report():
    # Ensure we have something to save
    global last_report
    if "last_report" not in globals():
        messagebox.showerror("No data", "Please click Compute before saving a PDF.")
        return

    data = last_report
    if not data or "lines" not in data:
        messagebox.showerror("No data", "Please click Compute before saving a PDF.")
        return

    # Build default filename: TRS398_Machine_Beam_YYYYMMDD_HHMM.pdf
    today_str = datetime.now().strftime("%Y%m%d_%H%M")
    machine_safe = (data.get("machine") or "Machine").replace(" ", "_")
    beam_safe = (data.get("beam") or "Beam").replace(" ", "_")
    default_name = f"TRS398_{machine_safe}_{beam_safe}_{today_str}.pdf"
    
    filename = filedialog.asksaveasfilename(
        defaultextension=".pdf",
        filetypes=[("PDF files", "*.pdf")],
        initialfile=default_name,
        title="Save calibration report as PDF",
    )
    if not filename:
        return  # user cancelled

    c = canvas.Canvas(filename, pagesize=A4)
    width, height = A4

    x_left = 50
    y = height - 50

    # Title
    c.setFont("Helvetica-Bold", 14)
    c.drawString(x_left, y, "TRS-398 Rev.1 Photon Calibration Report")
    y -= 24

    # Meta block
    c.setFont("Helvetica", 10)
    if data.get("user"):
        c.drawString(x_left, y, f"User: {data['user']}")
        y -= 14
    if data.get("machine"):
        c.drawString(x_left, y, f"Machine: {data['machine']}")
        y -= 14
    if data.get("beam"):
        c.drawString(x_left, y, f"Beam: {data['beam']}")
        y -= 14
    if data.get("chamber"):
        c.drawString(x_left, y, f"Chamber: {data['chamber']}")
        y -= 18

    # Summary lines from the Details box
    c.setFont("Helvetica", 9)
    for line in data["lines"]:
        if y < 60:
            c.showPage()
            c.setFont("Helvetica", 9)
            y = height - 50
        c.drawString(x_left, y, line)
        y -= 12

    c.showPage()
    c.save()

    messagebox.showinfo("Saved", f"PDF report saved to:\n{filename}")


# ===== UI =====

root = tk.Tk()
root.title("TRS-398 Rev.1 – Photon Worksheet 6.9 ")
def on_close():
    try:
        root.quit()
        root.destroy()
    finally:
        raise SystemExit

root.protocol("WM_DELETE_WINDOW", on_close)
# Scrollable container
container = ttk.Frame(root)
container.pack(fill="both", expand=True)

canvas = tk.Canvas(container)
vscroll = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
canvas.configure(yscrollcommand=vscroll.set)

vscroll.pack(side="right", fill="y")
canvas.pack(side="left", fill="both", expand=True)

# This is your old 'main' frame, now inside the scrollable canvas
main = ttk.Frame(canvas, padding=8)
canvas.create_window((0, 0), window=main, anchor="nw")

def _on_main_configure(event):
    canvas.configure(scrollregion=canvas.bbox("all"))

main.bind("<Configure>", _on_main_configure)

# keep your existing column stretch on main
for c in range(5):
    main.columnconfigure(c, weight=1)

r = 0
# User
ttk.Label(main, text="User:").grid(row=r, column=0, sticky="w", pady=(2,2))
entry_user = ttk.Entry(main, width=20)
entry_user.grid(row=r, column=1, sticky="w", pady=(2,2))
# Machine
ttk.Label(main, text="Machine:").grid(row=r, column=2, sticky="w", pady=(2,2))
entry_machine = ttk.Entry(main, width=20)
entry_machine.grid(row=r, column=3, sticky="w", pady=(2,2))
r += 1

# Beam
ttk.Label(main, text="Beam Energy:").grid(row=r, column=0, sticky="w", pady=(2,2))
entry_beam = ttk.Entry(main, width=14)
entry_beam.grid(row=r, column=1, sticky="w", pady=(2,2))
r += 1

# Chamber
ttk.Label(main, text="Chamber:").grid(row=r, column=2, sticky="w", pady=(2,2))
chamber_combo = ttk.Combobox(
    main,
    values=sorted(CHAMBERS.keys()),
    state="readonly",
    width=24,
)
chamber_combo.grid(row=r, column=3, columnspan=2, sticky="w", pady=(2,2))
r += 1
# NEW: Chamber SN and Cal Date
ttk.Label(main, text="Chamber SN#:").grid(row=r, column=2, sticky="w", pady=(2,2))
entry_ch_sn = ttk.Entry(main, width=20)
entry_ch_sn.grid(row=r, column=3, sticky="w", pady=(2,2))
r += 1
ttk.Label(main, text="Chamber Cal Date:").grid(row=r, column=2, sticky="w", pady=(2,2))
entry_ch_cal = ttk.Entry(main, width=20)
entry_ch_cal.grid(row=r, column=3, sticky="w", pady=(2,2))
r += 1
# N_D,w line
ttk.Label(main, text="N_D,w(60Co) [cGy/nC]:").grid(row=r, column=0, sticky="w", pady=(2,2))
entry_NDw = ttk.Entry(main, width=14)
entry_NDw.grid(row=r, column=1, sticky="w", pady=(2,2))
r += 1

# Electrometer factor
ttk.Label(main, text="k_elec:").grid(row=r, column=0, sticky="w", pady=(2,2))
entry_kelec = ttk.Entry(main, width=8)
entry_kelec.insert(0, "1.000")
entry_kelec.grid(row=r, column=1, sticky="w", pady=(2,2))
r += 1
# NEW: Electrometer model / SN / Cal Date
ttk.Label(main, text="Electrometer Model:").grid(row=r, column=0, sticky="w", pady=(2,2))
entry_em_model = ttk.Entry(main, width=20)
entry_em_model.grid(row=r, column=1, sticky="w", pady=(2,2))
r += 1

ttk.Label(main, text="Electrometer SN#:").grid(row=r, column=0, sticky="w", pady=(2,2))
entry_em_sn = ttk.Entry(main, width=20)
entry_em_sn.grid(row=r, column=1, sticky="w", pady=(2,2))
r += 1

ttk.Label(main, text="Electrometer Cal Date:").grid(row=r, column=0, sticky="w", pady=(2,2))
entry_em_cal = ttk.Entry(main, width=20)
entry_em_cal.grid(row=r, column=1, sticky="w", pady=(2,2))
r += 1
# T0
ttk.Label(main, text="T0 [°C] (from ND,w cert):").grid(row=r, column=2, sticky="w", pady=(2,2))
entry_T0 = ttk.Entry(main, width=8)
entry_T0.insert(0, "22.0")
entry_T0.grid(row=r, column=3, sticky="w", pady=(2,2))
r += 1

# P0
lbl_P0_unit = ttk.Label(main, text="P0 [hPa] (from ND,w cert):")
lbl_P0_unit.grid(row=r, column=2, sticky="w", pady=(2,2))
entry_P0 = ttk.Entry(main, width=8)
entry_P0.insert(0, "1013.25")
entry_P0.grid(row=r, column=3, sticky="w", pady=(2,2))
r += 1

# TPR / k_Q
ttk.Label(main, text="TPR20,10:").grid(row=r, column=0, sticky="w", pady=(2,2))
entry_TPR = ttk.Entry(main, width=8)
entry_TPR.grid(row=r, column=1, sticky="w", pady=(2,2))
ttk.Label(main, text="k_Q:").grid(row=r, column=2, sticky="e", pady=(2,2))
lbl_kQ = ttk.Label(main, text="—")
lbl_kQ.grid(row=r, column=3, sticky="w", pady=(2,2))
r += 1

# Ambient T, p
ttk.Label(main, text="T [°C]:").grid(row=r, column=0, sticky="w", pady=(2,2))
entry_T = ttk.Entry(main, width=8)
entry_T.grid(row=r, column=1, sticky="w", pady=(2,2))

pressure_unit_var = tk.StringVar(value="hPa")
lbl_p_unit = ttk.Label(main, text="p [hPa]:")
lbl_p_unit.grid(row=r, column=2, sticky="w", pady=(2,2))
entry_p = ttk.Entry(main, width=8)
entry_p.grid(row=r, column=3, sticky="w", pady=(2,2))
def _on_unit_change():
    u = pressure_unit_var.get()
    lbl_p_unit.config(text=f"p [{u}]:")
    lbl_P0_unit.config(text=f"P0 [{u}] (from ND,w cert):")
    entry_p.delete(0, tk.END)
    entry_p.insert(0, "760.00" if u == "mmHg" else "1013.25")
    entry_P0.delete(0, tk.END)
    entry_P0.insert(0, "760.00" if u == "mmHg" else "1013.25")
ttk.Radiobutton(main, text="hPa",  variable=pressure_unit_var, value="hPa",  command=_on_unit_change).grid(row=r, column=4, sticky="w", pady=(2,2))
r += 1
ttk.Radiobutton(main, text="mmHg", variable=pressure_unit_var, value="mmHg", command=_on_unit_change).grid(row=r, column=4, sticky="w", pady=(2,2))

ttk.Label(main, text="k_TP:").grid(row=r, column=2, sticky="e", pady=(2,2))
lbl_kTP = ttk.Label(main, text="—")
lbl_kTP.grid(row=r, column=3, sticky="w", pady=(2,2))
r += 1

# Voltages
ttk.Label(main, text="V_high [V] (signed):").grid(row=r, column=0, sticky="w", pady=(2,2))
entry_Vh = ttk.Entry(main, width=8)
entry_Vh.grid(row=r, column=1, sticky="w", pady=(2,2))

ttk.Label(main, text="V_low [V] (signed):").grid(row=r, column=2, sticky="w", pady=(2,2))
entry_Vl = ttk.Entry(main, width=8)
entry_Vl.grid(row=r, column=3, sticky="w", pady=(2,2))
r += 1

# Readings block header
ttk.Label(main, text="ABS Value of Readings [nC]").grid(
    row=r, column=1, columnspan=3, sticky="w", pady=(10, 0)
)
r += 1

# Column headers for readings
ttk.Label(main, text="").grid(row=r, column=0)
ttk.Label(main, text="Rdg 1").grid(row=r, column=1, sticky="w", pady=(2,2))
ttk.Label(main, text="Rdg 2").grid(row=r, column=2, sticky="w", pady=(2,2))
ttk.Label(main, text="Rdg 3").grid(row=r, column=3, sticky="w", pady=(2,2))
r += 1

# M(Vref , V_high)
ttk.Label(main, text="M(Vref , V_high) [nC]:").grid(
    row=r, column=0, sticky="w", pady=(2,2)
)
M_plus_Vh_entries = []
for c in range(1, 4):
    e = ttk.Entry(main, width=8)
    e.grid(row=r, column=c, sticky="w", pady=(2,2))
    M_plus_Vh_entries.append(e)
r += 1

# M(-1*Vref , -1*V_high)
ttk.Label(main, text="M(-1*Vref , -1*V_high) [nC]:").grid(
    row=r, column=0, sticky="w", pady=(2,2)
)
M_minus_Vh_entries = []
for c in range(1, 4):
    e = ttk.Entry(main, width=8)
    e.grid(row=r, column=c, sticky="w", pady=(2,2))
    M_minus_Vh_entries.append(e)
r += 1

# M(Vref/2 , V_low)
ttk.Label(main, text="M(Vref/2 , V_low) [nC]:").grid(
    row=r, column=0, sticky="w", pady=(2,2)
)
M_plus_Vl_entries = []
for c in range(1, 4):
    e = ttk.Entry(main, width=8)
    e.grid(row=r, column=c, sticky="w", pady=(2,2))
    M_plus_Vl_entries.append(e)
r += 1

# k_pol, k_s
ttk.Label(main, text="k_pol:").grid(row=r, column=2, sticky="e", pady=(2,2))
lbl_kpol = ttk.Label(main, text="—")
lbl_kpol.grid(row=r, column=3, sticky="w", pady=(2,2))
r += 1

ttk.Label(main, text="k_s:").grid(row=r, column=2, sticky="e", pady=(2,2))
lbl_ks = ttk.Label(main, text="—")
lbl_ks.grid(row=r, column=3, sticky="w", pady=(2,2))
r += 1

# k_vol
ttk.Label(main, text="SDD [cm]:").grid(row=r, column=0, sticky="w", pady=(2,2))
entry_SDD = ttk.Entry(main, width=8)
entry_SDD.insert(0, "110")
entry_SDD.grid(row=r, column=1, sticky="w", pady=(2,2))

def calc_kvol():
    try:
        chamber = chamber_combo.get().strip()
        tpr = float(entry_TPR.get().strip())
        sdd = float(entry_SDD.get().strip())
        kvol = get_kvol(chamber, tpr, sdd, True)
        entry_kvol.delete(0, tk.END)
        entry_kvol.insert(0, f"{kvol:.5f}")
    except Exception as e:
        messagebox.showerror("k_vol error", str(e))

ttk.Button(main, text="Calc k_vol", command=calc_kvol).grid(row=r, column=2, sticky="w", pady=(2,2))
r += 1

ttk.Label(main, text="k_vol:").grid(row=r, column=0, sticky="w", pady=(2,2))
entry_kvol = ttk.Entry(main, width=8)
entry_kvol.insert(0, "1.000")
entry_kvol.grid(row=r, column=1, sticky="w", pady=(2,2))
lbl_kvol_val = ttk.Label(main, text="—")
lbl_kvol_val.grid(row=r, column=2, sticky="w", pady=(2,2))
r += 1

# Clinical PDD for Dmax
ttk.Label(main, text="Clinical PDD or TMR at z_ref [%] (for Dmax):").grid(
    row=r, column=0, columnspan=2, sticky="w", pady=(2,2)
)
entry_PDD = ttk.Entry(main, width=8)
entry_PDD.grid(row=r, column=2, sticky="w", pady=(2,2))
r += 1

# Outputs
ttk.Label(main, text="M_corr:").grid(row=r, column=0, sticky="w", pady=(2,2))
lbl_Mcorr = ttk.Label(main, text="—")
lbl_Mcorr.grid(row=r, column=1, sticky="w", pady=(2,2))

ttk.Label(main, text="D_w,Q(z_ref):").grid(row=r, column=2, sticky="e", pady=(2,2))
lbl_DwQ = ttk.Label(main, text="—")
lbl_DwQ.grid(row=r, column=3, sticky="w", pady=(2,2))
r += 1

ttk.Label(main, text="D_w,Q(d_max):").grid(row=r, column=2, sticky="e", pady=(2,2))
lbl_Ddmax = ttk.Label(main, text="—")
lbl_Ddmax.grid(row=r, column=3, sticky="w", pady=(2,2))
r += 1

# Post-adjust
ttk.Label(main, text="Post-adjust M(Vref , V_high) [nC]:").grid(
    row=r, column=0, sticky="w", pady=(2,2)
)
post_M_entries = []
for c in range(1, 4):
    e = ttk.Entry(main, width=8)
    e.grid(row=r, column=c, sticky="w", pady=(2,2))
    post_M_entries.append(e)
r += 1

# Summary box
ttk.Label(main, text="Details:").grid(row=r, column=0, sticky="nw", pady=(2,2))
txt_summary = tk.Text(main, width=80, height=10, state=tk.DISABLED)
txt_summary.grid(row=r, column=1, columnspan=4, sticky="nsew", pady=(2,2))
r += 1

ttk.Button(main, text="Compute", command=compute).grid(
    row=r, column=0, pady=(5,5), sticky="w")

ttk.Button(main, text="Save PDF Report", command=save_pdf_report).grid(
    row=r, column=1, pady=(5,5), sticky="w")


if __name__ == "__main__":
    root.mainloop()
