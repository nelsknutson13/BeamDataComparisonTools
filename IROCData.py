import tkinter as tk
from tkinter import filedialog, ttk, messagebox
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# ── helpers ──────────────────────────────────────────────────────────────────

def _make_listbox(parent, height=6, **grid_kw):
    frame = ttk.Frame(parent)
    frame.grid(**grid_kw, sticky="nsew")
    sb = ttk.Scrollbar(frame, orient="vertical")
    lb = tk.Listbox(frame, selectmode="extended", height=height,
                    exportselection=False, yscrollcommand=sb.set)
    sb.config(command=lb.yview)
    lb.pack(side="left", fill="both", expand=True)
    sb.pack(side="right", fill="y")
    return lb

def _ratio_str(lo_var, hi_var):
    try:
        return f"{1 + float(lo_var.get())/100:.2f} – {1 + float(hi_var.get())/100:.2f}"
    except ValueError:
        return ""

def _toggle(lb):
    if lb.curselection() == tuple(range(lb.size())):
        lb.selection_clear(0, tk.END)
    else:
        lb.selection_set(0, tk.END)

def _selected(lb, values):
    sel = lb.curselection()
    return [values[i] for i in sel] if sel else list(values)

# ── App ───────────────────────────────────────────────────────────────────────

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("IROC TLD Plotter")
        self.resizable(True, True)

        self._excel_data = {}   # sheet_name -> DataFrame
        self._current_fig = None
        self._sheet_values = []
        self._site_values  = []
        self._type_values  = []

        pad = {"padx": 8, "pady": 4}

        # ── file row ──
        ttk.Label(self, text="Excel file:").grid(row=0, column=0, sticky="w", **pad)
        DEFAULT_PATH = r"C:\Users\nknutson\Box\Knutson\Research\Projects underway\Radformation research\RADMC Photon VALIDATION\IROC Data\IROC phantoms.xlsx"
        self._path_var = tk.StringVar(value=DEFAULT_PATH)
        ttk.Entry(self, textvariable=self._path_var, width=55).grid(
            row=0, column=1, columnspan=2, sticky="ew", **pad)
        ttk.Button(self, text="Browse…", command=self._browse).grid(
            row=0, column=3, **pad)

        # ── sheet selection ──
        ttk.Label(self, text="Sheets:").grid(row=1, column=0, sticky="nw", **pad)
        self._sheet_lb = _make_listbox(self, height=5, row=1, column=1, columnspan=2)
        self._sheet_lb.bind("<<ListboxSelect>>", lambda _: self._on_sheet_change())
        ttk.Button(self, text="Toggle all", command=lambda: (
            _toggle(self._sheet_lb), self._on_sheet_change()
        )).grid(row=2, column=1, sticky="w", **pad)

        # ── site selection ──
        ttk.Label(self, text="Sites (Report):").grid(row=3, column=0, sticky="nw", **pad)
        self._site_lb = _make_listbox(self, height=6, row=3, column=1)
        ttk.Button(self, text="Toggle all",
                   command=lambda: _toggle(self._site_lb)).grid(
            row=4, column=1, sticky="w", **pad)

        # ── phantom type / body site selection ──
        ttk.Label(self, text="Phantom types\n(Body Site):").grid(
            row=3, column=2, sticky="nw", **pad)
        self._type_lb = _make_listbox(self, height=6, row=3, column=2)
        ttk.Button(self, text="Toggle all",
                   command=lambda: _toggle(self._type_lb)).grid(
            row=4, column=2, sticky="w", **pad)

        # ── options ──
        opts = ttk.LabelFrame(self, text="Options", padding=5)
        opts.grid(row=5, column=0, columnspan=4, sticky="ew", **pad)

        self._add_all_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(opts, text="Include 'All' group",
                        variable=self._add_all_var).grid(row=0, column=0, sticky="w")

        self._show_legend_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(opts, text="Show site legend",
                        variable=self._show_legend_var).grid(row=0, column=1, sticky="w", padx=(20, 0))

        ttk.Label(opts, text="Font size:").grid(row=0, column=2, sticky="w", padx=(20, 4))
        self._font_var = tk.StringVar(value="16")
        ttk.Entry(opts, textvariable=self._font_var, width=5).grid(row=0, column=3, sticky="w")

        ttk.Label(opts, text="Jitter:").grid(row=0, column=4, sticky="w", padx=(20, 4))
        self._jitter_var = tk.StringVar(value="0")
        ttk.Entry(opts, textvariable=self._jitter_var, width=5).grid(row=0, column=5, sticky="w")

        ttk.Label(opts, text="Fig W×H:").grid(row=0, column=6, sticky="w", padx=(20, 4))
        self._fig_w_var = tk.StringVar(value="6")
        self._fig_h_var = tk.StringVar(value="8")
        ttk.Entry(opts, textvariable=self._fig_w_var, width=4).grid(row=0, column=7, sticky="w")
        ttk.Label(opts, text="×").grid(row=0, column=8, sticky="w")
        ttk.Entry(opts, textvariable=self._fig_h_var, width=4).grid(row=0, column=9, sticky="w")

        # row 1 — per-phantom tolerances
        self._show_tol_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(opts, text="Show IROC tolerances",
                        variable=self._show_tol_var).grid(row=1, column=0, sticky="w", pady=(4, 0))
        ttk.Button(opts, text="Edit tolerances…",
                   command=self._open_tol_editor).grid(row=1, column=1, sticky="w",
                                                        padx=(12, 0), pady=(4, 0))
        # {phantom_type: (lo_pct, hi_pct)}  — populated on sheet load, editable via dialog
        self._tol_map = {}

        ttk.Label(opts, text="Y-axis label:").grid(row=2, column=0, sticky="w", padx=(0, 4), pady=(4, 0))
        self._ylabel_var = tk.StringVar(value="Difference (MC vs IROC) [%RX]")
        ttk.Entry(opts, textvariable=self._ylabel_var, width=40).grid(
            row=2, column=1, columnspan=5, sticky="ew", pady=(4, 0))

        ttk.Label(opts, text="Y min:").grid(row=3, column=0, sticky="w", padx=(0, 4), pady=(4, 0))
        self._ymin_var = tk.StringVar(value="-10")
        ttk.Entry(opts, textvariable=self._ymin_var, width=6).grid(row=3, column=1, sticky="w", pady=(4, 0))
        ttk.Label(opts, text="Y max:").grid(row=3, column=2, sticky="w", padx=(12, 4), pady=(4, 0))
        self._ymax_var = tk.StringVar(value="10")
        ttk.Entry(opts, textvariable=self._ymax_var, width=6).grid(row=3, column=3, sticky="w", pady=(4, 0))

        # ── buttons ──
        btn_frame = ttk.Frame(self)
        btn_frame.grid(row=6, column=0, columnspan=4, pady=8)
        ttk.Button(btn_frame, text="Plot",       command=self._plot).pack(side="left", padx=6)
        ttk.Button(btn_frame, text="Save Plot…", command=self._save).pack(side="left", padx=6)

        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)

        # Auto-load default file if it exists
        import os
        if os.path.exists(DEFAULT_PATH):
            self.after(100, self._load_file)

    # ── file / population ─────────────────────────────────────────────────────

    def _browse(self):
        path = filedialog.askopenfilename(
            filetypes=[("Excel files", "*.xlsx *.xls")])
        if not path:
            return
        self._path_var.set(path)
        self._load_file()

    def _load_file(self):
        path = self._path_var.get().strip()
        if not path:
            return
        try:
            self._excel_data = pd.read_excel(path, sheet_name=None)
        except Exception as e:
            messagebox.showerror("Error loading file", str(e))
            return

        self._sheet_values = sorted(self._excel_data.keys())
        self._sheet_lb.delete(0, tk.END)
        for s in self._sheet_values:
            self._sheet_lb.insert(tk.END, s)
        self._sheet_lb.selection_set(0, tk.END)
        self._on_sheet_change()

    def _on_sheet_change(self):
        sheets = _selected(self._sheet_lb, self._sheet_values)
        if not sheets:
            return
        frames = []
        for s in sheets:
            df = self._excel_data.get(s)
            if df is not None:
                frames.append(df)
        if not frames:
            return
        combined = pd.concat(frames, ignore_index=True)

        # populate sites
        if "Report" in combined.columns:
            sites = sorted(combined["Report"].dropna().astype(str).unique())
        else:
            sites = []
        self._site_values = sites
        self._site_lb.delete(0, tk.END)
        for s in sites:
            self._site_lb.insert(tk.END, s)
        self._site_lb.selection_set(0, tk.END)

        # populate phantom types
        if "Body Site" in combined.columns:
            types = sorted(combined["Body Site"].dropna().astype(str).unique())
        else:
            types = []
        self._type_values = types
        self._type_lb.delete(0, tk.END)
        for t in types:
            self._type_lb.insert(tk.END, t)
        self._type_lb.selection_set(0, tk.END)

        # Pre-fill IROC tolerances for known phantom types (don't overwrite user edits)
        for t in types:
            if t in self._tol_map:
                continue
            tl = t.lower()
            if "lung" in tl:
                self._tol_map[t] = (-8.0, 5.0)   # 0.92 – 1.05
            elif any(k in tl for k in ("prostate", "h&n", "hn", "head")):
                self._tol_map[t] = (-7.0, 7.0)   # 0.93 – 1.07
            else:
                self._tol_map[t] = (-5.0, 5.0)   # generic default

    # ── plot ──────────────────────────────────────────────────────────────────

    def _open_tol_editor(self):
        types = self._type_values if self._type_values else list(self._tol_map.keys())
        if not types:
            messagebox.showinfo("No phantom types", "Load a file first.")
            return

        dlg = tk.Toplevel(self)
        dlg.title("IROC Tolerances per Phantom Type")
        dlg.resizable(False, False)
        dlg.grab_set()

        ttk.Label(dlg, text="Phantom Type", width=24).grid(row=0, column=0, padx=10, pady=4, sticky="w")
        ttk.Label(dlg, text="Lower %RX", width=10).grid(row=0, column=1, padx=6, pady=4)
        ttk.Label(dlg, text="Upper %RX", width=10).grid(row=0, column=2, padx=6, pady=4)
        ttk.Label(dlg, text="(ratio range)").grid(row=0, column=3, padx=6, pady=4, sticky="w")

        entries = {}
        for i, t in enumerate(types):
            lo, hi = self._tol_map.get(t, (-5.0, 5.0))
            ttk.Label(dlg, text=t, width=24).grid(row=i+1, column=0, padx=10, pady=3, sticky="w")
            lo_var = tk.StringVar(value=str(lo))
            hi_var = tk.StringVar(value=str(hi))
            ttk.Entry(dlg, textvariable=lo_var, width=8).grid(row=i+1, column=1, padx=6, pady=3)
            ttk.Entry(dlg, textvariable=hi_var, width=8).grid(row=i+1, column=2, padx=6, pady=3)

            lbl = ttk.Label(dlg, text="", width=12)
            lbl.grid(row=i+1, column=3, padx=6, pady=3, sticky="w")
            entries[t] = (lo_var, hi_var)

        ratio_lbls = {}
        for i, t in enumerate(types):
            lbl = dlg.grid_slaves(row=i+1, column=3)[0]
            ratio_lbls[i+1] = lbl
            lo_var, hi_var = entries[t]
            lo_var.trace_add("write", lambda *_, lv=lo_var, hv=hi_var, r=i+1: (
                ratio_lbls[r].config(text=_ratio_str(lv, hv))))
            hi_var.trace_add("write", lambda *_, lv=lo_var, hv=hi_var, r=i+1: (
                ratio_lbls[r].config(text=_ratio_str(lv, hv))))
            ratio_lbls[i+1].config(text=_ratio_str(lo_var, hi_var))

        def _save(*_):
            for t, (lo_v, hi_v) in entries.items():
                try:
                    self._tol_map[t] = (float(lo_v.get()), float(hi_v.get()))
                except ValueError:
                    pass
            dlg.destroy()

        n = len(types)
        btn = ttk.Frame(dlg)
        btn.grid(row=n+1, column=0, columnspan=4, pady=8)
        ttk.Button(btn, text="OK",     command=_save       ).pack(side="left", padx=6)
        ttk.Button(btn, text="Cancel", command=dlg.destroy ).pack(side="left", padx=6)

    def _plot(self):
        sheets = _selected(self._sheet_lb, self._sheet_values)
        sites  = _selected(self._site_lb,  self._site_values)
        types  = _selected(self._type_lb,  self._type_values)

        if not sheets:
            messagebox.showwarning("Nothing selected", "Select at least one sheet.")
            return

        try:
            font_size = float(self._font_var.get())
        except ValueError:
            font_size = 16
        try:
            jitter = float(self._jitter_var.get())
        except ValueError:
            jitter = 0.08

        frames = [self._excel_data[s] for s in sheets if s in self._excel_data]
        df = pd.concat(frames, ignore_index=True)
        df = df.dropna(subset=["Body Site", "Report", "%RX Difference"])
        df["Report"] = df["Report"].apply(
            lambda x: x.split("_")[0] if isinstance(x, str) else x)

        # apply filters
        df = df[df["Report"].astype(str).isin(sites)]
        df = df[df["Body Site"].astype(str).isin(types)]

        if df.empty:
            messagebox.showwarning("No data", "No data after applying filters.")
            return

        if self._add_all_var.get():
            df_all = df.copy()
            df_all["Body Site"] = "All"
            df = pd.concat([df, df_all], ignore_index=True)

        unique_sites = sorted(df["Report"].unique())
        markers = ["o", "s", "^", "D", "P", "X", "*", "v", "<", ">"]
        marker_map = {s: markers[i % len(markers)] for i, s in enumerate(unique_sites)}

        plt.rcParams.update({
            "font.size":             font_size,
            "axes.labelsize":        font_size,
            "xtick.labelsize":       font_size,
            "ytick.labelsize":       font_size,
            "legend.fontsize":       max(font_size - 2, 8),
            "legend.title_fontsize": max(font_size - 2, 8),
            "font.weight":           "bold",
            "axes.labelweight":      "bold",
        })

        try:
            fig_w = float(self._fig_w_var.get())
            fig_h = float(self._fig_h_var.get())
        except ValueError:
            fig_w, fig_h = 10, 8
        self._current_fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        sns.set_theme(style="whitegrid", palette="colorblind")
        flier_style = dict(marker="o", markerfacecolor="none",
                           markersize=5, linestyle="none", color="gray")
        sns.boxplot(data=df, x="Body Site", y="%RX Difference",
                    whis=1.5, flierprops=flier_style, ax=ax)

        body_sites = list(df["Body Site"].unique())
        for i, body_site in enumerate(body_sites):
            subset = df[df["Body Site"] == body_site]
            for site in unique_sites:
                site_data = subset[subset["Report"] == site]
                if not site_data.empty:
                    xs = np.random.normal(loc=i, scale=jitter, size=len(site_data))
                    ax.scatter(xs, site_data["%RX Difference"],
                               marker=marker_map[site],
                               edgecolor="black", facecolor="none",
                               s=70, alpha=0.8, linewidths=2.0)

        ax.axhline(0, color="gray", linestyle="--")

        # Per-phantom tolerance segments + legend
        tol_handles = []
        if self._show_tol_var.get():
            # Collect unique (lo, hi) pairs and assign colors
            tol_colors = {}
            color_cycle = ["#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4"]
            for body_site in body_sites:
                if body_site == "All":
                    continue
                pair = self._tol_map.get(body_site)
                if pair and pair not in tol_colors:
                    tol_colors[pair] = color_cycle[len(tol_colors) % len(color_cycle)]

            # Draw segments and build legend handles
            drawn_pairs = {}
            for i, body_site in enumerate(body_sites):
                if body_site == "All":
                    continue
                pair = self._tol_map.get(body_site)
                if not pair:
                    continue
                lo, hi = pair
                col = tol_colors[pair]
                ax.hlines([lo, hi], i - 0.45, i + 0.45,
                          colors=col, linestyles="--", linewidth=2.0, zorder=5)
                if pair not in drawn_pairs:
                    # find all phantom types using this tolerance
                    sites_for_pair = [b for b in body_sites
                                      if b != "All" and self._tol_map.get(b) == pair]
                    lo_r = 1 + lo / 100
                    hi_r = 1 + hi / 100
                    lbl = f"{', '.join(sites_for_pair)}: {lo_r:.2f}–{hi_r:.2f}"
                    tol_handles.append(
                        plt.Line2D([], [], linestyle="--", color=col,
                                   linewidth=2.0, label=lbl))
                    drawn_pairs[pair] = True

        ax.set_ylabel(self._ylabel_var.get() or "Difference [%RX]")
        ax.set_xlabel("Phantom Type")
        try:
            ax.set_ylim(float(self._ymin_var.get()), float(self._ymax_var.get()))
        except ValueError:
            pass

        all_handles = tol_handles.copy()
        if self._show_legend_var.get():
            all_handles += [
                plt.Line2D([], [], marker=marker_map[s], linestyle="None",
                           markersize=10, label=s,
                           markerfacecolor="none", markeredgecolor="black",
                           markeredgewidth=2.0)
                for s in unique_sites
            ]
        if all_handles:
            ax.legend(handles=all_handles, loc="upper right", frameon=True)
        plt.tight_layout(pad=0.5)
        self._current_fig.subplots_adjust(left=0.15)
        plt.show()

    # ── save ─────────────────────────────────────────────────────────────────

    def _save(self):
        if self._current_fig is None:
            messagebox.showwarning("No plot", "Generate a plot first.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".tiff",
            filetypes=[("TIFF", "*.tiff *.tif"), ("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")])
        if path:
            self._current_fig.savefig(path, dpi=300, bbox_inches="tight")


App().mainloop()
