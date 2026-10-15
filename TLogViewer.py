"""
Trajectory Log Viewer
Displays actual delivered parameters from a Varian trajectory log file.
Requires: pylinac, matplotlib, tkinter
"""

import tkinter as tk
from tkinter import filedialog, ttk, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np

try:
    from pylinac.log_analyzer import TrajectoryLog
except ImportError:
    TrajectoryLog = None


class TLogViewer(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Trajectory Log Viewer")
        self.resizable(True, True)
        self._log = None

        pad = {"padx": 8, "pady": 4}

        # ── file row ──
        ttk.Label(self, text="Trajectory log:").grid(row=0, column=0, sticky="w", **pad)
        self._path_var = tk.StringVar()
        ttk.Entry(self, textvariable=self._path_var, width=60).grid(
            row=0, column=1, sticky="ew", **pad)
        ttk.Button(self, text="Browse…", command=self._browse).grid(
            row=0, column=2, **pad)

        # ── x-axis choice ──
        ttk.Label(self, text="X axis:").grid(row=1, column=0, sticky="w", **pad)
        self._xaxis_var = tk.StringVar(value="time")
        xf = ttk.Frame(self)
        xf.grid(row=1, column=1, sticky="w", **pad)
        ttk.Radiobutton(xf, text="Time (s)", variable=self._xaxis_var,
                        value="time").pack(side="left", padx=(0, 12))
        ttk.Radiobutton(xf, text="MU", variable=self._xaxis_var,
                        value="mu").pack(side="left")

        # ── channel checkboxes ──
        ttk.Label(self, text="Show:").grid(row=2, column=0, sticky="nw", **pad)
        chk_frame = ttk.Frame(self)
        chk_frame.grid(row=2, column=1, sticky="w", **pad)
        self._show_gantry   = tk.BooleanVar(value=True)
        self._show_coll     = tk.BooleanVar(value=True)
        self._show_dose     = tk.BooleanVar(value=True)
        self._show_mlc      = tk.BooleanVar(value=True)
        self._show_jaw      = tk.BooleanVar(value=False)
        ttk.Checkbutton(chk_frame, text="Gantry",    variable=self._show_gantry).pack(side="left", padx=4)
        ttk.Checkbutton(chk_frame, text="Collimator",variable=self._show_coll  ).pack(side="left", padx=4)
        ttk.Checkbutton(chk_frame, text="Dose rate", variable=self._show_dose  ).pack(side="left", padx=4)
        ttk.Checkbutton(chk_frame, text="MLC RMS",   variable=self._show_mlc   ).pack(side="left", padx=4)
        ttk.Checkbutton(chk_frame, text="Jaws",      variable=self._show_jaw   ).pack(side="left", padx=4)

        # ── plot button ──
        ttk.Button(self, text="Plot", command=self._plot).grid(
            row=3, column=0, columnspan=3, pady=6)

        # ── canvas area ──
        self._canvas_frame = ttk.Frame(self)
        self._canvas_frame.grid(row=4, column=0, columnspan=3, sticky="nsew", padx=4, pady=4)
        self.rowconfigure(4, weight=1)
        self.columnconfigure(1, weight=1)

        if TrajectoryLog is None:
            messagebox.showwarning(
                "pylinac not found",
                "Install pylinac to use this viewer:\n  pip install pylinac")

    def _browse(self):
        path = filedialog.askopenfilename(
            filetypes=[("Trajectory logs", "*.bin *.tlog"), ("All files", "*.*")])
        if path:
            self._path_var.set(path)
            self._load(path)

    def _load(self, path):
        if TrajectoryLog is None:
            return
        try:
            self._log = TrajectoryLog(path)
        except Exception as e:
            messagebox.showerror("Error loading log", str(e))
            self._log = None

    def _plot(self):
        if TrajectoryLog is None:
            messagebox.showerror("Missing dependency", "pylinac is required.")
            return
        path = self._path_var.get().strip()
        if not path:
            messagebox.showwarning("No file", "Select a trajectory log file first.")
            return
        if self._log is None:
            self._load(path)
        if self._log is None:
            return

        log = self._log

        # ── determine panels to draw ──
        panels = []
        if self._show_gantry.get():
            panels.append("gantry")
        if self._show_coll.get():
            panels.append("collimator")
        if self._show_dose.get():
            panels.append("dose_rate")
        if self._show_mlc.get():
            panels.append("mlc_rms")
        if self._show_jaw.get():
            panels.append("jaws")
        if not panels:
            messagebox.showwarning("Nothing selected", "Select at least one channel.")
            return

        # clear old canvas
        for w in self._canvas_frame.winfo_children():
            w.destroy()

        use_mu = self._xaxis_var.get() == "mu"
        xlabel = "MU" if use_mu else "Time (s)"

        fig, axes = plt.subplots(len(panels), 1, figsize=(12, 2.8 * len(panels)),
                                 sharex=True)
        if len(panels) == 1:
            axes = [axes]
        fig.subplots_adjust(hspace=0.35, left=0.1, right=0.97, top=0.95, bottom=0.07)

        # All time-series data lives on log.axis_data (subbeams are metadata only)
        ad = log.axis_data
        n = log.header.num_snapshots
        dt_s = log.header.sampling_interval / 1000.0  # ms → s

        if use_mu:
            x = np.asarray(ad.mu.actual)
            xlabel = "MU"
        else:
            x = np.arange(n) * dt_s
            xlabel = "Time (s)"

        C_ACT = "#1f77b4"
        C_EXP = "#aec7e8"

        for ax, panel in zip(axes, panels):
            if panel == "gantry":
                ax.plot(x, ad.gantry.actual,   color=C_ACT, lw=1.0, label="actual")
                ax.plot(x, ad.gantry.expected,  color=C_EXP, lw=0.8, ls="--", label="expected")
                ax.set_ylabel("Gantry (°)")

            elif panel == "collimator":
                ax.plot(x, ad.collimator.actual,  color=C_ACT, lw=1.0, label="actual")
                ax.plot(x, ad.collimator.expected, color=C_EXP, lw=0.8, ls="--", label="expected")
                ax.set_ylabel("Collimator (°)")

            elif panel == "dose_rate":
                ax.plot(x, ad.mu.actual, color=C_ACT, lw=1.0, label="cumulative MU")
                ax.set_ylabel("Cumulative MU")

            elif panel == "mlc_rms":
                leaves = ad.mlc.leaf_axes  # dict: leaf_num -> Axis
                if leaves:
                    diffs = np.stack(
                        [np.asarray(lf.actual) - np.asarray(lf.expected)
                         for lf in leaves.values()], axis=1)   # (n_samples, n_leaves)
                    rms = np.sqrt(np.mean(diffs**2, axis=1)) * 10  # cm → mm
                    ax.plot(x, rms, color=C_ACT, lw=1.0, label="RMS error")
                ax.set_ylabel("MLC RMS error (mm)")

            elif panel == "jaws":
                for jaw_attr, name in [("y1", "Y1"), ("y2", "Y2"),
                                       ("x1", "X1"), ("x2", "X2")]:
                    jaw = getattr(ad.jaws, jaw_attr, None)
                    if jaw is not None:
                        ax.plot(x, np.asarray(jaw.actual), lw=1.0, label=name)
                ax.set_ylabel("Jaw position (cm)")

        for ax, panel in zip(axes, panels):
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=7, loc="upper right")

        axes[-1].set_xlabel(xlabel)

        canvas = FigureCanvasTkAgg(fig, master=self._canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        toolbar = NavigationToolbar2Tk(canvas, self._canvas_frame)
        toolbar.update()
        plt.close(fig)


TLogViewer().mainloop()
