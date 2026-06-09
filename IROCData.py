import tkinter as tk
from tkinter import filedialog
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def load_and_plot():
    file_path = filedialog.askopenfilename(filetypes=[("Excel files", "*.xlsx *.xls")])
    if not file_path:
        return

    try:
        # Load the Excel file and read the 'TLD Data' sheet
        df = pd.read_excel(file_path, sheet_name='TLD Data')

        # Drop rows with missing data in key columns
        df_clean = df.dropna(subset=['%RX Difference', 'Body Site', 'Report'])

        # Simplify 'Report' to just Site# before underscore
        df_clean['Report'] = df_clean['Report'].apply(lambda x: x.split('_')[0] if isinstance(x, str) else x)

        # Create 'All' group for overall stats
        df_all = df_clean.copy()
        df_all['Body Site'] = 'All'
        df_clean = pd.concat([df_clean, df_all], ignore_index=True)

        # Build unique site list and assign fixed markers
        unique_sites = sorted(df_clean['Report'].unique())
        markers = ['o', 's', '^', 'D', 'P', 'X', '*', 'v', '<', '>']
        marker_map = {site: markers[i % len(markers)] for i, site in enumerate(unique_sites)}

        # Setup plot
        plt.figure(figsize=(10, 6))
        sns.set(style="whitegrid", palette="colorblind")
        flier_style = dict(marker='o', markerfacecolor='none', markersize=5, linestyle='none', color='gray')
        ax = sns.boxplot(data=df_clean, x='Body Site', y='%RX Difference', whis=1.5, flierprops=flier_style)

        # Plot data by phantom type with consistent site markers
        for i, body_site in enumerate(df_clean['Body Site'].unique()):
            subset = df_clean[df_clean['Body Site'] == body_site]

            # Skip plotting markers for 'All' group
            #if body_site == 'All':
             #   continue

            for site in unique_sites:
                site_data = subset[subset['Report'] == site]
                if not site_data.empty:
                    x_jitter = np.random.normal(loc=i, scale=0.08, size=len(site_data))
                    plt.scatter(x_jitter, site_data['%RX Difference'],
                                marker=marker_map[site],
                                edgecolor='black', facecolor='none', s=70, alpha=0.8)

        plt.axhline(0, color='gray', linestyle='--')
        plt.ylabel('Difference (MC vs IROC)[%RX]')
        plt.xlabel('Phantom Type')

        # Build legend with consistent site-marker mapping
        legend_handles = []
        for site in unique_sites:
            handle = plt.Line2D([], [], marker=marker_map[site], linestyle='None',
                                markersize=7, label=site,
                                markerfacecolor='none', markeredgecolor='black')
            legend_handles.append(handle)

        plt.legend(handles=legend_handles,
                   title="Site:", loc='upper right',
                   fontsize='small', title_fontsize='small', frameon=True)

        plt.tight_layout()
        plt.show()

    except Exception as e:
        print(f"Error: {e}")

# Simple GUI
root = tk.Tk()
root.title("IROC TLD Plotter")
root.geometry("300x100")

btn_load = tk.Button(root, text="Load Excel and Plot", command=load_and_plot)
btn_load.pack(pady=30)

root.mainloop()
