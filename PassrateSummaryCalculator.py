# -*- coding: utf-8 -*-
"""
Created on Fri Aug  8 14:51:19 2025

@author: nknutson
"""

import re
import pandas as pd
from collections import defaultdict
import numpy as np
from tabulate import tabulate

LOG = r"""Starting comparison using analysis: Gamma
DD: 0.07, DTA: 0.5
HN | Site 1 | AntPost → Points with Γ > 1: 0/596   Pass Rate: 100.0%
HN | Site 1 | RightLeft → Points with Γ > 1: 22/929   Pass Rate: 97.6%
HN | Site 1 | SupInf → Points with Γ > 1: 0/783   Pass Rate: 100.0%
 HN | Site 2 | AntPost, Skipping — data missing in one of the sheets:
 HN | Site 2 | RightLeft, Skipping — data missing in one of the sheets:
 HN | Site 2 | SupInf, Skipping — data missing in one of the sheets:
 HN | Site 3 | AntPost, Skipping — data missing in one of the sheets:
 HN | Site 3 | RightLeft, Skipping — data missing in one of the sheets:
 HN | Site 3 | SupInf, Skipping — data missing in one of the sheets:
HN | Site 5 | AntPost → Points with Γ > 1: 30/658   Pass Rate: 95.4%
HN | Site 5 | RightLeft → Points with Γ > 1: 0/939   Pass Rate: 100.0%
HN | Site 5 | SupInf → Points with Γ > 1: 0/844   Pass Rate: 100.0%
 Head | Site 1 | AntPost, Skipping — data missing in one of the sheets:
 Head | Site 1 | RightLeft, Skipping — data missing in one of the sheets:
 Head | Site 1 | SupInf, Skipping — data missing in one of the sheets:
Head | Site 2 | AntPost → Points with Γ > 1: 0/563   Pass Rate: 100.0%
Head | Site 2 | RightLeft → Points with Γ > 1: 0/563   Pass Rate: 100.0%
Head | Site 2 | SupInf → Points with Γ > 1: 238/373   Pass Rate: 36.2%
 Head | Site 3 | AntPost, Skipping — data missing in one of the sheets:
 Head | Site 3 | RightLeft, Skipping — data missing in one of the sheets:
 Head | Site 3 | SupInf, Skipping — data missing in one of the sheets:
Head | Site 5 | AntPost → Points with Γ > 1: 0/578   Pass Rate: 100.0%
Head | Site 5 | RightLeft → Points with Γ > 1: 0/586   Pass Rate: 100.0%
Head | Site 5 | SupInf → Points with Γ > 1: 0/550   Pass Rate: 100.0%
Lung | Site 1 | AntPost → Points with Γ > 1: 8/891   Pass Rate: 99.1%
Lung | Site 1 | RightLeft → Points with Γ > 1: 0/1031   Pass Rate: 100.0%
Lung | Site 1 | SupInf → Points with Γ > 1: 68/1100   Pass Rate: 93.8%
Lung | Site 2 | AntPost → Points with Γ > 1: 5/1006   Pass Rate: 99.5%
Lung | Site 2 | RightLeft → Points with Γ > 1: 0/1111   Pass Rate: 100.0%
Lung | Site 2 | SupInf → Points with Γ > 1: 0/1140   Pass Rate: 100.0%
Lung | Site 3 | AntPost → Points with Γ > 1: 0/912   Pass Rate: 100.0%
Lung | Site 3 | RightLeft → Points with Γ > 1: 322/1116   Pass Rate: 71.1%
Lung | Site 3 | SupInf → Points with Γ > 1: 272/1172   Pass Rate: 76.8%
 Lung | Site 5 | AntPost, Skipping — data missing in one of the sheets:
 Lung | Site 5 | RightLeft, Skipping — data missing in one of the sheets:
 Lung | Site 5 | SupInf, Skipping — data missing in one of the sheets:
 Lung 10F | Site 1 | AntPost, Skipping — data missing in one of the sheets:
 Lung 10F | Site 1 | RightLeft, Skipping — data missing in one of the sheets:
 Lung 10F | Site 1 | SupInf, Skipping — data missing in one of the sheets:
Lung 10F | Site 2 | AntPost → Points with Γ > 1: 0/861   Pass Rate: 100.0%
Lung 10F | Site 2 | RightLeft → Points with Γ > 1: 0/1057   Pass Rate: 100.0%
 Lung 10F | Site 2 | SupInf, Skipping — data missing in one of the sheets:
 Lung 10F | Site 3 | AntPost, Skipping — data missing in one of the sheets:
 Lung 10F | Site 3 | RightLeft, Skipping — data missing in one of the sheets:
 Lung 10F | Site 3 | SupInf, Skipping — data missing in one of the sheets:
 Lung 10F | Site 5 | AntPost, Skipping — data missing in one of the sheets:
 Lung 10F | Site 5 | RightLeft, Skipping — data missing in one of the sheets:
 Lung 10F | Site 5 | SupInf, Skipping — data missing in one of the sheets:

"""

def _norm_prof(s):
    s = str(s).strip().lower()
    if s in ("supinf","sup-inf","si","sup_inf","sup–inf"): return "SupInf"
    if s in ("rightleft","right-left","rl","lr","left-right","left_right"): return "RightLeft"
    if s in ("antpost","anterior-posterior","ap","ant_post","ant-post"): return "AntPost"
    return s

# Match result lines like:
# HN | Site 1 | AntPost → Points with Γ > 1: 0/596   Pass Rate: 100.0%
pass_pat = re.compile(
    r'^\s*(?P<phantom>[^|]+)\|\s*(?P<site>[^|]+)\|\s*(?P<profile>\S+).*?'
    r'Points with (?:Γ|Gamma)\s*>\s*1:\s*(?P<fail>\d+)\s*/\s*(?P<total>\d+)',
    re.IGNORECASE
)

# Match skip lines in either order:
#  "Skipping — data missing...: HN | Site 2 | AntPost"
#  "HN | Site 2 | AntPost, Skipping — data missing..."
skip_pat_1 = re.compile(
    r'Skipping\s+—\s+data missing.*?:\s*(?P<phantom>[^|]+)\|\s*(?P<site>[^|]+)\|\s*(?P<profile>\S+)', re.IGNORECASE)
skip_pat_2 = re.compile(
    r'^\s*(?P<phantom>[^|]+)\|\s*(?P<site>[^|]+)\|\s*(?P<profile>\S+)\s*,\s*Skipping\s+—\s+data missing', re.IGNORECASE)

# De-dupe: keep the latest or the one with the largest total
best = {}  # (phantom, site, profile) -> (passed, total)
skipped = set()

for line in LOG.splitlines():
    m = pass_pat.search(line)
    if m:
        p = m.group('phantom').strip()
        s = m.group('site').strip()
        prof = _norm_prof(m.group('profile'))
        fail = int(m.group('fail')); total = int(m.group('total'))
        passed = total - fail
        key = (p, s, prof)
        if key in best and total <= best[key][1]:
            continue
        best[key] = (passed, total)
        continue
    # track skips (not used in totals, but handy if you want a skip report)
    m1 = skip_pat_1.search(line) or skip_pat_2.search(line)
    if m1:
        p = m1.group('phantom').strip()
        s = m1.group('site').strip()
        prof = _norm_prof(m1.group('profile'))
        skipped.add((p, s, prof))

# Build (Phantom|Site) table: Overall, Coronal (SI+RL), Sagittal (SI+AP)
overall_pass = defaultdict(int); overall_total = defaultdict(int)
prof_pass = defaultdict(int);    prof_total = defaultdict(int)

for (p, s, prof), (ps, tot) in best.items():
    overall_pass[(p, s)] += ps
    overall_total[(p, s)] += tot
    prof_pass[(p, s, prof)] += ps
    prof_total[(p, s, prof)] += tot

rows = []
for (p, s) in sorted(overall_total):
    o_ps, o_tt = overall_pass[(p, s)], overall_total[(p, s)]
    c_ps = prof_pass.get((p, s, "SupInf"), 0) + prof_pass.get((p, s, "RightLeft"), 0)
    c_tt = prof_total.get((p, s, "SupInf"), 0) + prof_total.get((p, s, "RightLeft"), 0)
    g_ps = prof_pass.get((p, s, "SupInf"), 0) + prof_pass.get((p, s, "AntPost"), 0)
    g_tt = prof_total.get((p, s, "SupInf"), 0) + prof_total.get((p, s, "AntPost"), 0)

    rows.append({
        "Phantom": p, "Site": s,
        "Overall_pass": o_ps, "Overall_total": o_tt, "Overall_%": round(100*o_ps/o_tt, 1) if o_tt else None,
        "Coronal_pass": c_ps, "Coronal_total": c_tt, "Coronal_%": round(100*c_ps/c_tt, 1) if c_tt else None,
        "Sagittal_pass": g_ps, "Sagittal_total": g_tt, "Sagittal_%": round(100*g_ps/g_tt, 1) if g_tt else None,
    })

df = pd.DataFrame(rows, columns=[
    "Phantom","Site",
    "Overall_pass","Overall_total","Overall_%",
    "Coronal_pass","Coronal_total","Coronal_%",
    "Sagittal_pass","Sagittal_total","Sagittal_%"
])

#print(df.to_string(index=False))
print(tabulate(df, headers="keys", tablefmt="github", showindex=False, floatfmt=".1f"))

df.to_excel("gamma_summary.xlsx", index=False)  # optional


# assume df is your table
df = df.copy()

# add fail counts
df["Overall_fail"]  = df["Overall_total"]  - df["Overall_pass"]
df["Coronal_fail"]  = df["Coronal_total"]  - df["Coronal_pass"]
df["Sagittal_fail"] = df["Sagittal_total"] - df["Sagittal_pass"]

# sort nicely by Phantom, then numeric Site
df["Site_no"] = df["Site"].str.extract(r'(\d+)').astype(float)
df = df.sort_values(["Phantom", "Site_no"], kind="mergesort").drop(columns="Site_no")

# flag anything below threshold
THRESH = 95.0
for col in ["Overall_%", "Coronal_%", "Sagittal_%"]:
    df[col + "_flag"] = df[col] < THRESH

print(df.to_string(index=False))

from tabulate import tabulate

# per-Phantom weighted summary
grp = (df.groupby("Phantom", as_index=False)
         .agg(Overall_pass=("Overall_pass","sum"),
              Overall_total=("Overall_total","sum"),
              Coronal_pass=("Coronal_pass","sum"),
              Coronal_total=("Coronal_total","sum"),
              Sagittal_pass=("Sagittal_pass","sum"),
              Sagittal_total=("Sagittal_total","sum")))

for a in ["Overall","Coronal","Sagittal"]:
    tot = grp[a + "_total"].replace(0, pd.NA)
    grp[a + "_%"] = (100.0 * grp[a + "_pass"] / tot).round(1)

cols = ["Phantom",
        "Overall_pass","Overall_total","Overall_%",
        "Coronal_pass","Coronal_total","Coronal_%",
        "Sagittal_pass","Sagittal_total","Sagittal_%"]

print("\nPer-Phantom summary (weighted by points):")
print(tabulate(grp[cols], headers="keys", tablefmt="github", showindex=False))
