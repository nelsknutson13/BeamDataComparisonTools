"""
normalize_sheet_names.py

Normalizes Excel sheet names in the consortium PDD/Profile data directory.

Rules:
  - "SNxx"      → "SN xx"      (e.g. SN19 → SN 19)
  - "TPSSNxx"   → "TPS SN xx"  (e.g. TPSSN20 → TPS SN 20, covers missing spaces)
  - "TPS SNxx"  → "TPS SN xx"  (e.g. TPS SN20 → TPS SN 20)

Sheets that already match the target format or do not match any known pattern
are left untouched (including SN20_MD and similar outliers).

Run in dry-run mode first (default) to preview changes, then pass --apply to
actually rename the sheets.

Usage:
    python normalize_sheet_names.py              # dry run (preview only)
    python normalize_sheet_names.py --apply      # write changes to files
"""

import re
import sys
from pathlib import Path

try:
    import openpyxl
except ImportError:
    sys.exit("openpyxl is required: pip install openpyxl")

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BASE_PATH = Path(
    r"C:\Users\nknutson\OneDrive - Washington University in St. Louis"
    r"\NGDS QA Consortium\Combined Consortium Data\Processed Combined Data"
)

# Regex patterns and their normalizers
# Each entry: (compiled_pattern, replacement_func)
# Replacement func receives the regex match and returns the normalized name.

def _normalize(name: str):
    """
    Return normalized sheet name, or None if no change needed.
    Handles:
      SN<digits>          → SN <digits>
      TPS SN<digits>      → TPS SN <digits>
      TPSSN<digits>       → TPS SN <digits>    (both spaces missing)
      TPS  SN <digits>    → TPS SN <digits>    (extra spaces)
    Returns None if the name already looks correct or is unrecognized.
    """
    # Normalize internal whitespace for comparison
    s = re.sub(r'\s+', ' ', name.strip())

    # Pattern: optional "TPS " prefix, then "SN", optional space, digits,
    # and an optional suffix like "_MD" (microdiamond detector dataset).
    m = re.fullmatch(
        r'(TPS\s+)?SN\s*(\d+)(_\w+)?',
        s,
        re.IGNORECASE
    )
    if not m:
        return None  # unrecognized pattern — leave unchanged

    tps_prefix = m.group(1)
    digits = m.group(2)
    suffix = m.group(3) or ""   # e.g. "_MD", or empty string

    if tps_prefix is not None:
        target = f"TPS SN {digits}{suffix}"
    else:
        target = f"SN {digits}{suffix}"

    return target if target != name else None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(apply: bool = False):
    xlsx_files = sorted(BASE_PATH.glob("*.xlsx"))
    if not xlsx_files:
        print(f"No .xlsx files found in:\n  {BASE_PATH}")
        return

    any_change = False

    for xlsx_path in xlsx_files:
        wb = openpyxl.load_workbook(xlsx_path)
        renames = {}  # old_name → new_name

        for sheet_name in wb.sheetnames:
            normalized = _normalize(sheet_name)
            if normalized is not None and normalized != sheet_name:
                renames[sheet_name] = normalized

        if not renames:
            print(f"  [ok]  {xlsx_path.name}")
            continue

        any_change = True
        print(f"  [{'APPLY' if apply else 'DRY RUN'}]  {xlsx_path.name}")
        for old, new in renames.items():
            print(f"         '{old}'  ->  '{new}'")

        if apply:
            for old, new in renames.items():
                wb[old].title = new
            wb.save(xlsx_path)
            print(f"         Saved.")

    if not any_change:
        print("\nAll sheet names are already normalized. Nothing to do.")
    elif not apply:
        print("\nDry run complete. Run with --apply to write changes.")
    else:
        print("\nAll renames applied.")


if __name__ == "__main__":
    apply_mode = "--apply" in sys.argv
    if apply_mode:
        print("=== APPLY MODE — files will be modified ===\n")
    else:
        print("=== DRY RUN — no files will be modified ===\n")
    run(apply=apply_mode)
