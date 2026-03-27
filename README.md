# Beam Data Comparison Tools

A collection of Python tools for radiation therapy beam data analysis, comparison, and QA.

**Download the latest version:** https://github.com/nelsknutson13/BeamDataComparisonTools

---

## Main GUI Tools

| Tool | Description |
|------|-------------|
| **ProfileCompare.py** | Compare measured vs. reference beam profiles with gamma index, DTA, and dose difference analysis. Generates PDF reports. |
| **PDDCompare.py** | Compare percent depth dose (PDD) curves with gamma index, DTA, and dose difference metrics. Generates PDF reports. |
| **TRS398Calculator.py** | TRS-398 absolute dose calculation with chamber-specific k_Q coefficients, volume-averaging corrections, and FFF beam support. |
| **PionCalculator.py** | Calculate pion correction factors and metrics for beam profiles, grouped by beam parameters. |
| **ProfileCorrectionFactorCalculator.py** | Calculate profile correction factors. |
| **OutputRoundRobinDataViewer.py** | Visualize round-robin QA consortium data as boxplots comparing institutions/systems by energy. |
| **DataDuplicateCleaner.py** | Detect and manage duplicate profile curves in grouped beam data using MD5 fingerprinting. |

## Data Readers / Converters

| Tool | Description |
|------|-------------|
| **IbaDataReader.py** | Convert IBA OPAX XML files to a standardized Excel format. |
| **PTWdatareader.py** | Convert PTW MCC files to Excel, with scan selection GUI. |
| **BlueEthosDataRipper.py** | Parse DICOM data from Blue Ethos (Ethos Enhanced RT Plan) downloads. |

## Batch / Automation Tools

| Tool | Description |
|------|-------------|
| **batch_profile_compare.py** | Batch process multiple profile comparison files and generate summary pass-rate statistics. |
| **batch_pdd_compare.py** | Batch process multiple PDD files across energy folders and generate summary statistics. |
| **scandadachecker.py** | Scan Excel files to inventory beam data (profiles, PDDs) by energy, SSD, field size, etc. Supports caching. |
| **make_summary_table.py** | Convert a PDD comparison summary CSV into a color-coded publication-ready table image. |
| **normalize_sheet_names.py** | Standardize Excel sheet names in consortium data files (e.g. "SN19" → "SN 19"). |

## Utility / Library Modules

| Module | Description |
|--------|-------------|
| **gamma.py** | Gamma index analysis for dose profile comparison. |
| **comp.py** | Dose difference and distance-to-agreement (DTA) calculation functions. |
| **center.py** | Profile centering based on 50% dose threshold. |

---

## Requirements

- Python 3.x
- Common packages: `pandas`, `numpy`, `matplotlib`, `tkinter`, `openpyxl`
- Some tools require: `pydicom`, `scipy`

## Getting Started

1. Download the ZIP from the link above or clone the repo
2. Install dependencies: `pip install pandas numpy matplotlib openpyxl pydicom scipy`
3. Run any GUI tool directly: `python ProfileCompare.py`
