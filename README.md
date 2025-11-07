# 🐈 Cat’s Eye  
**Lightweight Seeing Estimator for Astrophotography**

Cat’s Eye is a standalone and embeddable Python tool for measuring **astronomical seeing** (FWHM of stellar profiles) from raw FITS images.  
Originally developed as a companion to **[Astro_Cat](https://github.com/)**, it can be used on its own or integrated as a Python module.

---

## 🚀 Features

- 🔭 **Automated seeing estimation** from FITS/XISF images using `photutils` and `astropy`
- ⚙️ **Pixel scale detection** from FITS header metadata  
- 💻 **Parallel processing** using multiple CPU cores for rapid batch analysis  
- 📉 **Optional downsampling** for faster approximate results  
- 🧩 **JSON output** — ready for ingestion into Astro_Cat databases or dashboards  
- 🐍 **Dual mode:** command-line tool **or** importable Python module  

---

## 📦 Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/<your-username>/cats_eye.git
cd cats_eye
python -m venv venv
source venv/bin/activate    # or .\venv\Scripts\activate on Windows
pip install -r requirements.txt
```

---

## 🧭 Command Line Usage

Analyze one file or an entire folder of FITS images:

```bash
python cats_eye.py /path/to/folder
```

Optional arguments:

| Option | Description |
|--------|-------------|
| `--scale` | Override pixel scale (arcsec/pixel) |
| `--downsample` | Downsample factor (e.g. 2 = 2×2 binning) |
| `--min-stars` | Minimum stars required for valid result |
| `--threshold` | Detection threshold in sigma |
| `--max-procs` | Max CPU cores to use (default = all - 1) |
| `--json` | Output results to `cats_eye_results.json` |

Example:

```bash
python cats_eye.py ~/Astro/Images/NGC1333/ --downsample 2 --json
```

**Output (console):**
```
📂 Found 33 file(s) to process using up to 6 cores.
🔎 Pixel scale: 0.778"/px

📊 Summary:
   Files processed : 33 / 33
   Mean seeing     : 1.22"
   Min / Max       : 1.11" / 1.35"
   Std. deviation  : 0.06"

💾 JSON summary written to: cats_eye_results.json
```

**Example JSON Output:**
```json
{
  "summary": {
    "files_processed": 33,
    "mean_seeing": 1.22,
    "min_seeing": 1.11,
    "max_seeing": 1.35,
    "std_seeing": 0.06,
    "pixel_scale": 0.778,
    "runtime_sec": 66.4
  },
  "results": [
    {"filename": "NGC1333_LIGHT_0001.fits", "seeing_arcsec": 1.24},
    {"filename": "NGC1333_LIGHT_0002.fits", "seeing_arcsec": 1.27}
  ]
}
```

---

## 🧠 Python API Usage

You can also call Cat’s Eye directly from Python (e.g. inside **Astro_Cat**):

```python
from cats_eye import analyze_files

files = [
    "/data/NGC1333_LIGHT_0001.fits",
    "/data/NGC1333_LIGHT_0002.fits"
]

results = analyze_files(
    file_list=files,
    scale=0.778,
    downsample=2,
    max_procs=6
)

print("Mean seeing:", results["summary"]["mean_seeing"], "arcsec")
```

---

## ⚖️ License

MIT License © 2025 Corey Smart

---

## 🧭 Credits

- **Author:** Corey Smart  
- **Libraries:** [Astropy](https://www.astropy.org/), [Photutils](https://photutils.readthedocs.io/), [TQDM](https://github.com/tqdm/tqdm)  
- **Part of:** the [Astro_Cat](https://github.com/) ecosystem for astrophotography data automation  
