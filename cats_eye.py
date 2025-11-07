#!/usr/bin/env python3
"""
Cat's Eye — Lightweight seeing estimator for astrophotography frames.
v0.1 — dual mode: standalone CLI or importable function for astro_cat
"""

import argparse, os, time, json, numpy as np
from pathlib import Path
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from datetime import datetime

# ------------------------------------------------------------
def load_fits_data(path):
    try:
        with fits.open(path, memmap=True) as hdul:
            data = hdul[0].data.astype(float)
            header = hdul[0].header
    except Exception:
        with fits.open(path, memmap=False) as hdul:
            data = hdul[0].data.astype(float)
            header = hdul[0].header
    return data, header

# ------------------------------------------------------------
def infer_pixel_scale(header):
    for key in ("PIXSCALE", "SCALE", "CDELT1"):
        if key in header:
            try:
                val = abs(float(header[key]))
                if "CDELT" in key and abs(val) < 1:
                    val *= 3600.0
                return val
            except Exception:
                pass
    pixsz_keys = ("XPIXSZ", "PIXSIZE1", "PIXSZ1", "PIXSIZE")
    foclen_keys = ("FOCALLEN", "FOCLEN", "FOCALLENGTH")
    pix_size = next((header[k] for k in pixsz_keys if k in header), None)
    focal_len = next((header[k] for k in foclen_keys if k in header), None)
    if pix_size and focal_len:
        try:
            return 206.265 * float(pix_size) / float(focal_len)
        except Exception:
            pass
    return None

# ------------------------------------------------------------
def measure_fwhm(fits_path, pixel_scale=1.0, min_stars=10,
                 threshold_sigma=5.0, downsample=1):
    """Return dict with filename and seeing_arcsec."""
    try:
        data, _ = load_fits_data(fits_path)
        if data.ndim > 2:
            data = data[0]
        if downsample > 1:
            data = data[::downsample, ::downsample]
            pixel_scale *= downsample

        mean, median, std = sigma_clipped_stats(data, sigma=3.0)
        fwhm_kernel = max(3.0 / downsample, 1.5)
        daofind = DAOStarFinder(fwhm=fwhm_kernel, threshold=threshold_sigma * std)
        sources = daofind(data - median)
        if sources is None or len(sources) < min_stars:
            return None

        if 'fwhm' in sources.colnames:
            fwhm_pix = float(np.median(sources['fwhm']))
        elif 'peak' in sources.colnames and 'flux' in sources.colnames:
            rel = np.sqrt(np.abs(sources['flux'] / np.maximum(sources['peak'], 1)))
            fwhm_pix = float(np.median(rel) / 2.35)
        else:
            fwhm_pix = daofind.fwhm

        return {
            "filename": os.path.basename(fits_path),
            "seeing_arcsec": round(fwhm_pix * pixel_scale, 3)
        }
    except Exception:
        return None

# ------------------------------------------------------------
def analyze_files(file_list, scale=None, downsample=1,
                  min_stars=10, threshold=5.0, max_procs=max(cpu_count()-1,1)):
    """
    Core API function.
    Accepts list of file paths or JSON object with filenames.
    Returns JSON-like dict with 'summary' and 'results'.
    """
    if isinstance(file_list, str):
        try:
            file_list = json.loads(file_list)
            if isinstance(file_list, dict) and "files" in file_list:
                file_list = file_list["files"]
        except json.JSONDecodeError:
            file_list = [file_list]
    files = [Path(f) for f in file_list if Path(f).exists()]
    if not files:
        return {"error": "No valid FITS files provided"}

    # Infer pixel scale from first file if needed
    try:
        _, hdr = load_fits_data(files[0])
        scale = scale or infer_pixel_scale(hdr) or 1.0
    except Exception:
        scale = scale or 1.0

    start = time.time()
    per_file = []

    with ProcessPoolExecutor(max_workers=max_procs) as executor:
        futures = {executor.submit(measure_fwhm, f, scale, min_stars,
                                   threshold, downsample): f for f in files}
        for future in as_completed(futures):
            val = future.result()
            if val:
                per_file.append(val)

    elapsed = time.time() - start
    if not per_file:
        return {"error": "No valid measurements"}

    arr = np.array([r["seeing_arcsec"] for r in per_file])
    summary = {
        "files_processed": len(per_file),
        "total_files": len(files),
        "mean_seeing": round(arr.mean(), 3),
        "min_seeing": round(arr.min(), 3),
        "max_seeing": round(arr.max(), 3),
        "std_seeing": round(arr.std(), 3),
        "pixel_scale": round(scale, 4),
        "downsample": downsample,
        "runtime_sec": round(elapsed, 2),
        "avg_time_per_file": round(elapsed / len(files), 2),
        "timestamp": datetime.utcnow().isoformat() + "Z",
    }

    return {"summary": summary, "results": per_file}

# ------------------------------------------------------------
def gather_files(input_path):
    p = Path(input_path)
    if p.is_file():
        return [p]
    elif p.is_dir():
        return list(p.rglob("*.fits")) or list(p.rglob("*.fit"))
    else:
        raise FileNotFoundError(f"No such file or directory: {input_path}")

# ------------------------------------------------------------
def main():
    BANNER = r"""
█▓▒░  C A T ' S   E Y E  ░▒▓█
   🐈  Seeing Estimator for Astro_Cat
--------------------------------------
"""
    print(BANNER)
    parser = argparse.ArgumentParser(description="Measure seeing (FWHM) for one file or all FITS in a folder.")
    parser.add_argument("path", help="Path to FITS file or folder")
    parser.add_argument("--scale", type=float, default=None)
    parser.add_argument("--min-stars", type=int, default=10)
    parser.add_argument("--threshold", type=float, default=5.0)
    parser.add_argument("--downsample", type=int, default=1)
    parser.add_argument("--max-procs", type=int, default=max(cpu_count()-2,1))
    parser.add_argument("--json", action="store_true", help="Emit JSON summary to cats_eye_results.json")
    args = parser.parse_args()

    files = gather_files(args.path)
    print(f"📂 Found {len(files)} file(s) to process using up to {args.max_procs} cores.\n")

    results = analyze_files(
        file_list=files,
        scale=args.scale,
        downsample=args.downsample,
        min_stars=args.min_stars,
        threshold=args.threshold,
        max_procs=args.max_procs
    )

    if "error" in results:
        print(f"⚠️  {results['error']}")
        return

    summary = results["summary"]
    print(f"🔎 Pixel scale: {summary['pixel_scale']:.3f}\"/px\n")
    print("📊 Summary:")
    print(f"   Files processed : {summary['files_processed']} / {summary['total_files']}")
    print(f"   Mean seeing     : {summary['mean_seeing']:.2f}\"")
    print(f"   Min / Max       : {summary['min_seeing']:.2f}\" / {summary['max_seeing']:.2f}\"")
    print(f"   Std. deviation  : {summary['std_seeing']:.2f}\"")
    print("\n⏱️  Timing:")
    print(f"   Total time      : {summary['runtime_sec']:.1f} s")
    print(f"   Avg per file    : {summary['avg_time_per_file']:.2f} s/file")
    print(f"   Projected (100) : {(summary['avg_time_per_file']*100):.1f} s for 100 files")

    if args.json:
        out_path = Path(args.path) / "cats_eye_results.json"
        with open(out_path, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\n💾 JSON summary written to: {out_path}")

if __name__ == "__main__":
    main()
