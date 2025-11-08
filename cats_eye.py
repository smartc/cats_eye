#!/usr/bin/env python3
"""
Cat's Eye v0.2 — Multi-metric image profiler for Astro_Cat
Author: Corey Smart & ChatGPT
"""

import argparse, os, time, json, numpy as np
from pathlib import Path
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from photutils.segmentation import SourceCatalog
from photutils.aperture import CircularAperture


# ---------------------------------------------------------------------
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


def infer_pixel_scale(header):
    for key in ("PIXSCALE", "SCALE", "CDELT1"):
        if key in header:
            val = abs(float(header[key]))
            if "CDELT" in key and abs(val) < 1:
                val *= 3600.0
            return val
    for k1 in ("XPIXSZ", "PIXSIZE1", "PIXSZ1", "PIXSIZE"):
        for k2 in ("FOCALLEN", "FOCLEN", "FOCALLENGTH"):
            if k1 in header and k2 in header:
                try:
                    return 206.265 * float(header[k1]) / float(header[k2])
                except:
                    pass
    return 1.0


# ---------------------------------------------------------------------
def compute_star_moments(data, x, y, radius=10):
    """
    Compute second moments (covariance matrix) for a star.
    Returns (moments_cxx, moments_cyy, moments_cxy) or None if failed.
    """
    try:
        # Extract cutout around star
        x_int, y_int = int(np.round(x)), int(np.round(y))
        ymin = max(0, y_int - radius)
        ymax = min(data.shape[0], y_int + radius + 1)
        xmin = max(0, x_int - radius)
        xmax = min(data.shape[1], x_int + radius + 1)

        cutout = data[ymin:ymax, xmin:xmax]
        if cutout.size == 0:
            return None

        # Create coordinate grids relative to centroid
        y_grid, x_grid = np.ogrid[ymin:ymax, xmin:xmax]
        y_grid = y_grid.astype(float) - y
        x_grid = x_grid.astype(float) - x

        # Subtract local background
        cutout = cutout - np.median(cutout)
        cutout = np.maximum(cutout, 0)  # Clip negative values

        total_flux = np.sum(cutout)
        if total_flux <= 0:
            return None

        # Compute normalized second moments
        moments_xx = np.sum(cutout * x_grid * x_grid) / total_flux
        moments_yy = np.sum(cutout * y_grid * y_grid) / total_flux
        moments_xy = np.sum(cutout * x_grid * y_grid) / total_flux

        return moments_xx, moments_yy, moments_xy
    except:
        return None


# ---------------------------------------------------------------------
def measure_metrics(fits_path, pixel_scale=1.0, downsample=1,
                    min_stars=10, threshold_sigma=5.0):
    try:
        data, _ = load_fits_data(fits_path)
        if data.ndim > 2:
            data = data[0]

        # simple decimation downsample (fast and effective)
        if downsample > 1:
            data = data[::downsample, ::downsample]
            # Note: pixel_scale stays at original value because the FWHM calculation
            # from flux/peak ratio doesn't scale linearly with downsampling

        mean, median, std = sigma_clipped_stats(data, sigma=3.0)
        fwhm_kernel = max(3.0 / downsample, 1.5)

        daofind = DAOStarFinder(fwhm=fwhm_kernel, threshold=threshold_sigma * std)
        sources = daofind(data)
        if sources is None or len(sources) < min_stars:
            return None

        # --- Compute second moments for shape measurements ---
        moments_list = []
        for source in sources:
            moments = compute_star_moments(data, source['xcentroid'], source['ycentroid'],
                                          radius=int(3 * fwhm_kernel))
            moments_list.append(moments if moments else (np.nan, np.nan, np.nan))

        # Add moments to sources table
        moments_arr = np.array(moments_list)
        sources['moments_cxx'] = moments_arr[:, 0]
        sources['moments_cyy'] = moments_arr[:, 1]
        sources['moments_cxy'] = moments_arr[:, 2]

        # --- FWHM computation (restored flux/peak logic) ---
        if 'peak' in sources.colnames and 'flux' in sources.colnames:
            rel = np.sqrt(np.abs(sources['flux'] / np.maximum(sources['peak'], 1)))
            # For 2D Gaussian: flux = 2π * peak * σ_x * σ_y
            # So: sqrt(flux/peak) = sqrt(2π * σ_x * σ_y)
            # FWHM = 2.355 * sqrt(σ_x * σ_y) = 2.355 / sqrt(2π) * sqrt(flux/peak)
            fwhm_pix = np.median(rel) * (2.355 / np.sqrt(2 * np.pi))
        elif 'fwhm' in sources.colnames:
            fwhm_pix = np.median(sources['fwhm'])
        elif 'sigma_x' in sources.colnames and 'sigma_y' in sources.colnames:
            fwhm_pix = np.median(2.355 * np.sqrt(sources['sigma_x'] * sources['sigma_y']))
        else:
            fwhm_pix = fwhm_kernel

        seeing_arcsec = fwhm_pix * pixel_scale
        hfr_pix = 0.5 * fwhm_pix

        # --- Eccentricity / elongation ---
        if all(c in sources.colnames for c in ("a", "b", "theta")):
            ecc = np.sqrt(1 - (sources["b"] / sources["a"]) ** 2)
            ang = np.degrees(sources["theta"])
            eccentricity_med = np.median(ecc)
            elong_angle = np.median(ang)
        elif all(c in sources.colnames for c in ("moments_cxx", "moments_cyy", "moments_cxy")):
            cxx, cyy, cxy = sources["moments_cxx"], sources["moments_cyy"], sources["moments_cxy"]

            # Filter out NaN values and compute shape parameters
            valid_mask = np.isfinite(cxx) & np.isfinite(cyy) & np.isfinite(cxy)
            if np.sum(valid_mask) > 0:
                cxx, cyy, cxy = cxx[valid_mask], cyy[valid_mask], cxy[valid_mask]
                term = np.sqrt((cxx - cyy) ** 2 + 4 * cxy ** 2)
                a = np.sqrt(0.5 * (cxx + cyy + term))
                b = np.sqrt(0.5 * (cxx + cyy - term))

                # Avoid division by zero and invalid values
                valid_ab = (a > 0) & (b > 0) & (b <= a)
                if np.sum(valid_ab) > 0:
                    ecc = np.sqrt(1 - (b[valid_ab] / a[valid_ab]) ** 2)
                    theta = 0.5 * np.degrees(np.arctan2(2 * cxy[valid_ab], cxx[valid_ab] - cyy[valid_ab]))
                    eccentricity_med = float(np.median(ecc))
                    elong_angle = float(np.median(theta))
                else:
                    eccentricity_med, elong_angle = 0.0, 0.0
            else:
                eccentricity_med, elong_angle = 0.0, 0.0
        else:
            eccentricity_med, elong_angle = 0.0, 0.0

        star_count = len(sources)
        snr_proxy = float(np.median((sources["peak"] - median) / std))
        weight = (star_count / (1 + std)) / seeing_arcsec

        return {
            "filename": os.path.basename(fits_path),
            "fwhm_px": round(float(fwhm_pix), 3),
            "hfr_px": round(float(hfr_pix), 3),
            "seeing_arcsec": round(float(seeing_arcsec), 3),
            "eccentricity_med": round(float(eccentricity_med), 3),
            "elongation_angle_deg": round(float(elong_angle), 1),
            "star_count": int(star_count),
            "bkg_median_adu": round(float(median), 2),
            "bkg_sigma_adu": round(float(std), 2),
            "snr_proxy": round(float(snr_proxy), 2),
            "weight": round(float(weight), 5),
        }

    except Exception as e:
        print(f"⚠️  Error on {os.path.basename(fits_path)}: {e}")
        return None


# ---------------------------------------------------------------------
def analyze_files(file_list, scale=None, downsample=1,
                  min_stars=10, threshold=5.0,
                  max_procs=max(cpu_count() - 2, 1)):
    if isinstance(file_list, str):
        try:
            file_list = json.loads(file_list)
            if isinstance(file_list, dict) and "files" in file_list:
                file_list = file_list["files"]
        except:
            file_list = [file_list]
    files = [Path(f) for f in file_list if Path(f).exists()]
    if not files:
        return {"error": "No valid FITS files"}

    _, hdr = load_fits_data(files[0])
    scale = scale or infer_pixel_scale(hdr)

    start = time.time()
    results = []

    with ProcessPoolExecutor(max_workers=max_procs) as ex:
        futs = {ex.submit(measure_metrics, f, scale, downsample, min_stars, threshold): f for f in files}
        for fut in tqdm(as_completed(futs), total=len(futs), desc="Analyzing"):
            r = fut.result()
            if r:
                results.append(r)

    elapsed = time.time() - start
    if not results:
        return {"error": "No valid results"}

    # normalize weights (NumPy 2.0 safe)
    w = np.array([r["weight"] for r in results])
    rng = np.ptp(w)
    wnorm = (w - np.min(w)) / (rng if rng > 0 else 1)
    for r, v in zip(results, wnorm):
        r["weight"] = round(float(v), 3)

    arr = np.array([r["seeing_arcsec"] for r in results])
    summary = {
        "files_processed": len(results),
        "mean_seeing": round(arr.mean(), 3),
        "min_seeing": round(arr.min(), 3),
        "max_seeing": round(arr.max(), 3),
        "std_seeing": round(arr.std(), 3),
        "pixel_scale": round(scale, 4),
        "downsample": downsample,
        "runtime_sec": round(elapsed, 1),
    }
    return {"summary": summary, "results": results}


# ---------------------------------------------------------------------
def gather_files(path):
    p = Path(path)
    return [p] if p.is_file() else list(p.rglob("*.fits"))


# ---------------------------------------------------------------------
def main():
    print("█▓▒░  C A T ' S   E Y E  ░▒▓█")
    print("   🐈  Multi-metric image profiler\n--------------------------------------")

    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--scale", type=float)
    ap.add_argument("--downsample", type=int, default=1)
    ap.add_argument("--min-stars", type=int, default=10)
    ap.add_argument("--threshold", type=float, default=5.0)
    ap.add_argument("--max-procs", type=int, default=max(cpu_count() - 2, 1))
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()

    files = gather_files(a.path)
    print(f"📂 Found {len(files)} file(s) to process using up to {a.max_procs} cores.\n")

    res = analyze_files(files, a.scale, a.downsample, a.min_stars, a.threshold, a.max_procs)
    if "error" in res:
        print("⚠️", res["error"])
        return

    s = res["summary"]
    print(f"🔎 Pixel scale: {s['pixel_scale']:.3f}\"/px\n📊 Summary:")
    print(f"   Files processed : {s['files_processed']}")
    print(f"   Mean seeing     : {s['mean_seeing']:.2f}\"")
    print(f"   Min / Max       : {s['min_seeing']:.2f}\" / {s['max_seeing']:.2f}\"")
    print(f"   Std. deviation  : {s['std_seeing']:.2f}\"")
    print(f"   Runtime         : {s['runtime_sec']:.1f} s\n")

    if a.json:
        out = Path(a.path) / "cats_eye_results.json"
        with open(out, "w") as f:
            json.dump(res, f, indent=2)
        print(f"💾 JSON summary written to: {out}")


if __name__ == "__main__":
    main()
