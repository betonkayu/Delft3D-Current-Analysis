"""
© Aimar Rendra P.A. Teknik Kelautan — Universitas Hasanuddin
=============================================================================
  ANALISIS UTAMA FINAL — POTENSI ENERGI & TITIK DOMINAN ARUS LAUT
  Delft3D-FLOW NEFIS trim.dat/trim.def
=============================================================================

FITUR UTAMA
-----------
A. Baca output Delft3D NEFIS secara benar.
B. Baca GRD/DEP/ENC/BND untuk peta dan masking.
C. Hapus spin-up.
D. Analisis statistik domain dan power density.
E. Audit spike kecepatan.
F. Dua ranking Top 5:
   1) Top 5 Titik Dominan Arus
      - prioritas frekuensi V >= v_medium
      - tie-break V_mean dan V_P95
   2) Top 5 Potensi Energi
      - prioritas mean power density
G. Peta, grafik, CSV, ringkasan TXT.
H. Validasi opsional terhadap QuickPlot velocity pada titik kontrol.
I. Validasi RMSE muka air S1 terhadap CSV observasi pada grid OBS/manual.
J. Analisis Spring/Neap berbasis S1 referensi.
K. GIF opsional dengan interval frame harian, bukan setiap timestep.

LIBRARY
-------
pip install numpy pandas matplotlib imageio utide
(utide hanya perlu untuk analisis Formzahl; imageio hanya untuk GIF)

=============================================================================
"""

from __future__ import annotations

import os
import sys
import csv
import re
import math
import ctypes
import warnings
import traceback
from pathlib import Path
from ctypes import c_int, c_char, byref, create_string_buffer
from datetime import datetime

import numpy as np
import pandas as pd

try:
    from utide import solve
except ImportError:
    solve = None
    
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, BoundaryNorm

try:
    import imageio.v2 as imageio
except Exception:
    imageio = None

warnings.filterwarnings("ignore", category=RuntimeWarning)


# =============================================================================
# ★ KONFIGURASI — SESUAIKAN BAGIAN INI ★
#
# CATATAN UNTUK GITHUB:
# - Path NEFIS DLL Delft3D di bawah sengaja dipertahankan.
# - Isi semua path r"" pada bagian FILE OUTPUT, GRID/BATIMETRI/BOUNDARY,
#   VALIDASI, dan OUTPUT sesuai lokasi file di komputer Anda.
# - Jangan mengunggah file data mentah Delft3D (.dat/.def/.grd/.dep/.enc/.bnd/.obs)
#   ke repository jika data tersebut bersifat pribadi/internal.

# =============================================================================

CONFIG = {
    # -------------------------------------------------------------------------
    # IDENTITAS STUDI
    # -------------------------------------------------------------------------
    "nama_selat": "Selat Capalulu",
    "kode_selat": "Capalulu",

    # -------------------------------------------------------------------------
    # NEFIS DLL DELFT3D
    # -------------------------------------------------------------------------
    "nefis_dll": r"C:\Program Files\Deltares\Delft3D 4.05.01\x64\share\bin\nefis.dll",
    "nefis_dir": r"C:\Program Files\Deltares\Delft3D 4.05.01\x64\share\bin",

    # -------------------------------------------------------------------------
    # FILE OUTPUT DELFT3D — PASANGAN trim.dat / trim.def
    # -------------------------------------------------------------------------
    "trim_dat": r"",  # <-- ISI path file trim-capalulu.dat
    "trim_def": r"",  # <-- ISI path file trim-capalulu.def

    # -------------------------------------------------------------------------
    # FILE GRID / BATIMETRI / BOUNDARY
    # -------------------------------------------------------------------------
    'grd_file'      : r'',  # <-- ISI path file .grd
    'dep_file'      : r'',  # <-- ISI path file .dep
    'enc_file'      : r'',  # <-- ISI path file .enc
    'bnd_file'      : r'',  # <-- ISI path file .bnd


    # -------------------------------------------------------------------------
    # WAKTU SIMULASI
    # -------------------------------------------------------------------------
    "tanggal_mulai": "2026-03-01 00:00:00",  # <-- sesuaikan dengan tanggal mulai simulasi
    "interval_menit": 60,
    "spinup_hari": 1,

    # Jika None, timestep dihitung otomatis dari file NEFIS.
    # Isi angka jika ingin debug cepat.
    "maks_timestep": None,

    # -------------------------------------------------------------------------
    # AMBANG KECEPATAN
    # -------------------------------------------------------------------------
    "v_cutin": 0.40,     # m/s
    "v_medium": 0.75,    # m/s
    "v_tinggi": 1.00,    # m/s

    # -------------------------------------------------------------------------
    # KRITERIA KEDALAMAN PENEMPATAN TURBIN
    # -------------------------------------------------------------------------
    "depth_min": 10.0,   # m
    "depth_max": 40.0,   # m

    # -------------------------------------------------------------------------
    # PARAMETER FISIK
    # -------------------------------------------------------------------------
    "rho": 1025.0,       # kg/m^3

    # -------------------------------------------------------------------------
    # KONTROL MASKING DAN SPIKE
    # -------------------------------------------------------------------------
    "missing_abs_threshold": 900.0,   # |-999.999| dll dianggap missing
    "spike_alert_mps": 5.0,           # audit, tidak dihapus otomatis
    "max_spike_rows_csv": 200000,     # batasi CSV spike agar tidak terlalu besar

    # -------------------------------------------------------------------------
    # RANKING TOP 5
    # -------------------------------------------------------------------------
    "jumlah_top": 5,
    "exclusion_radius_sel": 3,

    # Jika True, ranking hanya dari zona kedalaman layak + P95 >= cut-in.
    # False = ranking seluruh sel valid, tetapi status layak tetap dilaporkan.
    "ranking_dominan_arus_hanya_zona_layak": False,
    "ranking_potensi_energi_hanya_zona_layak": False,

    # -------------------------------------------------------------------------
    # VALIDASI VELOCITY vs QUICKPLOT — OPSIONAL
    # -------------------------------------------------------------------------
    "validasi_quickplot_aktif": True,
    "quickplot_velocity_csv": r"",  # <-- isi path CSV hasil QuickPlot jika validasi diaktifkan
    "quickplot_M": 176,
    "quickplot_N": 79,

    # -------------------------------------------------------------------------
    # VALIDASI RMSE ELEVASI MUKA AIR S1 — OPSIONAL
    # -------------------------------------------------------------------------
    "validasi_rmse_s1_aktif": False,

    # Mode titik S1 untuk RMSE:
    #   "manual" -> pakai rmse_s1_M / rmse_s1_N
    #   "obs_file" -> cari nama stasiun di file OBS
    "rmse_s1_mode": "obs_file",
    "obs_file": r"",  # <-- isi path file OBS jika validasi RMSE S1 diaktifkan
    "obs_station_name": "BITUNG",
    "rmse_s1_M": None,
    "rmse_s1_N": None,

    # CSV observasi S1. Header yang dikenali:
    #   waktu,eta_obs
    # atau
    #   date and time,water level (m)
    "obs_s1_csv": r"",
    "rmse_resample_freq": "10min",

    # -------------------------------------------------------------------------
    # ANALISIS SPRING/NEAP BERBASIS S1 — OPSIONAL
    # -------------------------------------------------------------------------
    "analisis_spring_neap_aktif": True,

    # Referensi S1:
    #   "p1_dominan_arus" -> P1 dari ranking dominan arus
    #   "p1_potensi_energi" -> P1 dari ranking energi
    #   "manual" -> gunakan spring_neap_M/N
    "spring_neap_ref": "p1_dominan_arus",
    "spring_neap_M": None,
    "spring_neap_N": None,
    "spring_neap_window_jam": 25,
    "spring_neap_percentile_spring": 80,
    "spring_neap_percentile_neap": 20,
    # -------------------------------------------------------------------------
    # ANALISIS HARMONIK PASUT / FORMZAHL — OPSIONAL
    # -------------------------------------------------------------------------
    "analisis_formzahl_aktif": True,

    # Data minimal yang disarankan agar pemisahan komponen utama lebih stabil.
    # Data Anda 29 hari setelah spin-up sudah mencukupi.
    "formzahl_min_durasi_hari": 15,
    # -------------------------------------------------------------------------
    # VISUALISASI & OUTPUT
    # -------------------------------------------------------------------------
    "output_dir": r"",  # <-- isi folder output hasil analisis di komputer Anda
    "colormap": "jet",
    "skip_quiver": 4,
    "skala_quiver": 5,
    "dpi": 200,

    # GIF — dinonaktifkan default agar tidak berat.
    # Jika aktif, frame diambil per hari, bukan setiap timestep.
    "buat_gif": False,
    "gif_frame_interval_hari": 1.0,
    "gif_fps": 3,

    # Sampel maksimum histogram agar hemat memori
    "hist_max_sample": 2000000,
}

C = CONFIG

# =============================================================================
# KONSTANTA GLOBAL
# =============================================================================

NAME_BUF = 17
WARNA_TOP = ["#e63946", "#f4a261", "#2a9d8f", "#457b9d", "#8338ec",
             "#9b5de5", "#00bbf9", "#f15bb5", "#fee440", "#00f5d4"]

if not str(C["output_dir"]).strip():
    raise SystemExit(
        "CONFIG['output_dir'] masih kosong. Isi folder output di bagian KONFIGURASI."
    )
OUTPUT_DIR = Path(C["output_dir"])
SUBDIRS = {
    "peta": OUTPUT_DIR / "peta",
    "grafik": OUTPUT_DIR / "grafik",
    "csv": OUTPUT_DIR / "csv",
    "gif": OUTPUT_DIR / "animasi",
}
for d in [OUTPUT_DIR, *SUBDIRS.values()]:
    d.mkdir(parents=True, exist_ok=True)

LOG_LINES: list[str] = []



# =============================================================================
# WATERMARK
# =============================================================================

WATERMARK_TEXT = "© Aimar Rendra P.A. | Teknik Kelautan — Universitas Hasanuddin"


def savefig_wm(fp, **kwargs) -> None:
    """Simpan figure aktif dengan watermark di pojok kanan bawah."""
    fig = plt.gcf()
    fig.text(
        0.99, 0.005, WATERMARK_TEXT,
        ha="right", va="bottom", fontsize=7, color="0.35", alpha=0.75,
        zorder=100,
    )
    plt.savefig(fp, **kwargs)

# =============================================================================
# UTILITAS UMUM
# =============================================================================

def log(text: str = "") -> None:
    print(text)
    LOG_LINES.append(str(text))


def section(title: str) -> None:
    log("")
    log("─" * 78)
    log(f"  {title}")
    log("─" * 78)


def safe_name(text: str) -> str:
    return str(text).strip().replace(" ", "_").replace("/", "_").replace("\\", "_")


def finite_count(arr: np.ndarray) -> int:
    return int(np.sum(np.isfinite(arr)))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if fieldnames:
            writer.writeheader()
        for row in rows:
            writer.writerow(row)


def to_float_or_nan(v):
    try:
        return float(v)
    except Exception:
        return np.nan


def nan_stat(arr: np.ndarray, fn, default=np.nan):
    if np.any(np.isfinite(arr)):
        return float(fn(arr))
    return default


def flatten_valid_sample(arr: np.ndarray, max_sample: int, seed: int = 42) -> np.ndarray:
    vals = arr[np.isfinite(arr)]
    if vals.size <= max_sample:
        return vals
    rng = np.random.default_rng(seed)
    idx = rng.choice(vals.size, size=max_sample, replace=False)
    return vals[idx]


def create_model_time(n: int) -> pd.DatetimeIndex:
    start = pd.Timestamp(C["tanggal_mulai"]) + pd.Timedelta(days=C["spinup_hari"])
    return pd.date_range(
        start=start,
        periods=n,
        freq=pd.to_timedelta(C["interval_menit"], unit="m"),
    )


# =============================================================================
# NEFIS DLL
# =============================================================================

def load_nefis():
    os.environ["PATH"] = C["nefis_dir"] + os.pathsep + os.environ.get("PATH", "")
    if hasattr(os, "add_dll_directory"):
        try:
            os.add_dll_directory(C["nefis_dir"])
        except Exception:
            pass

    try:
        nef = ctypes.cdll.LoadLibrary(C["nefis_dll"])
    except Exception as exc:
        raise RuntimeError(f"Gagal load NEFIS DLL: {exc}") from exc

    nef.Crenef.restype = c_int
    nef.Clsnef.restype = c_int
    nef.Getelt.restype = c_int
    return nef


def nefis_open(nef, dat_path: str, def_path: str):
    fd = c_int(0)
    ret = nef.Crenef(
        byref(fd),
        dat_path.encode(),
        def_path.encode(),
        c_char(b"N"),
        c_char(b"r"),
    )
    if ret != 0:
        raise RuntimeError(
            f"Gagal buka NEFIS kode={ret}. "
            f"Pastikan DAT/DEF berpasangan.\nDAT={dat_path}\nDEF={def_path}"
        )
    return fd


def nefis_close(nef, fd) -> None:
    try:
        nef.Clsnef(byref(fd))
    except Exception:
        pass


def read_int_scalar(nef, fd, group: str, element: str) -> int:
    g = create_string_buffer(group.encode(), NAME_BUF)
    e = create_string_buffer(element.encode(), NAME_BUF)
    ui = (c_int * 3)(1, 1, 1)
    uo = (c_int * 1)(1)
    bl = c_int(4)
    buf = (c_int * 1)()
    ret = nef.Getelt(byref(fd), g, e, ui, uo, byref(bl), buf)
    if ret != 0:
        raise RuntimeError(f"Getelt gagal {group}/{element}, kode={ret}")
    return int(buf[0])


def read_series_element_time_only(
    nef,
    fd,
    element: str,
    t0: int,
    count_float: int,
) -> np.ndarray:
    """
    Membaca satu timestep penuh dari group map-series.
    TIME-only sesuai metadata:
      UINDEX = (t+1, t+1, 1)
    """
    g = create_string_buffer(b"map-series", NAME_BUF)
    e = create_string_buffer(element.encode(), NAME_BUF)
    ui = (c_int * 3)(t0 + 1, t0 + 1, 1)
    uo = (c_int * 1)(1)
    bl = c_int(count_float * 4)
    buf = (ctypes.c_float * count_float)()
    ret = nef.Getelt(byref(fd), g, e, ui, uo, byref(bl), buf)
    if ret != 0:
        raise EOFError(f"Getelt map-series/{element} gagal pada timestep {t0+1}, kode={ret}")
    return np.array(buf[:count_float], dtype=np.float32)


def read_const_element_time_only(
    nef,
    fd,
    element: str,
    count_float: int,
) -> np.ndarray:
    """
    Membaca element map-const berukuran array penuh. Tidak selalu dipakai;
    fallback jika DEP/GRD tidak tersedia.
    """
    g = create_string_buffer(b"map-const", NAME_BUF)
    e = create_string_buffer(element.encode(), NAME_BUF)
    ui = (c_int * 3)(1, 1, 1)
    uo = (c_int * 1)(1)
    bl = c_int(count_float * 4)
    buf = (ctypes.c_float * count_float)()
    ret = nef.Getelt(byref(fd), g, e, ui, uo, byref(bl), buf)
    if ret != 0:
        raise RuntimeError(f"Getelt map-const/{element} gagal, kode={ret}")
    return np.array(buf[:count_float], dtype=np.float32)


def count_timesteps(nef, fd, count_float_uv: int) -> int:
    section("BAGIAN 3: Menghitung timestep map-series")
    n = 0
    max_ts = C["maks_timestep"]
    while True:
        if max_ts is not None and n >= int(max_ts):
            log(f"  ○ Penghitungan dihentikan oleh CONFIG maks_timestep={max_ts}")
            break
        try:
            _ = read_series_element_time_only(nef, fd, "U1", n, count_float_uv)
            n += 1
            if n % 500 == 0:
                log(f"    {n} timestep terdeteksi...")
        except EOFError:
            break
    log(f"  ✓ Timestep ditemukan: {n}")
    return n


# =============================================================================
# PARSER GRD / DEP / ENC / BND / OBS
# =============================================================================

def parse_grd_file(filepath: str | Path | None):
    """
    Parser GRD Delft3D/RGFGRID yang mendukung format seperti:
      Coordinate System = Spherical
      Missing Value     = -9.99999E+02
           MMAX NMAX
       0 0 0
       ETA= 1  X_1 X_2 ...
              X_...
       ...
       ETA= NMAX ...
       ETA= 1  Y_1 Y_2 ...
       ...
       ETA= NMAX ...

    Pada format ini TIDAK terdapat label 'X=' dan 'Y='.
    Blok ETA pertama sebanyak NMAX baris adalah koordinat X,
    blok ETA kedua sebanyak NMAX baris adalah koordinat Y.
    """
    if not filepath or not Path(filepath).exists():
        log("  ○ GRD tidak tersedia.")
        return None, None

    path = Path(filepath)
    log(f"  → Membaca GRD: {path.name}")
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()

    # ---------------------------------------------------------------------
    # 1. Deteksi dimensi MMAX, NMAX
    # ---------------------------------------------------------------------
    mmax = nmax = None
    dim_line_idx = None
    for i, line in enumerate(lines):
        parts = line.strip().split()
        if len(parts) == 2:
            try:
                a, b = int(parts[0]), int(parts[1])
                if a > 1 and b > 1:
                    mmax, nmax = a, b
                    dim_line_idx = i
                    break
            except Exception:
                continue

    if mmax is None or nmax is None:
        log("  ⚠ Dimensi GRD tidak terdeteksi.")
        return None, None

    # ---------------------------------------------------------------------
    # 2. Ambil semua blok ETA
    # ---------------------------------------------------------------------
    eta_blocks: list[tuple[int, list[float]]] = []
    i = dim_line_idx + 1
    while i < len(lines):
        line = lines[i]
        if "ETA=" not in line.upper():
            i += 1
            continue

        match = re.search(r"ETA=\s*(\d+)", line, flags=re.IGNORECASE)
        if not match:
            i += 1
            continue

        eta = int(match.group(1))
        values: list[float] = []

        # Nilai setelah ETA= ...
        tail = line[match.end():]
        for tok in tail.replace(",", " ").split():
            try:
                values.append(float(tok))
            except Exception:
                pass

        # Baris lanjutan hingga ETA berikutnya
        j = i + 1
        while j < len(lines) and "ETA=" not in lines[j].upper():
            for tok in lines[j].replace(",", " ").split():
                try:
                    values.append(float(tok))
                except Exception:
                    pass
            j += 1

        eta_blocks.append((eta, values))
        i = j

    if len(eta_blocks) < 2 * nmax:
        log(
            f"  ⚠ Blok ETA GRD tidak cukup: ditemukan {len(eta_blocks)}, "
            f"dibutuhkan minimal {2*nmax} untuk X dan Y."
        )
        return None, None

    # ---------------------------------------------------------------------
    # 3. Blok pertama = X, blok kedua = Y
    # ---------------------------------------------------------------------
    X = np.full((nmax, mmax), np.nan, dtype=float)
    Y = np.full((nmax, mmax), np.nan, dtype=float)

    for eta, vals in eta_blocks[:nmax]:
        r = eta - 1
        if 0 <= r < nmax:
            nfill = min(mmax, len(vals))
            X[r, :nfill] = vals[:nfill]

    for eta, vals in eta_blocks[nmax:2*nmax]:
        r = eta - 1
        if 0 <= r < nmax:
            nfill = min(mmax, len(vals))
            Y[r, :nfill] = vals[:nfill]

    # Missing value Delphi/Delft3D
    X = np.where(np.abs(X) >= C["missing_abs_threshold"], np.nan, X)
    Y = np.where(np.abs(Y) >= C["missing_abs_threshold"], np.nan, Y)

    valid = np.isfinite(X) & np.isfinite(Y)
    if not np.any(valid):
        log("  ⚠ GRD terbaca tetapi tidak ada koordinat valid.")
        return None, None

    log(
        f"    ✓ GRD dimensi MMAX={mmax}, NMAX={nmax} | "
        f"koordinat valid={int(np.sum(valid))} | "
        f"X={np.nanmin(X):.6f}→{np.nanmax(X):.6f} | "
        f"Y={np.nanmin(Y):.6f}→{np.nanmax(Y):.6f}"
    )
    return X, Y


def parse_dep_file(filepath: str | Path | None, target_shape: tuple[int, int]):
    if not filepath or not Path(filepath).exists():
        log("  ○ DEP tidak tersedia.")
        return None

    path = Path(filepath)
    log(f"  → Membaca DEP: {path.name}")

    nums: list[float] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        for tok in line.replace(",", " ").split():
            try:
                nums.append(float(tok))
            except Exception:
                pass

    n_target, m_target = target_shape
    candidate_shapes = [
        (n_target, m_target),
        (n_target + 1, m_target + 1),
        (n_target + 1, m_target),
        (n_target, m_target + 1),
    ]
    arr = None
    shape_found = None
    for shp in candidate_shapes:
        if len(nums) == shp[0] * shp[1]:
            arr = np.array(nums, dtype=float).reshape(shp)
            shape_found = shp
            break

    if arr is None:
        log(
            f"  ⚠ Jumlah nilai DEP={len(nums)} tidak cocok dengan target {target_shape} "
            f"atau varian +1 baris/kolom."
        )
        return None

    if arr.shape != target_shape:
        arr = arr[:n_target, :m_target]
        log(f"    ✓ DEP shape {shape_found} diselaraskan menjadi {arr.shape}")
    else:
        log(f"    ✓ DEP shape {arr.shape}")

    arr = np.abs(arr)
    arr = np.where(np.abs(arr) >= C["missing_abs_threshold"], np.nan, arr)
    valid = np.isfinite(arr)
    if np.any(valid):
        log(
            f"    ✓ DEP valid={int(np.sum(valid))} | "
            f"min={np.nanmin(arr):.3f} m | max={np.nanmax(arr):.3f} m"
        )
    else:
        log("    ⚠ DEP tidak memiliki nilai valid.")
    return arr


def parse_enc_file(filepath: str | Path | None, X: np.ndarray | None, Y: np.ndarray | None):
    if not filepath or not Path(filepath).exists():
        log("  ○ ENC tidak tersedia.")
        return None

    path = Path(filepath)
    coords = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith("*"):
            continue
        parts = s.split()
        if len(parts) < 2:
            continue
        try:
            m = int(parts[0]) - 1
            n = int(parts[1]) - 1
            if X is not None and Y is not None and 0 <= n < X.shape[0] and 0 <= m < X.shape[1]:
                x, y = X[n, m], Y[n, m]
                if np.isfinite(x) and np.isfinite(y):
                    coords.append((float(x), float(y)))
        except Exception:
            continue
    log(f"  ✓ ENC: {len(coords)} titik closed boundary")
    return coords or None


def parse_bnd_file(filepath: str | Path | None, X: np.ndarray | None, Y: np.ndarray | None):
    if not filepath or not Path(filepath).exists():
        log("  ○ BND tidak tersedia.")
        return None

    path = Path(filepath)
    segs = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith("*"):
            continue
        ints = []
        for tok in s.split():
            try:
                ints.append(int(tok))
            except Exception:
                pass
        if len(ints) >= 4 and X is not None and Y is not None:
            m1, n1, m2, n2 = ints[0]-1, ints[1]-1, ints[2]-1, ints[3]-1
            if (0 <= n1 < X.shape[0] and 0 <= m1 < X.shape[1] and
                0 <= n2 < X.shape[0] and 0 <= m2 < X.shape[1]):
                x1, y1 = X[n1, m1], Y[n1, m1]
                x2, y2 = X[n2, m2], Y[n2, m2]
                if all(np.isfinite(v) for v in [x1, y1, x2, y2]):
                    segs.append(((float(x1), float(y1)), (float(x2), float(y2))))
    log(f"  ✓ BND: {len(segs)} segmen open boundary")
    return segs or None


def parse_obs_file(filepath: str | Path, station_name: str):
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"File OBS tidak ditemukan: {path}")
    station_norm = station_name.strip().lower()
    candidates = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith("*"):
            continue
        parts = s.split()
        if len(parts) < 3:
            continue
        try:
            m = int(parts[-2])
            n = int(parts[-1])
            name = " ".join(parts[:-2]).strip()
            candidates.append((name, m, n))
            if name.strip().lower() == station_norm:
                return m, n, candidates
        except Exception:
            continue
    names = ", ".join([c[0] for c in candidates[:20]])
    raise ValueError(
        f"Stasiun '{station_name}' tidak ditemukan di OBS. "
        f"Contoh stasiun tersedia: {names}"
    )


# =============================================================================
# VELOCITY RECONSTRUCTION — METODE TERVALIDASI
# =============================================================================

def clean_values(arr: np.ndarray) -> np.ndarray:
    arr = arr.astype(np.float32, copy=False)
    arr = np.where(np.abs(arr) >= C["missing_abs_threshold"], np.nan, arr)
    arr = np.where(np.abs(arr) > 1e10, np.nan, arr)
    return arr.astype(np.float32, copy=False)


def reshape_uv_f_order(flat: np.ndarray, nmax: int, mmax: int, kmax: int) -> np.ndarray:
    return flat.reshape((nmax, mmax, kmax), order="F")[:, :, 0]


def reshape_s1_f_order(flat: np.ndarray, nmax: int, mmax: int) -> np.ndarray:
    return flat.reshape((nmax, mmax), order="F")


def reconstruct_backward_velocity(U: np.ndarray, V: np.ndarray):
    """
    Metode tervalidasi terhadap QuickPlot:
      Uc[n,m] = 0.5 * (U[n,m-1] + U[n,m])
      Vc[n,m] = 0.5 * (V[n-1,m] + V[n,m])
    """
    Uc = np.full_like(U, np.nan, dtype=np.float32)
    Vc = np.full_like(V, np.nan, dtype=np.float32)

    Uc[:, 1:] = 0.5 * (U[:, :-1] + U[:, 1:])
    Vc[1:, :] = 0.5 * (V[:-1, :] + V[1:, :])

    MAG = np.sqrt(Uc**2 + Vc**2).astype(np.float32)
    return Uc, Vc, MAG


# =============================================================================
# PLOTTING
# =============================================================================

def coord_mode(X: np.ndarray | None, Y: np.ndarray | None) -> str:
    if X is None or Y is None:
        return "index"
    vals = X[np.isfinite(X)]
    if vals.size and np.nanmax(np.abs(vals)) < 360:
        return "geo"
    return "index"


def axes_labels(X, Y):
    if coord_mode(X, Y) == "geo":
        return "Longitude (°)", "Latitude (°)"
    return "Indeks Grid M", "Indeks Grid N"


def plot_field(ax, X, Y, Z, cmap, vmin=None, vmax=None, label=None, s=8):
    """
    Aman terhadap koordinat NaN:
    - Jika X/Y seluruhnya finite -> pcolormesh
    - Jika ada NaN -> scatter hanya pada titik valid
    """
    if X is None or Y is None:
        n, m = Z.shape
        XX, YY = np.meshgrid(np.arange(1, m+1), np.arange(1, n+1))
        X, Y = XX, YY

    if np.all(np.isfinite(X)) and np.all(np.isfinite(Y)):
        artist = ax.pcolormesh(X, Y, Z, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    else:
        valid = np.isfinite(X) & np.isfinite(Y) & np.isfinite(Z)
        artist = ax.scatter(X[valid], Y[valid], c=Z[valid], s=s, cmap=cmap, vmin=vmin, vmax=vmax, marker="s")
    return artist


def add_boundaries(ax, enc_coords, bnd_segs):
    if enc_coords and len(enc_coords) >= 2:
        xs = [p[0] for p in enc_coords] + [enc_coords[0][0]]
        ys = [p[1] for p in enc_coords] + [enc_coords[0][1]]
        ax.plot(xs, ys, color="black", lw=1.8, alpha=0.9, zorder=4, label="Closed Boundary")
    if bnd_segs:
        first = True
        for (x1, y1), (x2, y2) in bnd_segs:
            ax.plot([x1, x2], [y1, y2], color="red", lw=2.0, ls="--",
                    label="Open Boundary" if first else None, zorder=5)
            first = False


def plot_top_points(ax, points: list[dict], label_field="v_mean"):
    for p in points:
        color = WARNA_TOP[(int(p["no"]) - 1) % len(WARNA_TOP)]
        ax.plot(p["x"], p["y"], marker="*", markersize=16, color=color, zorder=8)
        text = f"P{p['no']}\n{p[label_field]:.3f}"
        ax.annotate(
            text,
            (p["x"], p["y"]),
            xytext=(6, 5),
            textcoords="offset points",
            fontsize=8,
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.4", alpha=0.8),
            zorder=9,
        )


# =============================================================================
# STATISTIK DAN RANKING
# =============================================================================

def safe_percentile_axis(arr: np.ndarray, q: float) -> np.ndarray:
    return np.nanpercentile(arr, q, axis=0)


def build_grid_statistics(
    MAG: np.ndarray,
    POWER_mean: np.ndarray,
    POWER_p95: np.ndarray,
    DEPTH: np.ndarray,
    X: np.ndarray,
    Y: np.ndarray,
    valid_mask: np.ndarray,
):
    v_min = np.nanmin(MAG, axis=0)
    v_max = np.nanmax(MAG, axis=0)
    v_mean = np.nanmean(MAG, axis=0)
    v_median = np.nanpercentile(MAG, 50, axis=0)
    v_p90 = np.nanpercentile(MAG, 90, axis=0)
    v_p95 = np.nanpercentile(MAG, 95, axis=0)
    v_std = np.nanstd(MAG, axis=0)

    valid_ts_count = np.sum(np.isfinite(MAG), axis=0)
    freq_cutin = np.sum(MAG >= C["v_cutin"], axis=0) / np.maximum(valid_ts_count, 1) * 100.0
    freq_medium = np.sum(MAG >= C["v_medium"], axis=0) / np.maximum(valid_ts_count, 1) * 100.0
    freq_tinggi = np.sum(MAG >= C["v_tinggi"], axis=0) / np.maximum(valid_ts_count, 1) * 100.0

    stats = {
        "v_min": v_min,
        "v_max": v_max,
        "v_mean": v_mean,
        "v_median": v_median,
        "v_p90": v_p90,
        "v_p95": v_p95,
        "v_std": v_std,
        "freq_ge_cutin_pct": freq_cutin,
        "freq_ge_medium_pct": freq_medium,
        "freq_ge_tinggi_pct": freq_tinggi,
        "power_mean": POWER_mean,
        "power_p95": POWER_p95,
        "valid_ts_count": valid_ts_count,
    }

    rows = []
    ns, ms = np.where(valid_mask)
    for n0, m0 in zip(ns, ms):
        rows.append({
            "M": int(m0 + 1),
            "N": int(n0 + 1),
            "python_n0": int(n0),
            "python_m0": int(m0),
            "longitude_x": to_float_or_nan(X[n0, m0]) if X is not None else np.nan,
            "latitude_y": to_float_or_nan(Y[n0, m0]) if Y is not None else np.nan,
            "depth_m": to_float_or_nan(DEPTH[n0, m0]) if DEPTH is not None else np.nan,
            **{k: to_float_or_nan(v[n0, m0]) for k, v in stats.items() if k != "valid_ts_count"},
            "valid_ts_count": int(valid_ts_count[n0, m0]),
        })
    return stats, rows


def define_layak_mask(valid_mask: np.ndarray, v_p95: np.ndarray, DEPTH: np.ndarray | None):
    velocity_ok = v_p95 >= C["v_cutin"]
    if DEPTH is None or not np.any(np.isfinite(DEPTH)):
        depth_ok = np.ones_like(valid_mask, dtype=bool)
        depth_note = "DEP tidak tersedia; syarat kedalaman tidak diterapkan."
    else:
        depth_ok = np.isfinite(DEPTH) & (DEPTH >= C["depth_min"]) & (DEPTH <= C["depth_max"])
        depth_note = "Syarat kedalaman diterapkan."
    layak = valid_mask & velocity_ok & depth_ok
    optimal = layak & (v_p95 >= C["v_medium"])
    return layak, optimal, depth_note


def exclusion_mask(mask: np.ndarray, n0: int, m0: int, radius: int) -> None:
    nmin = max(0, n0 - radius)
    nmax = min(mask.shape[0], n0 + radius + 1)
    mmin = max(0, m0 - radius)
    mmax = min(mask.shape[1], m0 + radius + 1)
    mask[nmin:nmax, mmin:mmax] = False


def point_from_cell(
    no: int,
    n0: int,
    m0: int,
    MAG: np.ndarray,
    DEPTH: np.ndarray | None,
    X: np.ndarray | None,
    Y: np.ndarray | None,
    stats: dict[str, np.ndarray],
    layak_mask: np.ndarray,
):
    ts = MAG[:, n0, m0].astype(float)
    ts_valid = ts[np.isfinite(ts)]
    pwr_ts = 0.5 * C["rho"] * ts**3

    return {
        "no": int(no),
        "M": int(m0 + 1),
        "N": int(n0 + 1),
        "python_n0": int(n0),
        "python_m0": int(m0),
        "x": to_float_or_nan(X[n0, m0]) if X is not None else float(m0 + 1),
        "y": to_float_or_nan(Y[n0, m0]) if Y is not None else float(n0 + 1),
        "depth_m": to_float_or_nan(DEPTH[n0, m0]) if DEPTH is not None else np.nan,
        "v_min": nan_stat(ts_valid, np.nanmin),
        "v_max": nan_stat(ts_valid, np.nanmax),
        "v_mean": nan_stat(ts_valid, np.nanmean),
        "v_median": nan_stat(ts_valid, lambda a: np.nanpercentile(a, 50)),
        "v_p90": nan_stat(ts_valid, lambda a: np.nanpercentile(a, 90)),
        "v_p95": nan_stat(ts_valid, lambda a: np.nanpercentile(a, 95)),
        "v_std": nan_stat(ts_valid, np.nanstd),
        "freq_ge_cutin_pct": to_float_or_nan(stats["freq_ge_cutin_pct"][n0, m0]),
        "freq_ge_medium_pct": to_float_or_nan(stats["freq_ge_medium_pct"][n0, m0]),
        "freq_ge_tinggi_pct": to_float_or_nan(stats["freq_ge_tinggi_pct"][n0, m0]),
        "power_mean_wm2": nan_stat(pwr_ts, np.nanmean),
        "power_p95_wm2": nan_stat(pwr_ts, lambda a: np.nanpercentile(a, 95)),
        "power_max_wm2": nan_stat(pwr_ts, np.nanmax),
        "layak_zona": bool(layak_mask[n0, m0]),
        "_ts_mag": ts,
        "_ts_pwr": pwr_ts,
    }


def rank_top_dominant_current(
    MAG: np.ndarray,
    valid_mask: np.ndarray,
    layak_mask: np.ndarray,
    stats: dict[str, np.ndarray],
    DEPTH: np.ndarray | None,
    X: np.ndarray | None,
    Y: np.ndarray | None,
):
    candidate = layak_mask.copy() if C["ranking_dominan_arus_hanya_zona_layak"] else valid_mask.copy()
    rows = np.argwhere(candidate)
    df = pd.DataFrame({
        "n0": rows[:, 0] if len(rows) else [],
        "m0": rows[:, 1] if len(rows) else [],
    })
    if df.empty:
        return []

    df["freq_medium"] = [stats["freq_ge_medium_pct"][int(n), int(m)] for n, m in zip(df["n0"], df["m0"])]
    df["v_mean"] = [stats["v_mean"][int(n), int(m)] for n, m in zip(df["n0"], df["m0"])]
    df["v_p95"] = [stats["v_p95"][int(n), int(m)] for n, m in zip(df["n0"], df["m0"])]
    df = df.sort_values(["freq_medium", "v_mean", "v_p95"], ascending=[False, False, False])

    active = candidate.copy()
    points = []
    for _, row in df.iterrows():
        n0, m0 = int(row["n0"]), int(row["m0"])
        if not active[n0, m0]:
            continue
        points.append(point_from_cell(len(points)+1, n0, m0, MAG, DEPTH, X, Y, stats, layak_mask))
        exclusion_mask(active, n0, m0, int(C["exclusion_radius_sel"]))
        if len(points) >= int(C["jumlah_top"]):
            break
    return points


def rank_top_energy(
    MAG: np.ndarray,
    valid_mask: np.ndarray,
    layak_mask: np.ndarray,
    stats: dict[str, np.ndarray],
    DEPTH: np.ndarray | None,
    X: np.ndarray | None,
    Y: np.ndarray | None,
):
    candidate = layak_mask.copy() if C["ranking_potensi_energi_hanya_zona_layak"] else valid_mask.copy()
    rows = np.argwhere(candidate)
    df = pd.DataFrame({
        "n0": rows[:, 0] if len(rows) else [],
        "m0": rows[:, 1] if len(rows) else [],
    })
    if df.empty:
        return []

    df["power_mean"] = [stats["power_mean"][int(n), int(m)] for n, m in zip(df["n0"], df["m0"])]
    df["v_mean"] = [stats["v_mean"][int(n), int(m)] for n, m in zip(df["n0"], df["m0"])]
    df = df.sort_values(["power_mean", "v_mean"], ascending=[False, False])

    active = candidate.copy()
    points = []
    for _, row in df.iterrows():
        n0, m0 = int(row["n0"]), int(row["m0"])
        if not active[n0, m0]:
            continue
        points.append(point_from_cell(len(points)+1, n0, m0, MAG, DEPTH, X, Y, stats, layak_mask))
        exclusion_mask(active, n0, m0, int(C["exclusion_radius_sel"]))
        if len(points) >= int(C["jumlah_top"]):
            break
    return points


def public_point_rows(points: list[dict]) -> list[dict]:
    rows = []
    for p in points:
        rows.append({k: v for k, v in p.items() if not k.startswith("_")})
    return rows


# =============================================================================
# VALIDASI QUICKPLOT VELOCITY
# =============================================================================

def load_velocity_quickplot_csv(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    df = pd.read_csv(path)

    time_col = None
    vel_col = None
    for col in df.columns:
        cl = col.strip().lower()
        if cl in {"date and time", "waktu", "datetime", "time", "date_time"}:
            time_col = col
        if "magnitude of depth averaged velocity" in cl:
            vel_col = col
        if cl in {"velocity", "vel", "mag", "magnitude"}:
            vel_col = vel_col or col

    if time_col is None or vel_col is None:
        raise ValueError(
            "Kolom QuickPlot velocity tidak dikenali. "
            "Butuh kolom waktu dan magnitude velocity."
        )

    out = pd.DataFrame({
        "waktu": pd.to_datetime(df[time_col], errors="coerce"),
        "velocity_quickplot": pd.to_numeric(df[vel_col], errors="coerce"),
    }).dropna(subset=["waktu"]).sort_values("waktu").reset_index(drop=True)
    return out


def compare_model_quickplot_velocity(
    MAG: np.ndarray,
    WAKTU: pd.DatetimeIndex,
    M: int,
    N: int,
    csv_path: str | Path,
):
    n0, m0 = N - 1, M - 1
    if n0 < 0 or m0 < 0 or n0 >= MAG.shape[1] or m0 >= MAG.shape[2]:
        raise ValueError(f"Titik validasi QuickPlot M={M}, N={N} di luar grid analisis.")

    model_df = pd.DataFrame({
        "waktu": pd.to_datetime(WAKTU),
        "velocity_model": MAG[:, n0, m0].astype(float),
    })
    qp = load_velocity_quickplot_csv(csv_path)

    merged = model_df.merge(qp, on="waktu", how="inner").dropna()
    if merged.empty:
        raise ValueError("Tidak ada overlap waktu antara model dan CSV QuickPlot.")

    resid = merged["velocity_model"] - merged["velocity_quickplot"]
    rmse = float(np.sqrt(np.mean(resid**2)))
    mae = float(np.mean(np.abs(resid)))
    bias = float(np.mean(resid))
    r = float(np.corrcoef(merged["velocity_model"], merged["velocity_quickplot"])[0, 1]) if len(merged) > 1 else np.nan

    summary = {
        "M": M,
        "N": N,
        "N_data": int(len(merged)),
        "RMSE_mps": rmse,
        "MAE_mps": mae,
        "Bias_model_minus_quickplot_mps": bias,
        "r": r,
        "mean_model_mps": float(merged["velocity_model"].mean()),
        "mean_quickplot_mps": float(merged["velocity_quickplot"].mean()),
        "max_model_mps": float(merged["velocity_model"].max()),
        "max_quickplot_mps": float(merged["velocity_quickplot"].max()),
    }
    return summary, merged


# =============================================================================
# S1: EKSTRAKSI, RMSE OBSERVASI, SPRING/NEAP
# =============================================================================

def get_s1_request_points(
    top_current: list[dict],
    top_energy: list[dict],
):
    requests: dict[str, tuple[int, int]] = {}

    if C["validasi_rmse_s1_aktif"]:
        if C["rmse_s1_mode"] == "manual":
            if C["rmse_s1_M"] is None or C["rmse_s1_N"] is None:
                raise ValueError("rmse_s1_mode='manual' tetapi M/N belum diisi.")
            requests["rmse_obs"] = (int(C["rmse_s1_M"]), int(C["rmse_s1_N"]))
        elif C["rmse_s1_mode"] == "obs_file":
            m, n, _ = parse_obs_file(C["obs_file"], C["obs_station_name"])
            requests["rmse_obs"] = (m, n)
        else:
            raise ValueError("rmse_s1_mode harus 'manual' atau 'obs_file'.")

    if C["analisis_spring_neap_aktif"]:
        ref = C["spring_neap_ref"]
        if ref == "p1_dominan_arus":
            if top_current:
                requests["spring_neap"] = (int(top_current[0]["M"]), int(top_current[0]["N"]))
        elif ref == "p1_potensi_energi":
            if top_energy:
                requests["spring_neap"] = (int(top_energy[0]["M"]), int(top_energy[0]["N"]))
        elif ref == "manual":
            if C["spring_neap_M"] is None or C["spring_neap_N"] is None:
                raise ValueError("spring_neap_ref='manual' tetapi M/N belum diisi.")
            requests["spring_neap"] = (int(C["spring_neap_M"]), int(C["spring_neap_N"]))
        else:
            raise ValueError("spring_neap_ref tidak dikenali.")

    return requests


def extract_s1_timeseries(
    nef,
    fd,
    n_time: int,
    n_spinup: int,
    nmax_full: int,
    mmax_full: int,
    requests: dict[str, tuple[int, int]],
):
    if not requests:
        return {}

    section("BAGIAN 12: Ekstraksi S1 untuk RMSE / Spring-Neap")
    count_s1 = nmax_full * mmax_full
    series = {k: [] for k in requests}

    for t in range(n_time):
        try:
            flat = read_series_element_time_only(nef, fd, "S1", t, count_s1)
        except EOFError as exc:
            log(f"  ⚠ S1 berhenti pada timestep {t+1}: {exc}")
            break
        if t < n_spinup:
            continue
        s1 = clean_values(reshape_s1_f_order(flat, nmax_full, mmax_full))
        for key, (m, n) in requests.items():
            n0, m0 = n - 1, m - 1
            if 0 <= n0 < s1.shape[0] and 0 <= m0 < s1.shape[1]:
                series[key].append(float(s1[n0, m0]))
            else:
                series[key].append(np.nan)

        if (t + 1) % 1000 == 0 or (t + 1) == n_time:
            log(f"    S1 {t+1}/{n_time} timestep ✓")

    out = {k: np.asarray(v, dtype=float) for k, v in series.items()}
    for key, arr in out.items():
        log(f"  ✓ S1 '{key}' diekstrak: {len(arr)} timestep, valid={int(np.sum(np.isfinite(arr)))}")
    return out


def load_obs_s1_csv(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    df = pd.read_csv(path)

    time_col = None
    eta_col = None
    for col in df.columns:
        cl = col.strip().lower()
        if cl in {"waktu", "date and time", "datetime", "date_time", "time"}:
            time_col = col
        if cl in {"eta_obs", "water level (m)", "water_level", "water level", "eta"}:
            eta_col = col
    if time_col is None or eta_col is None:
        raise ValueError(
            "CSV observasi S1 harus memiliki kolom waktu dan eta_obs/water level (m)."
        )

    out = pd.DataFrame({
        "waktu": pd.to_datetime(df[time_col], errors="coerce"),
        "eta_obs": pd.to_numeric(df[eta_col], errors="coerce"),
    }).dropna(subset=["waktu"]).sort_values("waktu").reset_index(drop=True)
    return out


def rmse_s1_validation(s1_model: np.ndarray, WAKTU: pd.DatetimeIndex, obs_csv: str | Path):
    model = pd.DataFrame({
        "waktu": pd.to_datetime(WAKTU[:len(s1_model)]),
        "eta_model": s1_model,
    })
    obs = load_obs_s1_csv(obs_csv)

    freq = C["rmse_resample_freq"]
    model_r = model.set_index("waktu").resample(freq).mean()
    obs_r = obs.set_index("waktu")["eta_obs"].resample(freq).mean()
    comb = model_r.join(obs_r, how="inner").dropna()
    if comb.empty:
        raise ValueError("Tidak ada overlap waktu S1 model–observasi.")

    resid = comb["eta_model"] - comb["eta_obs"]
    rmse = float(np.sqrt(np.mean(resid**2)))
    mae = float(np.mean(np.abs(resid)))
    bias = float(np.mean(resid))
    r = float(np.corrcoef(comb["eta_model"], comb["eta_obs"])[0, 1]) if len(comb) > 1 else np.nan
    denom = float(np.sum((comb["eta_obs"] - comb["eta_obs"].mean()) ** 2))
    nse = 1.0 - float(np.sum(resid**2)) / denom if denom > 0 else np.nan

    summary = {
        "N_data": int(len(comb)),
        "RMSE_m": rmse,
        "MAE_m": mae,
        "Bias_model_minus_obs_m": bias,
        "r": r,
        "NSE": nse,
    }
    return summary, comb.reset_index()


def identify_spring_neap(eta: np.ndarray):
    eta = np.asarray(eta, dtype=float)
    w = int(round(C["spring_neap_window_jam"] * 60 / C["interval_menit"]))
    w = max(w, 3)
    step = max(1, w // 2)

    mids = []
    amps = []
    for i in range(0, len(eta) - w + 1, step):
        seg = eta[i:i+w]
        if np.sum(np.isfinite(seg)) >= max(3, w // 2):
            amp = (np.nanmax(seg) - np.nanmin(seg)) / 2.0
            mids.append(i + w // 2)
            amps.append(float(amp))

    if not amps:
        return None

    amps_arr = np.asarray(amps, dtype=float)
    mids_arr = np.asarray(mids, dtype=int)
    th_s = float(np.percentile(amps_arr, C["spring_neap_percentile_spring"]))
    th_n = float(np.percentile(amps_arr, C["spring_neap_percentile_neap"]))

    return {
        "t_mid": mids_arr,
        "amplitudo": amps_arr,
        "threshold_spring": th_s,
        "threshold_neap": th_n,
        "idx_spring": mids_arr[amps_arr >= th_s],
        "idx_neap": mids_arr[amps_arr <= th_n],
    }
def hitung_formzahl_utide(
    eta: np.ndarray,
    waktu: pd.DatetimeIndex,
    latitude_deg: float,
):
    """
    Menghitung Formzahl dari time series elevasi muka air S1.

    Rumus:
        F = (A_K1 + A_O1) / (A_M2 + A_S2)

    Analisis harmonik dilakukan dengan UTide terhadap 4 komponen utama:
        M2, S2, K1, O1
    """

    if solve is None:
        raise ImportError(
            "Library UTide belum terpasang. Jalankan: pip install utide"
        )

    eta = np.asarray(eta, dtype=float)
    waktu = pd.to_datetime(waktu)

    # Hanya ambil data valid
    mask = np.isfinite(eta) & waktu.notna()
    eta_valid = eta[mask]
    waktu_valid = waktu[mask]

    if len(eta_valid) < 10:
        raise ValueError("Data S1 valid terlalu sedikit untuk analisis Formzahl.")

    durasi_hari = (
        waktu_valid.max() - waktu_valid.min()
    ).total_seconds() / 86400.0

    if durasi_hari < C["formzahl_min_durasi_hari"]:
        raise ValueError(
            f"Durasi data hanya {durasi_hari:.2f} hari. "
            f"Minimal disarankan {C['formzahl_min_durasi_hari']} hari."
        )

    # Analisis harmonik dengan UTide
    coef = solve(
        waktu_valid,
        eta_valid,
        lat=latitude_deg,
        constit=["M2", "S2", "K1", "O1"],
        nodal=False,
        trend=False,
        method="ols",
        conf_int="linear",
        Rayleigh_min=0.95,
        verbose=False,
    )

    # Normalisasi nama konstituen
    def norm_name(name):
        if isinstance(name, bytes):
            return name.decode(errors="ignore").strip().upper()
        return str(name).strip().upper()

    amplitudo = {
        norm_name(name): float(amp)
        for name, amp in zip(coef.name, coef.A)
    }

    A_M2 = amplitudo.get("M2", np.nan)
    A_S2 = amplitudo.get("S2", np.nan)
    A_K1 = amplitudo.get("K1", np.nan)
    A_O1 = amplitudo.get("O1", np.nan)

    if not all(np.isfinite(v) for v in [A_M2, A_S2, A_K1, A_O1]):
        raise ValueError(
            "Salah satu amplitudo M2, S2, K1, atau O1 tidak berhasil dihitung."
        )

    penyebut = A_M2 + A_S2
    if penyebut <= 0:
        raise ValueError("Penyebut Formzahl M2+S2 tidak valid.")

    formzahl = (A_K1 + A_O1) / penyebut

    # Klasifikasi tipe pasut
    if formzahl <= 0.25:
        tipe_pasut = "Harian ganda / Semidiurnal"
    elif formzahl <= 1.50:
        tipe_pasut = "Campuran condong harian ganda"
    elif formzahl <= 3.00:
        tipe_pasut = "Campuran condong harian tunggal"
    else:
        tipe_pasut = "Harian tunggal / Diurnal"

    return {
        "durasi_data_hari": float(durasi_hari),
        "jumlah_data_valid": int(len(eta_valid)),
        "latitude_deg": float(latitude_deg),

        "A_M2_m": float(A_M2),
        "A_S2_m": float(A_S2),
        "A_K1_m": float(A_K1),
        "A_O1_m": float(A_O1),

        "Formzahl_F": float(formzahl),
        "Tipe_Pasut": tipe_pasut,
    }

def gather_window_indices(mid_indices: np.ndarray, n_t: int, half_window_ts: int):
    s = set()
    for mid in mid_indices:
        for t in range(max(0, int(mid) - half_window_ts), min(n_t, int(mid) + half_window_ts + 1)):
            s.add(t)
    return sorted(s)

# =============================================================================
# OUTPUT GAMBAR
# =============================================================================

def save_maps_and_graphs(
    WAKTU,
    MAG,
    U_mean,
    V_mean,
    stats,
    POWER_mean,
    POWER_p95,
    valid_mask,
    layak_mask,
    optimal_mask,
    X,
    Y,
    DEPTH,
    enc_coords,
    bnd_segs,
    top_current,
    top_energy,
    quickplot_validation=None,
    rmse_validation=None,
    spring_neap_result=None,
    spring_neap_maps=None,
):
    dpi = int(C["dpi"])
    xlabel, ylabel = axes_labels(X, Y)

    section("BAGIAN 14: Membuat peta dan grafik")

    # 1. Peta V mean
    fig, ax = plt.subplots(figsize=(12, 9))
    vmax = float(np.nanpercentile(stats["v_mean"], 99)) if np.any(np.isfinite(stats["v_mean"])) else None
    art = plot_field(ax, X, Y, stats["v_mean"], C["colormap"], vmin=0, vmax=vmax)
    cb = plt.colorbar(art, ax=ax, pad=0.02)
    cb.set_label("Kecepatan Mean (m/s)")
    add_boundaries(ax, enc_coords, bnd_segs)
    # Quiver hanya pada koordinat finite
    if X is not None and Y is not None:
        skip = int(C["skip_quiver"])
        Xq = X[::skip, ::skip]
        Yq = Y[::skip, ::skip]
        Uq = U_mean[::skip, ::skip]
        Vq = V_mean[::skip, ::skip]
        qmask = np.isfinite(Xq) & np.isfinite(Yq) & np.isfinite(Uq) & np.isfinite(Vq)
        ax.quiver(Xq[qmask], Yq[qmask], Uq[qmask], Vq[qmask],
                  color="white", scale=C["skala_quiver"], width=0.002, alpha=0.75)
    ax.set_title(f"Distribusi Kecepatan Arus Mean — {C['nama_selat']}", fontweight="bold")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    plt.tight_layout()
    fp = SUBDIRS["peta"] / f"01_kecepatan_mean_{safe_name(C['kode_selat'])}.png"
    savefig_wm(fp, dpi=dpi, bbox_inches="tight"); plt.close()
    log(f"  ✓ {fp.name}")

    # 2. Peta frekuensi medium
    fig, ax = plt.subplots(figsize=(12, 9))
    art = plot_field(ax, X, Y, stats["freq_ge_medium_pct"], "viridis", vmin=0, vmax=100)
    cb = plt.colorbar(art, ax=ax, pad=0.02)
    cb.set_label(f"Frekuensi V ≥ {C['v_medium']} m/s (%)")
    add_boundaries(ax, enc_coords, bnd_segs)
    plot_top_points(ax, top_current, label_field="freq_ge_medium_pct")
    ax.set_title(f"Frekuensi Arus Kuat & Top Dominan Arus — {C['nama_selat']}", fontweight="bold")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    plt.tight_layout()
    fp = SUBDIRS["peta"] / f"02_frekuensi_medium_top_dominan_{safe_name(C['kode_selat'])}.png"
    savefig_wm(fp, dpi=dpi, bbox_inches="tight"); plt.close()
    log(f"  ✓ {fp.name}")

    # 3. Peta power density mean
    fig, ax = plt.subplots(figsize=(12, 9))
    vmax = float(np.nanpercentile(POWER_mean, 99)) if np.any(np.isfinite(POWER_mean)) else None
    art = plot_field(ax, X, Y, POWER_mean, "hot_r", vmin=0, vmax=vmax)
    cb = plt.colorbar(art, ax=ax, pad=0.02)
    cb.set_label("Mean Power Density (W/m²)")
    add_boundaries(ax, enc_coords, bnd_segs)
    plot_top_points(ax, top_energy, label_field="power_mean_wm2")
    ax.set_title(f"Mean Power Density & Top Potensi Energi — {C['nama_selat']}", fontweight="bold")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    plt.tight_layout()
    fp = SUBDIRS["peta"] / f"03_power_density_top_energi_{safe_name(C['kode_selat'])}.png"
    savefig_wm(fp, dpi=dpi, bbox_inches="tight"); plt.close()
    log(f"  ✓ {fp.name}")

    # 4. Peta kelayakan
    zona = np.zeros(valid_mask.shape, dtype=float)
    zona[valid_mask] = 0
    zona[layak_mask] = 1
    zona[optimal_mask] = 2
    zona[~valid_mask] = np.nan
    cmap = ListedColormap(["#d9d9d9", "#a8d5a2", "#2d6a4f"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    fig, ax = plt.subplots(figsize=(12, 9))
    art = plot_field(ax, X, Y, zona, cmap, vmin=0, vmax=2)
    add_boundaries(ax, enc_coords, bnd_segs)
    patches = [
        mpatches.Patch(color="#d9d9d9", label="Valid tetapi belum layak"),
        mpatches.Patch(color="#a8d5a2", label="Layak"),
        mpatches.Patch(color="#2d6a4f", label="Optimal"),
    ]
    ax.legend(handles=patches, loc="best", fontsize=9)
    ax.set_title(f"Peta Kelayakan Turbin — {C['nama_selat']}", fontweight="bold")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    plt.tight_layout()
    fp = SUBDIRS["peta"] / f"04_kelayakan_turbin_{safe_name(C['kode_selat'])}.png"
    savefig_wm(fp, dpi=dpi, bbox_inches="tight"); plt.close()
    log(f"  ✓ {fp.name}")

    # 5. Histogram sample
    sample = flatten_valid_sample(MAG, int(C["hist_max_sample"]))
    fig, ax = plt.subplots(figsize=(12, 6))
    if sample.size:
        ax.hist(sample, bins=80, alpha=0.85, edgecolor="white")
    for thr, label in [
        (C["v_cutin"], f"Cut-in {C['v_cutin']} m/s"),
        (C["v_medium"], f"Medium {C['v_medium']} m/s"),
        (C["v_tinggi"], f"Tinggi {C['v_tinggi']} m/s"),
    ]:
        ax.axvline(thr, ls="--", lw=1.5, label=label)
    ax.set_xlabel("Kecepatan Arus (m/s)")
    ax.set_ylabel("Frekuensi sampel")
    ax.set_title(f"Histogram Kecepatan Arus — {C['nama_selat']}", fontweight="bold")
    ax.legend(); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fp = SUBDIRS["grafik"] / f"05_histogram_kecepatan_{safe_name(C['kode_selat'])}.png"
    savefig_wm(fp, dpi=dpi, bbox_inches="tight"); plt.close()
    log(f"  ✓ {fp.name}")

    # 6. Timeseries P1 dominan arus
    if top_current:
        p = top_current[0]
        ts = p["_ts_mag"]
        fig, ax = plt.subplots(figsize=(16, 6))
        ax.plot(WAKTU, ts, lw=0.9, label="P1 Dominan Arus")
        for thr, label in [
            (C["v_cutin"], f"{C['v_cutin']} m/s"),
            (C["v_medium"], f"{C['v_medium']} m/s"),
            (C["v_tinggi"], f"{C['v_tinggi']} m/s"),
        ]:
            ax.axhline(thr, ls="--", lw=1.2, label=label)
        ax.set_title(
            f"Time Series P1 Dominan Arus M={p['M']}, N={p['N']} — {C['nama_selat']}\n"
            f"Freq ≥ {C['v_medium']} m/s = {p['freq_ge_medium_pct']:.2f}% | Vmean={p['v_mean']:.3f} m/s",
            fontweight="bold"
        )
        ax.set_xlabel("Waktu"); ax.set_ylabel("Kecepatan (m/s)")
        ax.legend(ncol=4, fontsize=9); ax.grid(True, alpha=0.3)
        fig.autofmt_xdate()
        plt.tight_layout()
        fp = SUBDIRS["grafik"] / f"06_timeseries_P1_dominan_arus_{safe_name(C['kode_selat'])}.png"
        savefig_wm(fp, dpi=dpi, bbox_inches="tight"); plt.close()
        log(f"  ✓ {fp.name}")

    # 7. Timeseries P1 energi
    if top_energy:
        p = top_energy[0]
        ts = p["_ts_mag"]
        pwr = p["_ts_pwr"]
        fig, axes = plt.subplots(2, 1, figsize=(16, 9), sharex=True)
        axes[0].plot(WAKTU, ts, lw=0.9)
        axes[0].set_ylabel("V (m/s)")
        axes[0].set_title(
            f"P1 Potensi Energi M={p['M']}, N={p['N']} | Vmean={p['v_mean']:.3f} m/s",
            fontweight="bold"
        )
        axes[0].grid(True, alpha=0.3)
        axes[1].plot(WAKTU, pwr, lw=0.9)
        axes[1].set_ylabel("Power Density (W/m²)")
        axes[1].set_xlabel("Waktu")
        axes[1].set_title(f"Mean Power={p['power_mean_wm2']:.2f} W/m²", fontweight="bold")
        axes[1].grid(True, alpha=0.3)
        fig.autofmt_xdate()
        plt.tight_layout()
        fp = SUBDIRS["grafik"] / f"07_timeseries_P1_potensi_energi_{safe_name(C['kode_selat'])}.png"
        savefig_wm(fp, dpi=dpi, bbox_inches="tight"); plt.close()
        log(f"  ✓ {fp.name}")

    # 8. Validasi QuickPlot
    if quickplot_validation is not None:
        summary, merged = quickplot_validation
        fig, ax = plt.subplots(figsize=(16, 6))
        ax.plot(merged["waktu"], merged["velocity_quickplot"], lw=0.9, label="QuickPlot")
        ax.plot(merged["waktu"], merged["velocity_model"], lw=0.9, label="Model Script Final")
        ax.set_title(
            f"Validasi Velocity vs QuickPlot M={summary['M']}, N={summary['N']} | "
            f"RMSE={summary['RMSE_mps']:.8f} m/s | r={summary['r']:.6f}",
            fontweight="bold"
        )
        ax.set_xlabel("Waktu"); ax.set_ylabel("Velocity (m/s)")
        ax.legend(); ax.grid(True, alpha=0.3)
        fig.autofmt_xdate()
        plt.tight_layout()
        fp = SUBDIRS["grafik"] / f"08_validasi_velocity_quickplot_{safe_name(C['kode_selat'])}.png"
        savefig_wm(fp, dpi=dpi, bbox_inches="tight"); plt.close()
        log(f"  ✓ {fp.name}")

    # 9. Validasi RMSE S1
    if rmse_validation is not None:
        summary, comb = rmse_validation
        fig, axes = plt.subplots(2, 1, figsize=(16, 10))
        axes[0].plot(comb["waktu"], comb["eta_model"], label="S1 Model", lw=0.9)
        axes[0].plot(comb["waktu"], comb["eta_obs"], label="Observasi", lw=0.9)
        axes[0].set_ylabel("Elevasi muka air (m)")
        axes[0].set_title(
            f"Validasi S1 | RMSE={summary['RMSE_m']:.4f} m | r={summary['r']:.4f} | NSE={summary['NSE']:.4f}",
            fontweight="bold"
        )
        axes[0].legend(); axes[0].grid(True, alpha=0.3)
        axes[1].scatter(comb["eta_obs"], comb["eta_model"], s=10, alpha=0.6)
        mn = float(min(comb["eta_obs"].min(), comb["eta_model"].min()))
        mx = float(max(comb["eta_obs"].max(), comb["eta_model"].max()))
        axes[1].plot([mn, mx], [mn, mx], ls="--", lw=1.2)
        axes[1].set_xlabel("Observasi (m)")
        axes[1].set_ylabel("Model (m)")
        axes[1].set_title("Scatter Model vs Observasi", fontweight="bold")
        axes[1].grid(True, alpha=0.3)
        fig.autofmt_xdate()
        plt.tight_layout()
        fp = SUBDIRS["grafik"] / f"09_validasi_RMSE_S1_{safe_name(C['kode_selat'])}.png"
        savefig_wm(fp, dpi=dpi, bbox_inches="tight"); plt.close()
        log(f"  ✓ {fp.name}")

    # 10. Spring/Neap
    if spring_neap_result is not None:
        sn = spring_neap_result
        tm_days = sn["t_mid"] * C["interval_menit"] / (24 * 60)
        fig, ax = plt.subplots(figsize=(16, 5))
        ax.plot(tm_days, sn["amplitudo"], lw=1.5, label="Amplitudo pasut lokal")
        ax.axhline(sn["threshold_spring"], ls="--", lw=1.2,
                   label=f"Spring P{C['spring_neap_percentile_spring']} = {sn['threshold_spring']:.3f} m")
        ax.axhline(sn["threshold_neap"], ls="--", lw=1.2,
                   label=f"Neap P{C['spring_neap_percentile_neap']} = {sn['threshold_neap']:.3f} m")
        ax.fill_between(tm_days, sn["amplitudo"], sn["threshold_spring"],
                        where=sn["amplitudo"] >= sn["threshold_spring"], alpha=0.25,
                        label="Periode Spring")
        ax.fill_between(tm_days, sn["amplitudo"], sn["threshold_neap"],
                        where=sn["amplitudo"] <= sn["threshold_neap"], alpha=0.25,
                        label="Periode Neap")
        ax.set_xlabel("Hari ke- sejak awal analisis")
        ax.set_ylabel("Amplitudo pasut (m)")
        ax.set_title(f"Identifikasi Spring/Neap — {C['nama_selat']}", fontweight="bold")
        ax.legend(ncol=2, fontsize=9); ax.grid(True, alpha=0.3)
        plt.tight_layout()
        fp = SUBDIRS["grafik"] / f"10_identifikasi_spring_neap_{safe_name(C['kode_selat'])}.png"
        savefig_wm(fp, dpi=dpi, bbox_inches="tight"); plt.close()
        log(f"  ✓ {fp.name}")

        if spring_neap_maps is not None:
            map_s, map_n, map_mean = spring_neap_maps
            vmax = float(np.nanpercentile(map_s, 99)) if np.any(np.isfinite(map_s)) else None
            fig, axes = plt.subplots(1, 3, figsize=(20, 7))
            for ax, data, title in zip(
                axes,
                [map_s, map_n, map_mean],
                ["Spring Tide", "Neap Tide", "Mean Keseluruhan"],
            ):
                art = plot_field(ax, X, Y, data, C["colormap"], vmin=0, vmax=vmax)
                add_boundaries(ax, enc_coords, bnd_segs)
                ax.set_title(title, fontweight="bold")
                ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
            fig.colorbar(art, ax=axes.ravel().tolist(), pad=0.02, label="Kecepatan Mean (m/s)")
            plt.suptitle(f"Perbandingan Spring / Neap — {C['nama_selat']}", fontweight="bold")
            plt.tight_layout()
            fp = SUBDIRS["peta"] / f"05_peta_spring_neap_{safe_name(C['kode_selat'])}.png"
            savefig_wm(fp, dpi=dpi, bbox_inches="tight"); plt.close()
            log(f"  ✓ {fp.name}")


# =============================================================================
# GIF OPSIONAL
# =============================================================================

def create_gif_daily(WAKTU, MAG, X, Y, enc_coords, bnd_segs):
    if not C["buat_gif"]:
        return
    if imageio is None:
        log("  ⚠ GIF dilewati: imageio tidak tersedia.")
        return

    section("BAGIAN 15: Membuat GIF harian")
    stride = max(1, int(round(C["gif_frame_interval_hari"] * 24 * 60 / C["interval_menit"])))
    frame_indices = list(range(0, MAG.shape[0], stride))
    out_gif = SUBDIRS["gif"] / f"animasi_arus_harian_{safe_name(C['kode_selat'])}.gif"
    tmp_dir = SUBDIRS["gif"] / "_frames_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    frame_paths = []

    vmax = float(np.nanpercentile(MAG, 99)) if np.any(np.isfinite(MAG)) else None
    xlabel, ylabel = axes_labels(X, Y)

    for i, t in enumerate(frame_indices, start=1):
        fig, ax = plt.subplots(figsize=(11, 8))
        art = plot_field(ax, X, Y, MAG[t], C["colormap"], vmin=0, vmax=vmax)
        plt.colorbar(art, ax=ax, pad=0.02, label="Kecepatan (m/s)")
        add_boundaries(ax, enc_coords, bnd_segs)
        ax.set_title(f"{C['nama_selat']} | {WAKTU[t]}", fontweight="bold")
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
        plt.tight_layout()
        fp = tmp_dir / f"frame_{i:04d}.png"
        savefig_wm(fp, dpi=140, bbox_inches="tight")
        plt.close(fig)
        frame_paths.append(fp)
        log(f"    Frame {i}/{len(frame_indices)} ✓")

    with imageio.get_writer(out_gif, mode="I", fps=int(C["gif_fps"])) as writer:
        for fp in frame_paths:
            writer.append_data(imageio.imread(fp))

    for fp in frame_paths:
        try:
            fp.unlink()
        except Exception:
            pass
    try:
        tmp_dir.rmdir()
    except Exception:
        pass

    log(f"  ✓ GIF tersimpan: {out_gif}")


# =============================================================================
# EKSPOR CSV & RINGKASAN
# =============================================================================

def export_outputs(
    domain_summary: dict,
    grid_rows: list[dict],
    top_current: list[dict],
    top_energy: list[dict],
    spike_rows: list[dict],
    quickplot_validation=None,
    rmse_validation=None,
    spring_neap_result=None,
):
    section("BAGIAN 16: Ekspor CSV")

    pd.DataFrame([domain_summary]).to_csv(
        SUBDIRS["csv"] / f"statistik_domain_{safe_name(C['kode_selat'])}.csv",
        index=False,
        encoding="utf-8-sig",
    )
    pd.DataFrame(grid_rows).to_csv(
        SUBDIRS["csv"] / f"statistik_grid_{safe_name(C['kode_selat'])}.csv",
        index=False,
        encoding="utf-8-sig",
    )
    pd.DataFrame(public_point_rows(top_current)).to_csv(
        SUBDIRS["csv"] / f"top5_dominan_arus_{safe_name(C['kode_selat'])}.csv",
        index=False,
        encoding="utf-8-sig",
    )
    pd.DataFrame(public_point_rows(top_energy)).to_csv(
        SUBDIRS["csv"] / f"top5_potensi_energi_{safe_name(C['kode_selat'])}.csv",
        index=False,
        encoding="utf-8-sig",
    )
    pd.DataFrame(spike_rows).to_csv(
        SUBDIRS["csv"] / f"audit_spike_kecepatan_{safe_name(C['kode_selat'])}.csv",
        index=False,
        encoding="utf-8-sig",
    )

    if top_current:
        pd.DataFrame({
            "waktu": domain_summary["_WAKTU"],
            "velocity_P1_dominan_arus_mps": top_current[0]["_ts_mag"],
            "power_P1_dominan_arus_wm2": top_current[0]["_ts_pwr"],
        }).to_csv(
            SUBDIRS["csv"] / f"timeseries_P1_dominan_arus_{safe_name(C['kode_selat'])}.csv",
            index=False,
            encoding="utf-8-sig",
        )

    if top_energy:
        pd.DataFrame({
            "waktu": domain_summary["_WAKTU"],
            "velocity_P1_potensi_energi_mps": top_energy[0]["_ts_mag"],
            "power_P1_potensi_energi_wm2": top_energy[0]["_ts_pwr"],
        }).to_csv(
            SUBDIRS["csv"] / f"timeseries_P1_potensi_energi_{safe_name(C['kode_selat'])}.csv",
            index=False,
            encoding="utf-8-sig",
        )

    if quickplot_validation is not None:
        summary, merged = quickplot_validation
        pd.DataFrame([summary]).to_csv(
            SUBDIRS["csv"] / f"validasi_velocity_quickplot_ringkas_{safe_name(C['kode_selat'])}.csv",
            index=False,
            encoding="utf-8-sig",
        )
        merged.to_csv(
            SUBDIRS["csv"] / f"validasi_velocity_quickplot_series_{safe_name(C['kode_selat'])}.csv",
            index=False,
            encoding="utf-8-sig",
        )

    if rmse_validation is not None:
        summary, comb = rmse_validation
        pd.DataFrame([summary]).to_csv(
            SUBDIRS["csv"] / f"validasi_RMSE_S1_ringkas_{safe_name(C['kode_selat'])}.csv",
            index=False,
            encoding="utf-8-sig",
        )
        comb.to_csv(
            SUBDIRS["csv"] / f"validasi_RMSE_S1_series_{safe_name(C['kode_selat'])}.csv",
            index=False,
            encoding="utf-8-sig",
        )

    if spring_neap_result is not None:
        sn = spring_neap_result
        pd.DataFrame({
            "t_mid_index": sn["t_mid"],
            "amplitudo_m": sn["amplitudo"],
        }).to_csv(
            SUBDIRS["csv"] / f"spring_neap_amplitudo_{safe_name(C['kode_selat'])}.csv",
            index=False,
            encoding="utf-8-sig",
        )

    log("  ✓ Semua CSV selesai diekspor.")


def write_summary_txt(
    domain_summary: dict,
    top_current: list[dict],
    top_energy: list[dict],
    depth_note: str,
    quickplot_validation=None,
    rmse_validation=None,
    spring_neap_result=None,
    formzahl_result=None,
):
    section("BAGIAN 17: Menulis ringkasan TXT")
    out = OUTPUT_DIR / f"RINGKASAN_ANALISIS_FINAL_{safe_name(C['kode_selat'])}.txt"

    L: list[str] = []
    def a(s=""):
        L.append(str(s))

    a(WATERMARK_TEXT)
    a("=" * 90)
    a(f"RINGKASAN ANALISIS FINAL — {C['nama_selat']}")
    a("=" * 90)
    a(f"Metode velocity final:")
    a("  - Getelt map-series TIME-only")
    a("  - Reshape order='F'")
    a("  - BACKWARD_Ubackward_Vbackward")
    a("")
    a("RINGKASAN DATA")
    a("-" * 90)
    for k in [
        "n_timestep_total",
        "n_timestep_spinup",
        "n_timestep_analisis",
        "n_sel_valid",
        "n_sel_layak",
        "n_sel_optimal",
        "v_min_mps",
        "v_max_mps",
        "v_mean_mps",
        "v_median_mps",
        "v_p90_mps",
        "v_p95_mps",
        "power_mean_wm2",
        "power_p95_wm2",
        "power_max_wm2",
    ]:
        a(f"{k:<28}: {domain_summary.get(k)}")
    a(f"Catatan depth              : {depth_note}")
    a("")

    a("TOP 5 TITIK DOMINAN ARUS")
    a("-" * 90)
    a("Definisi: frekuensi V >= v_medium tertinggi, lalu V_mean, lalu V_P95.")
    for p in top_current:
        a(
            f"P{p['no']} M={p['M']}, N={p['N']} | "
            f"X={p['x']:.6f}, Y={p['y']:.6f} | "
            f"Depth={p['depth_m']:.2f} m | "
            f"Freq≥{C['v_medium']}={p['freq_ge_medium_pct']:.2f}% | "
            f"Vmean={p['v_mean']:.4f} m/s | Vp95={p['v_p95']:.4f} m/s | "
            f"PowerMean={p['power_mean_wm2']:.2f} W/m² | Layak={p['layak_zona']}"
        )
    a("")

    a("TOP 5 POTENSI ENERGI")
    a("-" * 90)
    a("Definisi: mean power density tertinggi.")
    for p in top_energy:
        a(
            f"P{p['no']} M={p['M']}, N={p['N']} | "
            f"X={p['x']:.6f}, Y={p['y']:.6f} | "
            f"Depth={p['depth_m']:.2f} m | "
            f"PowerMean={p['power_mean_wm2']:.2f} W/m² | "
            f"Vmean={p['v_mean']:.4f} m/s | Vp95={p['v_p95']:.4f} m/s | "
            f"Layak={p['layak_zona']}"
        )
    a("")

    if quickplot_validation is not None:
        summary, _ = quickplot_validation
        a("VALIDASI VELOCITY vs QUICKPLOT")
        a("-" * 90)
        for k, v in summary.items():
            a(f"{k:<40}: {v}")
        a("")

    if rmse_validation is not None:
        summary, _ = rmse_validation
        a("VALIDASI RMSE S1")
        a("-" * 90)
        for k, v in summary.items():
            a(f"{k:<40}: {v}")
        a("")

    if spring_neap_result is not None:
        sn = spring_neap_result
        a("ANALISIS SPRING/NEAP")
        a("-" * 90)
        a(f"Threshold spring P{C['spring_neap_percentile_spring']}: {sn['threshold_spring']:.4f} m")
        a(f"Threshold neap P{C['spring_neap_percentile_neap']}: {sn['threshold_neap']:.4f} m")
        a(f"Jumlah pusat window spring: {len(sn['idx_spring'])}")
        a(f"Jumlah pusat window neap  : {len(sn['idx_neap'])}")
        a("")

    a("OUTPUT")
    a("-" * 90)
    a(f"Folder output: {OUTPUT_DIR}")
    a("Subfolder: peta/, grafik/, csv/, animasi/")
    a("=" * 90)

    out.write_text("\n".join(L), encoding="utf-8")
    log(f"  ✓ {out.name}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    for _key in ("trim_dat", "trim_def"):
        if not str(C[_key]).strip():
            raise SystemExit(f"CONFIG['{_key}'] masih kosong. Isi path file trim di bagian KONFIGURASI.")
    t_global = datetime.now()
    log("=" * 78)
    log("  ANALISIS UTAMA FINAL — POTENSI ENERGI & TITIK DOMINAN ARUS")
    log(f"  {WATERMARK_TEXT}")
    log(f"  {C['nama_selat'].upper()}")
    log("=" * 78)
    log(f"  Mulai  : {t_global.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"  Output : {OUTPUT_DIR}")

    # -------------------------------------------------------------------------
    # 1. Baca GRD lebih awal agar diketahui shape analisis
    # -------------------------------------------------------------------------
    section("BAGIAN 1: Parser GRD / DEP / ENC / BND")
    X_grd, Y_grd = parse_grd_file(C["grd_file"])

    # -------------------------------------------------------------------------
    # 2. Buka NEFIS dan baca dimensi
    # -------------------------------------------------------------------------
    section("BAGIAN 2: NEFIS DLL dan dimensi trim")
    nef = load_nefis()
    fd = nefis_open(nef, C["trim_dat"], C["trim_def"])
    try:
        mmax_full = read_int_scalar(nef, fd, "map-const", "MMAX")
        nmax_full = read_int_scalar(nef, fd, "map-const", "NMAX")
        try:
            kmax = read_int_scalar(nef, fd, "map-const", "KMAX")
        except Exception:
            kmax = 1
        log(f"  ✓ MMAX={mmax_full}, NMAX={nmax_full}, KMAX={kmax}")

        # Shape analisis: jika GRD tersedia dan ukurannya <= trim, gunakan GRD.
        if X_grd is not None and Y_grd is not None:
            an_n, an_m = X_grd.shape
            if an_n <= nmax_full and an_m <= mmax_full:
                analysis_shape = (an_n, an_m)
                X, Y = X_grd, Y_grd
                if analysis_shape != (nmax_full, mmax_full):
                    log(
                        f"  ✓ Sinkronisasi shape: trim ({nmax_full},{mmax_full}) "
                        f"→ GRD {analysis_shape} dengan crop baris/kolom akhir."
                    )
            else:
                log("  ⚠ Shape GRD lebih besar dari trim; fallback ke indeks grid trim.")
                analysis_shape = (nmax_full, mmax_full)
                XX, YY = np.meshgrid(np.arange(1, mmax_full+1), np.arange(1, nmax_full+1))
                X, Y = XX.astype(float), YY.astype(float)
        else:
            analysis_shape = (nmax_full, mmax_full)
            XX, YY = np.meshgrid(np.arange(1, mmax_full+1), np.arange(1, nmax_full+1))
            X, Y = XX.astype(float), YY.astype(float)
            log("  ⚠ GRD tidak tersedia; koordinat memakai indeks grid.")

        DEPTH = parse_dep_file(C["dep_file"], analysis_shape)
        enc_coords = parse_enc_file(C["enc_file"], X, Y)
        bnd_segs = parse_bnd_file(C["bnd_file"], X, Y)

        # ---------------------------------------------------------------------
        # 3. Hitung timestep
        # ---------------------------------------------------------------------
        count_uv = nmax_full * mmax_full * kmax
        n_time = count_timesteps(nef, fd, count_uv)
        if n_time <= 0:
            raise RuntimeError("Tidak ada timestep U1/V1 yang terbaca.")

        n_spinup = int(round(C["spinup_hari"] * 24 * 60 / C["interval_menit"]))
        if n_spinup >= n_time:
            raise ValueError(
                f"Spin-up {n_spinup} timestep >= total timestep {n_time}. "
                "Kurangi spinup_hari atau periksa interval."
            )
        n_analysis = n_time - n_spinup
        log(f"  ✓ Spin-up: {C['spinup_hari']} hari = {n_spinup} timestep")
        log(f"  ✓ Timestep analisis: {n_analysis}")

        # ---------------------------------------------------------------------
        # 4. Baca U/V dan susun MAG_all secara streaming
        # ---------------------------------------------------------------------
        section("BAGIAN 4: Membaca U1/V1 TIME-only dan merekonstruksi velocity")
        an_n, an_m = analysis_shape
        MAG = np.full((n_analysis, an_n, an_m), np.nan, dtype=np.float32)

        sum_u = np.zeros((an_n, an_m), dtype=np.float64)
        sum_v = np.zeros((an_n, an_m), dtype=np.float64)
        count_uv_map = np.zeros((an_n, an_m), dtype=np.int64)
        sum_power = np.zeros((an_n, an_m), dtype=np.float64)
        count_power = np.zeros((an_n, an_m), dtype=np.int64)

        spike_rows = []
        spike_total = 0

        analysis_idx = 0
        for t in range(n_time):
            U_flat = read_series_element_time_only(nef, fd, "U1", t, count_uv)
            V_flat = read_series_element_time_only(nef, fd, "V1", t, count_uv)
            if t < n_spinup:
                continue

            U = clean_values(reshape_uv_f_order(U_flat, nmax_full, mmax_full, kmax))
            V = clean_values(reshape_uv_f_order(V_flat, nmax_full, mmax_full, kmax))
            Uc, Vc, Mag_full = reconstruct_backward_velocity(U, V)

            Uc = Uc[:an_n, :an_m]
            Vc = Vc[:an_n, :an_m]
            Mag = Mag_full[:an_n, :an_m]

            MAG[analysis_idx] = Mag

            finite_uv = np.isfinite(Uc) & np.isfinite(Vc)
            sum_u[finite_uv] += Uc[finite_uv]
            sum_v[finite_uv] += Vc[finite_uv]
            count_uv_map[finite_uv] += 1

            power = 0.5 * C["rho"] * Mag.astype(np.float64)**3
            finite_power = np.isfinite(power)
            sum_power[finite_power] += power[finite_power]
            count_power[finite_power] += 1

            spikes = np.argwhere(np.isfinite(Mag) & (Mag > C["spike_alert_mps"]))
            spike_total += int(len(spikes))
            if spikes.size and len(spike_rows) < int(C["max_spike_rows_csv"]):
                waktu = (
                    pd.Timestamp(C["tanggal_mulai"])
                    + pd.to_timedelta(t * C["interval_menit"], unit="m")
                )
                for n0, m0 in spikes:
                    if len(spike_rows) >= int(C["max_spike_rows_csv"]):
                        break
                    spike_rows.append({
                        "timestep_asli_1based": int(t + 1),
                        "timestep_analisis_1based": int(analysis_idx + 1),
                        "waktu": waktu,
                        "M": int(m0 + 1),
                        "N": int(n0 + 1),
                        "velocity_mps": float(Mag[n0, m0]),
                    })

            analysis_idx += 1
            if analysis_idx % 250 == 0 or analysis_idx == n_analysis:
                log(f"    {analysis_idx}/{n_analysis} timestep analisis ✓")

        if spike_total > 0:
            log(
                f"  ⚠ Audit spike: {spike_total} nilai > {C['spike_alert_mps']} m/s. "
                f"CSV menyimpan maksimal {C['max_spike_rows_csv']} baris."
            )
        else:
            log(f"  ✓ Tidak ada spike > {C['spike_alert_mps']} m/s.")

    finally:
        nefis_close(nef, fd)

    # -------------------------------------------------------------------------
    # 5. Waktu analisis dan masking akhir
    # -------------------------------------------------------------------------
    section("BAGIAN 5: Masking akhir dan statistik peta")
    WAKTU = create_model_time(n_analysis)

    finite_any = np.any(np.isfinite(MAG), axis=0)
    max_map = np.nanmax(MAG, axis=0)
    dynamic_valid = finite_any & np.isfinite(max_map) & (max_map > 0)

    coord_valid = np.isfinite(X) & np.isfinite(Y) if X is not None and Y is not None else np.ones_like(dynamic_valid, dtype=bool)
    if DEPTH is not None and np.any(np.isfinite(DEPTH)):
        depth_valid = np.isfinite(DEPTH)
    else:
        depth_valid = np.ones_like(dynamic_valid, dtype=bool)

    valid_mask = dynamic_valid & coord_valid & depth_valid
    MAG[:, ~valid_mask] = np.nan

    U_mean = np.divide(sum_u, count_uv_map, out=np.full_like(sum_u, np.nan), where=count_uv_map > 0)
    V_mean = np.divide(sum_v, count_uv_map, out=np.full_like(sum_v, np.nan), where=count_uv_map > 0)
    U_mean[~valid_mask] = np.nan
    V_mean[~valid_mask] = np.nan

    POWER_mean = np.divide(sum_power, count_power, out=np.full_like(sum_power, np.nan), where=count_power > 0)
    POWER_mean[~valid_mask] = np.nan

    MAG_p95 = np.nanpercentile(MAG, 95, axis=0)
    POWER_p95 = 0.5 * C["rho"] * MAG_p95**3
    POWER_p95[~valid_mask] = np.nan

    log(f"  ✓ Sel valid analisis: {int(np.sum(valid_mask))}")
    log(f"  ✓ Sel coordinate invalid: {int(np.sum(~coord_valid))}")
    if DEPTH is not None:
        log(f"  ✓ Sel depth invalid: {int(np.sum(~depth_valid))}")

    # -------------------------------------------------------------------------
    # 6. Statistik domain dan grid
    # -------------------------------------------------------------------------
    section("BAGIAN 6: Statistik domain & power density")

    stats, grid_rows = build_grid_statistics(
        MAG, POWER_mean, POWER_p95, DEPTH, X, Y, valid_mask
    )
    layak_mask, optimal_mask, depth_note = define_layak_mask(valid_mask, stats["v_p95"], DEPTH)

    valid_total = int(np.sum(np.isfinite(MAG)))
    prob_cutin = float(np.sum(MAG >= C["v_cutin"]) / max(valid_total, 1) * 100.0)
    prob_medium = float(np.sum(MAG >= C["v_medium"]) / max(valid_total, 1) * 100.0)
    prob_tinggi = float(np.sum(MAG >= C["v_tinggi"]) / max(valid_total, 1) * 100.0)

    v_eff_cutin = float(np.nanmean(np.where(MAG >= C["v_cutin"], MAG, np.nan)))
    v_eff_medium = float(np.nanmean(np.where(MAG >= C["v_medium"], MAG, np.nan)))
    v_eff_tinggi = float(np.nanmean(np.where(MAG >= C["v_tinggi"], MAG, np.nan)))

    v_global_p95 = float(np.nanpercentile(MAG, 95))
    power_global_mean = float(np.nansum(sum_power[valid_mask]) / max(int(np.nansum(count_power[valid_mask])), 1))
    power_global_p95 = float(0.5 * C["rho"] * v_global_p95**3)
    power_global_max = float(0.5 * C["rho"] * float(np.nanmax(MAG))**3)

    domain_summary = {
        "nama_selat": C["nama_selat"],
        "kode_selat": C["kode_selat"],
        "tanggal_awal_analisis": str(WAKTU[0]) if len(WAKTU) else "",
        "tanggal_akhir_analisis": str(WAKTU[-1]) if len(WAKTU) else "",
        "interval_menit": C["interval_menit"],
        "n_timestep_total": int(n_time),
        "n_timestep_spinup": int(n_spinup),
        "n_timestep_analisis": int(n_analysis),
        "n_sel_valid": int(np.sum(valid_mask)),
        "n_sel_layak": int(np.sum(layak_mask)),
        "n_sel_optimal": int(np.sum(optimal_mask)),
        "v_min_mps": float(np.nanmin(MAG)),
        "v_max_mps": float(np.nanmax(MAG)),
        "v_mean_mps": float(np.nanmean(MAG)),
        "v_median_mps": float(np.nanpercentile(MAG, 50)),
        "v_p90_mps": float(np.nanpercentile(MAG, 90)),
        "v_p95_mps": v_global_p95,
        "v_std_mps": float(np.nanstd(MAG)),
        "prob_ge_cutin_pct": prob_cutin,
        "prob_ge_medium_pct": prob_medium,
        "prob_ge_tinggi_pct": prob_tinggi,
        "v_eff_ge_cutin_mps": v_eff_cutin,
        "v_eff_ge_medium_mps": v_eff_medium,
        "v_eff_ge_tinggi_mps": v_eff_tinggi,
        "power_mean_wm2": power_global_mean,
        "power_p95_wm2": power_global_p95,
        "power_max_wm2": power_global_max,
        "spike_total_count": int(spike_total),
        "_WAKTU": WAKTU,
    }

    for k in [
        "v_min_mps", "v_max_mps", "v_mean_mps", "v_median_mps",
        "v_p90_mps", "v_p95_mps", "power_mean_wm2", "power_p95_wm2",
        "power_max_wm2"
    ]:
        log(f"  {k:<24}: {domain_summary[k]:.6f}")
    log(f"  Prob V ≥ {C['v_cutin']:.2f} m/s   : {prob_cutin:.2f}%")
    log(f"  Prob V ≥ {C['v_medium']:.2f} m/s  : {prob_medium:.2f}%")
    log(f"  Prob V ≥ {C['v_tinggi']:.2f} m/s  : {prob_tinggi:.2f}%")
    log(f"  Zona layak / optimal      : {int(np.sum(layak_mask))} / {int(np.sum(optimal_mask))} sel")
    log(f"  {depth_note}")

    # -------------------------------------------------------------------------
    # 7. Ranking top points
    # -------------------------------------------------------------------------
    section("BAGIAN 7: Ranking Top 5")
    top_current = rank_top_dominant_current(MAG, valid_mask, layak_mask, stats, DEPTH, X, Y)
    top_energy = rank_top_energy(MAG, valid_mask, layak_mask, stats, DEPTH, X, Y)

    log("  TOP DOMINAN ARUS:")
    for p in top_current:
        log(
            f"    P{p['no']} M={p['M']},N={p['N']} | "
            f"Freq≥{C['v_medium']}={p['freq_ge_medium_pct']:.2f}% | "
            f"Vmean={p['v_mean']:.4f} | PowerMean={p['power_mean_wm2']:.2f}"
        )

    log("  TOP POTENSI ENERGI:")
    for p in top_energy:
        log(
            f"    P{p['no']} M={p['M']},N={p['N']} | "
            f"PowerMean={p['power_mean_wm2']:.2f} | "
            f"Vmean={p['v_mean']:.4f} | Freq≥{C['v_medium']}={p['freq_ge_medium_pct']:.2f}%"
        )

    # -------------------------------------------------------------------------
    # 8. Validasi QuickPlot
    # -------------------------------------------------------------------------
    quickplot_validation = None
    if C["validasi_quickplot_aktif"]:
        section("BAGIAN 8: Validasi velocity vs QuickPlot")
        try:
            quickplot_validation = compare_model_quickplot_velocity(
                MAG, WAKTU,
                int(C["quickplot_M"]),
                int(C["quickplot_N"]),
                C["quickplot_velocity_csv"],
            )
            qsum, _ = quickplot_validation
            log(
                f"  ✓ QuickPlot M={qsum['M']}, N={qsum['N']} | "
                f"RMSE={qsum['RMSE_mps']:.8f} m/s | r={qsum['r']:.8f}"
            )
        except Exception as exc:
            log(f"  ⚠ Validasi QuickPlot gagal/dilewati: {exc}")

    # -------------------------------------------------------------------------
    # 9. S1 extractions
    # -------------------------------------------------------------------------
    requests = {}
    try:
        requests = get_s1_request_points(top_current, top_energy)
    except Exception as exc:
        log(f"  ⚠ Penyusunan request S1 gagal: {exc}")

    s1_series = {}
    if requests:
        nef = load_nefis()
        fd = nefis_open(nef, C["trim_dat"], C["trim_def"])
        try:
            s1_series = extract_s1_timeseries(nef, fd, n_time, n_spinup, nmax_full, mmax_full, requests)
        finally:
            nefis_close(nef, fd)

    # -------------------------------------------------------------------------
    # 10. RMSE S1
    # -------------------------------------------------------------------------
    rmse_validation = None
    if C["validasi_rmse_s1_aktif"]:
        section("BAGIAN 10: Validasi RMSE S1")
        try:
            if "rmse_obs" not in s1_series:
                raise ValueError("S1 series untuk RMSE tidak tersedia.")
            rmse_validation = rmse_s1_validation(s1_series["rmse_obs"], WAKTU, C["obs_s1_csv"])
            rsum, _ = rmse_validation
            log(
                f"  ✓ RMSE S1={rsum['RMSE_m']:.6f} m | "
                f"MAE={rsum['MAE_m']:.6f} m | r={rsum['r']:.6f} | NSE={rsum['NSE']:.6f}"
            )
        except Exception as exc:
            log(f"  ⚠ Validasi RMSE S1 gagal/dilewati: {exc}")

    # -------------------------------------------------------------------------
    # 11. Spring/Neap
    # -------------------------------------------------------------------------
    spring_neap_result = None
    spring_neap_maps = None
    if C["analisis_spring_neap_aktif"]:
        section("BAGIAN 11: Analisis Spring/Neap")
        try:
            if "spring_neap" not in s1_series:
                raise ValueError("S1 referensi Spring/Neap tidak tersedia.")
            spring_neap_result = identify_spring_neap(s1_series["spring_neap"])
            if spring_neap_result is None:
                raise ValueError("Amplitudo Spring/Neap tidak dapat dihitung.")
            sn = spring_neap_result
            half_window_ts = int(round(12 * 60 / C["interval_menit"]))
            idx_s = gather_window_indices(sn["idx_spring"], MAG.shape[0], half_window_ts)
            idx_n = gather_window_indices(sn["idx_neap"], MAG.shape[0], half_window_ts)
            map_s = np.nanmean(MAG[idx_s], axis=0) if idx_s else np.nanmean(MAG, axis=0)
            map_n = np.nanmean(MAG[idx_n], axis=0) if idx_n else np.nanmean(MAG, axis=0)
            map_m = np.nanmean(MAG, axis=0)
            spring_neap_maps = (map_s, map_n, map_m)
            log(f"  ✓ Spring P{C['spring_neap_percentile_spring']} = {sn['threshold_spring']:.4f} m")
            log(f"  ✓ Neap P{C['spring_neap_percentile_neap']} = {sn['threshold_neap']:.4f} m")
            log(f"  ✓ Timestep spring/neap digunakan = {len(idx_s)} / {len(idx_n)}")
        except Exception as exc:
            log(f"  ⚠ Analisis Spring/Neap gagal/dilewati: {exc}")

    # -------------------------------------------------------------------------
    # 11B. Analisis Harmonik Pasut / Formzahl
    # -------------------------------------------------------------------------
    formzahl_result = None

    if C["analisis_formzahl_aktif"]:
        section("BAGIAN 11B: Analisis Harmonik Pasut / Formzahl")

        try:
            # Formzahl memakai S1 yang sama dengan Spring/Neap
            if "spring_neap" not in s1_series:
                raise ValueError(
                    "S1 referensi Spring/Neap tidak tersedia. "
                    "Pastikan analisis_spring_neap_aktif=True."
                )

            if "spring_neap" not in requests:
                raise ValueError(
                    "Titik referensi Spring/Neap tidak ditemukan."
                )

            formzahl_M, formzahl_N = requests["spring_neap"]

            # Ambil latitude dari grid GRD
            latitude_formzahl = 0.0

            if Y is not None:
                n0 = int(formzahl_N) - 1
                m0 = int(formzahl_M) - 1

                if (
                    0 <= n0 < Y.shape[0]
                    and 0 <= m0 < Y.shape[1]
                    and np.isfinite(Y[n0, m0])
                ):
                    latitude_formzahl = float(Y[n0, m0])
                else:
                    log(
                        "  ⚠ Latitude titik Formzahl tidak valid; "
                        "menggunakan latitude 0.0°."
                    )

            formzahl_result = hitung_formzahl_utide(
                eta=s1_series["spring_neap"],
                waktu=WAKTU[:len(s1_series["spring_neap"])],
                latitude_deg=latitude_formzahl,
            )

            # Simpan CSV
            fp_formzahl = (
                SUBDIRS["csv"]
                / f"formzahl_{safe_name(C['kode_selat'])}.csv"
            )

            pd.DataFrame([formzahl_result]).to_csv(
                fp_formzahl,
                index=False,
                encoding="utf-8-sig",
            )

            log(
                f"  ✓ Formzahl F = {formzahl_result['Formzahl_F']:.4f}"
            )
            log(
                f"  ✓ Tipe pasut = {formzahl_result['Tipe_Pasut']}"
            )
            log(
                f"  ✓ M2={formzahl_result['A_M2_m']:.4f} m | "
                f"S2={formzahl_result['A_S2_m']:.4f} m | "
                f"K1={formzahl_result['A_K1_m']:.4f} m | "
                f"O1={formzahl_result['A_O1_m']:.4f} m"
            )
            log(f"  ✓ CSV Formzahl: {fp_formzahl.name}")

        except Exception as exc:
            log(f"  ⚠ Analisis Formzahl gagal/dilewati: {exc}")
    # -------------------------------------------------------------------------
    # 12. Outputs
    # -------------------------------------------------------------------------
    save_maps_and_graphs(
        WAKTU, MAG, U_mean, V_mean, stats, POWER_mean, POWER_p95,
        valid_mask, layak_mask, optimal_mask,
        X, Y, DEPTH, enc_coords, bnd_segs,
        top_current, top_energy,
        quickplot_validation=quickplot_validation,
        rmse_validation=rmse_validation,
        spring_neap_result=spring_neap_result,
        spring_neap_maps=spring_neap_maps,
    )

    create_gif_daily(WAKTU, MAG, X, Y, enc_coords, bnd_segs)

    export_outputs(
        domain_summary, grid_rows, top_current, top_energy, spike_rows,
        quickplot_validation=quickplot_validation,
        rmse_validation=rmse_validation,
        spring_neap_result=spring_neap_result,
    )

    write_summary_txt(
    domain_summary, top_current, top_energy, depth_note,
    quickplot_validation=quickplot_validation,
    rmse_validation=rmse_validation,
    spring_neap_result=spring_neap_result,
    formzahl_result=formzahl_result,
)
    

    # Simpan log terminal
    log_path = OUTPUT_DIR / f"LOG_ANALISIS_FINAL_{safe_name(C['kode_selat'])}.txt"
    log_path.write_text("\n".join(LOG_LINES), encoding="utf-8")
    log("")
    log("=" * 78)
    log("  ANALISIS FINAL SELESAI")
    log("=" * 78)
    log(f"  Durasi : {(datetime.now() - t_global).total_seconds():.1f} detik")
    log(f"  Output : {OUTPUT_DIR}")
    log(f"  Log    : {log_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log("")
        log("!" * 78)
        log(f"ERROR FATAL: {exc}")
        log(traceback.format_exc())
        log("!" * 78)
        error_path = OUTPUT_DIR / f"ERROR_ANALISIS_FINAL_{safe_name(C['kode_selat'])}.txt"
        error_path.write_text("\n".join(LOG_LINES), encoding="utf-8")
        raise

# © Aimar Rendra P.A. Teknik Kelautan — Universitas Hasanuddin
