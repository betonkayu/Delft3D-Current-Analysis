r"""
© Aimar Rendra P.A. Teknik Kelautan — Universitas Hasanuddin
===============================================================================
ONE-RUN DATA PIPELINE DELFT3D-FLOW CLASSIC -> CSV GEOREFERENSI SIAP QGIS
STUDI KASUS: LEMBEH / CAPALULU / ALAS (pilih lewat variabel SELAT)
===============================================================================

Fungsi script
-------------
Script ini menggabungkan dua tahap:

1) Membaca output Delft3D-FLOW klasik NEFIS:
   trim-*.dat + trim-*.def -> statistik arus.

2) Membaca koordinat asli dari file grid Delft3D (*.grd), lalu langsung
   menyimpan CSV akhir dengan kolom x/y berupa koordinat geografis sebenarnya.

Dengan script ini Anda tidak perlu menjalankan terpisah:
- delft3d_to_qgis_NEFIS_fixed.py
- capalulu_final_qgis_pipeline.py

Output utama mengikuti path `output_final` pada CONFIG.

Cara menjalankan dari PowerShell
--------------------------------
python delft3d_nefis_to_csv.py

Catatan penting
---------------
- Script ini dijalankan di Python Windows / PowerShell, bukan di QGIS.
- Script membutuhkan numpy dan pandas.
- Script membutuhkan nefis.dll dari instalasi Delft3D 4.
- Script ini hanya untuk model barotropik KMAX=1.
===============================================================================
"""

from __future__ import annotations

import ctypes
import math
import os
import re
import sys
from ctypes import byref, c_char, c_int, create_string_buffer
from pathlib import Path

import numpy as np
import pandas as pd


# =============================================================================
# USER CONFIGURATION — EDIT PATH DI BAGIAN INI SAJA
# =============================================================================
# 1. NEFIS_DLL / NEFIS_DIR: lokasi instalasi Delft3D di komputer Anda.
# 2. CONFIG["lembeh"], CONFIG["capalulu"], CONFIG["alas"]:
#    isi path file input dan folder output sesuai model Anda.
#
# Contoh struktur folder yang mudah dipahami:
#   D:/PROJECT/alas/trim-alas.dat
#   D:/PROJECT/alas/trim-alas.def
#   D:/PROJECT/alas/alass4.grd
#   D:/PROJECT/alas/output
#
# Jangan mengubah bagian algoritma di bawah konfigurasi jika hanya ingin
# menjalankan analisis dengan data/model Anda sendiri.
# =============================================================================

# Lokasi NEFIS Delft3D 4. Sesuaikan jika instalasi Delft3D Anda berbeda.
NEFIS_DLL = r"C:\Program Files\Deltares\Delft3D 4.05.01\x64\share\bin\nefis.dll"
NEFIS_DIR = r"C:\Program Files\Deltares\Delft3D 4.05.01\x64\share\bin"

# Pilih selat yang diproses.
SELAT = "alas"  # pilihan: "lembeh", "capalulu", "alas"
SELAT = SELAT.strip().lower()  # agar "ALAS" / "Alas" tetap terbaca

CONFIG = {
    "lembeh": {
        "trim_dat": r"",  # <-- ISI path file trim.dat
        "trim_def": r"",  # <-- ISI path file trim.def
        "grid_file": r"",  # <-- ISI path file .grd
        "output_final": r"",  # <-- ISI folder output
        "tanggal_mulai": "2026-03-01 00:00:00",  # <-- sesuaikan dengan tanggal mulai simulasi
        "interval_menit": 10,
        "spinup_hari": 1,
        "spring_timestep": 3002,  # zero-based; sesuaikan dengan timestep tanggal spring yang digunakan
        "neap_timestep": 1751,    # zero-based
        "cartesian_crs": "",
    },
    "capalulu": {
        "trim_dat": r"",  # <-- ISI path file trim-capalulu.dat
        "trim_def": r"",  # <-- ISI path file trim-capalulu.def
        "grid_file": r"",  # <-- ISI path file grid .grd
        "output_final": r"",  # <-- ISI folder output
        "tanggal_mulai": "2026-03-01 00:00:00",
        "interval_menit": 60,
        "spinup_hari": 1,
        "spring_timestep": 538,  # zero-based
        "neap_timestep": 313,    # zero-based
        "cartesian_crs": "",
    },
    "alas": {
        "trim_dat": r"",  # <-- ISI path file trim-alas.dat
        "trim_def": r"",  # <-- ISI path file trim-alas.def
        "grid_file": r"",  # <-- ISI path file grid .grd
        "output_final": r"",  # <-- ISI folder output
        "tanggal_mulai": "2026-03-01 00:00:00",
        "interval_menit": 10,
        "spinup_hari": 1,
        "spring_timestep": 3240,  # zero-based
        "neap_timestep": 1450,    # zero-based
        "cartesian_crs": "",
    },
}

RHO = 1025.0
NAME_BUF = 17
MISSING_ABS_THRESHOLD_VELOCITY = 900.0
ONLY_ACTIVE_WATER_CELLS = True
DEFAULT_GRID_MISSING_VALUE = -999.0
FLOAT_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?")


# =============================================================================
# UTILITAS CETAK DAN VALIDASI
# =============================================================================

WATERMARK_TEXT = "© Aimar Rendra P.A. | Teknik Kelautan — Universitas Hasanuddin"


def section(title: str) -> None:
    print("\n" + "=" * 78)
    print(f"  {title}")
    print("=" * 78)


def info(message: str) -> None:
    print(f"  [OK] {message}")


def warning(message: str) -> None:
    print(f"  [PERINGATAN] {message}")


def fail(message: str) -> None:
    raise RuntimeError(message)


def require_file(path: Path, label: str) -> None:
    if not path.exists():
        fail(f"{label} tidak ditemukan: {path}")
    info(f"{label}: {path}")


def parse_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def numbers(text: str) -> list[float]:
    return [parse_float(token) for token in FLOAT_RE.findall(text)]


def timestamp_for_timestep(cfg: dict, t0: int) -> pd.Timestamp:
    return pd.Timestamp(cfg["tanggal_mulai"]) + pd.to_timedelta(
        t0 * cfg["interval_menit"], unit="m"
    )


def safe_divide(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    out = np.full_like(num, np.nan, dtype=np.float64)
    np.divide(num, den, out=out, where=den > 0)
    return out


def clean_velocity_values(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr, dtype=np.float32)
    arr = np.where(np.abs(arr) >= MISSING_ABS_THRESHOLD_VELOCITY, np.nan, arr)
    arr = np.where(np.abs(arr) > 1.0e10, np.nan, arr)
    return arr.astype(np.float32, copy=False)


def is_valid_coordinate(value: float, missing_value: float) -> bool:
    if not math.isfinite(value):
        return False
    if abs(value) > 1.0e20:
        return False
    tolerance = max(1.0e-8, abs(missing_value) * 1.0e-6)
    if math.isclose(value, missing_value, abs_tol=tolerance):
        return False
    return True


# =============================================================================
# PEMBACA FILE GRID DELFT3D RGFGRID (*.grd)
# =============================================================================

class GridData:
    def __init__(self, path: Path, coordinate_system: str, missing_value: float,
                 mmax: int, nmax: int, x: list[list[float]], y: list[list[float]]):
        self.path = path
        self.coordinate_system = coordinate_system
        self.missing_value = missing_value
        self.mmax = mmax
        self.nmax = nmax
        self.x = x
        self.y = y


def find_dimension_line(lines: list[str]) -> tuple[int, int, int]:
    for index, raw in enumerate(lines):
        stripped = raw.strip()
        if not stripped or stripped.startswith("*") or "=" in stripped:
            continue
        tokens = stripped.split()
        if len(tokens) not in (2, 3):
            continue
        try:
            ints = [int(token) for token in tokens]
        except ValueError:
            continue
        if ints[0] > 1 and ints[1] > 1:
            return index, ints[0], ints[1]
    fail("Baris dimensi MMAX NMAX tidak ditemukan dalam file .grd.")


def collect_eta_rows(lines: list[str], start_index: int) -> list[tuple[int, list[float]]]:
    rows: list[tuple[int, list[float]]] = []
    current_eta = None
    current_values: list[float] = []

    def flush() -> None:
        nonlocal current_eta, current_values
        if current_eta is not None:
            rows.append((current_eta, current_values))
        current_eta = None
        current_values = []

    for raw in lines[start_index:]:
        stripped = raw.strip()
        if not stripped or stripped.startswith("*"):
            continue

        match = re.match(r"^ETA\s*=\s*([+-]?\d+)\s*(.*)$", stripped, flags=re.I)
        if match:
            flush()
            current_eta = int(match.group(1))
            current_values = numbers(match.group(2))
            continue

        if current_eta is not None:
            current_values.extend(numbers(stripped))

    flush()
    return rows


def read_grd(path: Path) -> GridData:
    require_file(path, "GRID_FILE (.grd)")
    text = path.read_text(encoding="utf-8", errors="ignore")
    lines = text.splitlines()

    coord_match = re.search(r"Coordinate\s+System\s*=\s*([^\r\n]+)", text, flags=re.I)
    coordinate_system = coord_match.group(1).strip() if coord_match else "Cartesian (default)"

    miss_match = re.search(r"Missing\s+Value\s*=\s*([^\s\r\n]+)", text, flags=re.I)
    try:
        missing_value = parse_float(miss_match.group(1)) if miss_match else DEFAULT_GRID_MISSING_VALUE
    except ValueError:
        missing_value = DEFAULT_GRID_MISSING_VALUE

    dimension_line, mmax, nmax = find_dimension_line(lines)
    eta_rows = collect_eta_rows(lines, dimension_line + 1)

    if len(eta_rows) < 2 * nmax:
        fail(
            "File .grd belum dapat dibaca sebagai format RGFGRID. "
            f"Ditemukan {len(eta_rows)} blok ETA, minimal {2 * nmax} dibutuhkan."
        )

    x_rows = eta_rows[:nmax]
    y_rows = eta_rows[nmax: 2 * nmax]

    def build(rows: list[tuple[int, list[float]]], label: str) -> list[list[float]]:
        output = [[math.nan for _ in range(mmax)] for _ in range(nmax)]
        for fallback_row, (eta_index, vals) in enumerate(rows, start=1):
            if len(vals) < mmax:
                fail(f"Baris {label} ETA={eta_index} hanya memiliki {len(vals)} nilai; dibutuhkan {mmax}.")
            row_number = eta_index if 1 <= eta_index <= nmax else fallback_row
            output[row_number - 1] = vals[:mmax]
        return output

    x = build(x_rows, "X")
    y = build(y_rows, "Y")

    valid_pairs = []
    for n_index in range(nmax):
        for m_index in range(mmax):
            xv = x[n_index][m_index]
            yv = y[n_index][m_index]
            if is_valid_coordinate(xv, missing_value) and is_valid_coordinate(yv, missing_value):
                valid_pairs.append((xv, yv))

    if not valid_pairs:
        fail("Tidak ditemukan koordinat valid pada file .grd.")

    xs = [p[0] for p in valid_pairs]
    ys = [p[1] for p in valid_pairs]
    info(f"Coordinate System: {coordinate_system}")
    info(f"MMAX, NMAX dari .grd: {mmax}, {nmax}")
    info(f"Titik koordinat valid: {len(valid_pairs)}")
    info(f"Rentang X: {min(xs):.10f} s.d. {max(xs):.10f}")
    info(f"Rentang Y: {min(ys):.10f} s.d. {max(ys):.10f}")

    return GridData(path, coordinate_system, missing_value, mmax, nmax, x, y)


def determine_crs(grid: GridData, cartesian_crs: str) -> str:
    if "spher" in grid.coordinate_system.lower():
        return "EPSG:4326"
    if cartesian_crs.strip():
        return cartesian_crs.strip()
    return "BELUM DIISI — grid Cartesian membutuhkan cartesian_crs"


# =============================================================================
# NEFIS DLL DAN PEMBACAAN TRIM
# =============================================================================

def load_nefis():
    require_file(Path(NEFIS_DLL), "NEFIS DLL")
    os.environ["PATH"] = NEFIS_DIR + os.pathsep + os.environ.get("PATH", "")
    if hasattr(os, "add_dll_directory"):
        try:
            os.add_dll_directory(NEFIS_DIR)
        except Exception:
            pass

    try:
        nef = ctypes.cdll.LoadLibrary(NEFIS_DLL)
    except Exception as exc:
        raise RuntimeError(f"Gagal memuat NEFIS DLL: {exc}") from exc

    nef.Crenef.restype = c_int
    nef.Clsnef.restype = c_int
    nef.Getelt.restype = c_int
    return nef


def nefis_open(nef, dat_path: str, def_path: str):
    require_file(Path(dat_path), "Trim DAT")
    require_file(Path(def_path), "Trim DEF")

    fd = c_int(0)
    ret = nef.Crenef(
        byref(fd),
        dat_path.encode(),
        def_path.encode(),
        c_char(b"N"),
        c_char(b"r"),
    )
    if ret != 0:
        fail(f"Gagal membuka pasangan NEFIS, kode={ret}. DAT={dat_path} DEF={def_path}")
    return fd


def nefis_close(nef, fd) -> None:
    if fd is None:
        return
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
        fail(f"Getelt gagal untuk {group}/{element}, kode={ret}")
    return int(buf[0])


def read_const_int(nef, fd, element: str, count: int) -> np.ndarray:
    g = create_string_buffer(b"map-const", NAME_BUF)
    e = create_string_buffer(element.encode(), NAME_BUF)
    ui = (c_int * 3)(1, 1, 1)
    uo = (c_int * 1)(1)
    bl = c_int(count * 4)
    buf = (ctypes.c_int * count)()
    ret = nef.Getelt(byref(fd), g, e, ui, uo, byref(bl), buf)
    if ret != 0:
        fail(f"Getelt gagal untuk map-const/{element}, kode={ret}")
    return np.asarray(buf[:count], dtype=np.int32)


def read_series_time_only(nef, fd, element: str, t0: int, count: int) -> np.ndarray:
    g = create_string_buffer(b"map-series", NAME_BUF)
    e = create_string_buffer(element.encode(), NAME_BUF)
    ui = (c_int * 3)(t0 + 1, t0 + 1, 1)
    uo = (c_int * 1)(1)
    bl = c_int(count * 4)
    buf = (ctypes.c_float * count)()
    ret = nef.Getelt(byref(fd), g, e, ui, uo, byref(bl), buf)
    if ret != 0:
        raise EOFError(f"Getelt map-series/{element} berhenti pada timestep {t0 + 1}, kode={ret}")
    return np.asarray(buf[:count], dtype=np.float32)


def reshape_f(flat: np.ndarray, nmax: int, mmax: int, kmax: int = 1) -> np.ndarray:
    return flat.reshape((nmax, mmax, kmax), order="F")[:, :, 0]


def reshape_const_f(flat: np.ndarray, nmax: int, mmax: int) -> np.ndarray:
    return flat.reshape((nmax, mmax), order="F")


def reconstruct_backward_velocity(U: np.ndarray, V: np.ndarray):
    Uc = np.full_like(U, np.nan, dtype=np.float32)
    Vc = np.full_like(V, np.nan, dtype=np.float32)
    Uc[:, 1:] = 0.5 * (U[:, :-1] + U[:, 1:])
    Vc[1:, :] = 0.5 * (V[:-1, :] + V[1:, :])
    speed = np.sqrt(Uc**2 + Vc**2).astype(np.float32)
    return Uc, Vc, speed


def read_kcs_mask(nef, fd, nmax: int, mmax: int) -> np.ndarray:
    count = nmax * mmax
    try:
        kcs = reshape_const_f(read_const_int(nef, fd, "KCS", count), nmax, mmax)
        if ONLY_ACTIVE_WATER_CELLS:
            mask = kcs == 1
            info(f"Mask KCS: hanya sel aktif KCS == 1 ({int(mask.sum())} sel)")
        else:
            mask = kcs > 0
            info(f"Mask KCS: seluruh sel air KCS > 0 ({int(mask.sum())} sel)")
        return mask
    except Exception as exc:
        warning(f"KCS tidak terbaca ({exc}); mask akan ditentukan dari data valid.")
        return np.ones((nmax, mmax), dtype=bool)


def count_timesteps(nef, fd, count_uv: int) -> int:
    info("Menghitung jumlah timestep dari map-series/U1 ...")
    n = 0
    while True:
        try:
            _ = read_series_time_only(nef, fd, "U1", n, count_uv)
            n += 1
            if n % 500 == 0:
                print(f"    {n} timestep terdeteksi ...")
        except EOFError:
            break
    if n <= 0:
        fail("Tidak ada timestep map-series/U1 yang dapat dibaca.")
    info(f"Jumlah timestep: {n}")
    return n


# =============================================================================
# EKSPOR CSV FINAL GEOREFERENSI
# =============================================================================

def grid_arrays_from_grd(grid: GridData) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray(grid.x, dtype=np.float64)
    y = np.asarray(grid.y, dtype=np.float64)
    valid = np.zeros((grid.nmax, grid.mmax), dtype=bool)
    for n in range(grid.nmax):
        for m in range(grid.mmax):
            valid[n, m] = is_valid_coordinate(float(x[n, m]), grid.missing_value) and is_valid_coordinate(float(y[n, m]), grid.missing_value)
    return x, y, valid




def align_grid_arrays_to_trim(
    grid: GridData,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    grid_coord_valid: np.ndarray,
    trim_mmax: int,
    trim_nmax: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Menyamakan ukuran koordinat .grd dengan dimensi TRIM.

    Pada beberapa output Delft3D-FLOW klasik, MMAX/NMAX pada TRIM dapat lebih
    besar satu baris/kolom dibanding file .grd karena baris/kolom dummy/boundary
    pada susunan NEFIS/staggered grid. Kondisi ini aman ditangani dengan padding
    NaN pada sisi luar sehingga sel dummy tidak ikut diekspor.
    """
    if grid.mmax == trim_mmax and grid.nmax == trim_nmax:
        info("Dimensi .grd sama dengan trim.")
        return x_grid, y_grid, grid_coord_valid

    if grid.mmax == trim_mmax - 1 and grid.nmax == trim_nmax - 1:
        warning(
            f"Dimensi .grd ({grid.mmax}, {grid.nmax}) lebih kecil 1 dari trim "
            f"({trim_mmax}, {trim_nmax}). Ini umum pada grid Delft3D tertentu "
            "karena baris/kolom dummy/boundary di TRIM. Script akan menambahkan "
            "padding NaN pada sisi luar dan tidak mengekspor sel dummy."
        )
        x_pad = np.full((trim_nmax, trim_mmax), np.nan, dtype=np.float64)
        y_pad = np.full((trim_nmax, trim_mmax), np.nan, dtype=np.float64)
        valid_pad = np.zeros((trim_nmax, trim_mmax), dtype=bool)
        x_pad[:grid.nmax, :grid.mmax] = x_grid
        y_pad[:grid.nmax, :grid.mmax] = y_grid
        valid_pad[:grid.nmax, :grid.mmax] = grid_coord_valid
        return x_pad, y_pad, valid_pad

    fail(
        f"Dimensi .grd ({grid.mmax}, {grid.nmax}) tidak cocok dengan trim "
        f"({trim_mmax}, {trim_nmax}). Jika selisihnya bukan +1, kemungkinan "
        "GRID_FILE bukan grid yang dipakai oleh trim tersebut."
    )


def build_dataframe(
    x: np.ndarray,
    y: np.ndarray,
    valid_mask: np.ndarray,
    speed_mean: np.ndarray,
    speed_max: np.ndarray,
    speed_p90: np.ndarray,
    speed_spring: np.ndarray,
    speed_neap: np.ndarray,
    power_mean: np.ndarray,
    u_mean: np.ndarray,
    v_mean: np.ndarray,
    dir_mean: np.ndarray,
    u_spring: np.ndarray,
    v_spring: np.ndarray,
    u_neap: np.ndarray,
    v_neap: np.ndarray,
    prob_v05: np.ndarray,
    prob_v075: np.ndarray,
    prob_v10: np.ndarray,
) -> pd.DataFrame:
    n0, m0 = np.where(valid_mask)
    return pd.DataFrame({
        "M": m0 + 1,
        "N": n0 + 1,
        "x": x[n0, m0],
        "y": y[n0, m0],
        "speed_mean": speed_mean[n0, m0],
        "speed_max": speed_max[n0, m0],
        "speed_p90": speed_p90[n0, m0],
        "speed_spring": speed_spring[n0, m0],
        "speed_neap": speed_neap[n0, m0],
        "power_mean": power_mean[n0, m0],
        "u_mean": u_mean[n0, m0],
        "v_mean": v_mean[n0, m0],
        "dir_mean_deg": dir_mean[n0, m0],
        "u_spring": u_spring[n0, m0],
        "v_spring": v_spring[n0, m0],
        "u_neap": u_neap[n0, m0],
        "v_neap": v_neap[n0, m0],
        "prob_v05_pct": prob_v05[n0, m0],
        "prob_v075_pct": prob_v075[n0, m0],
        "prob_v10_pct": prob_v10[n0, m0],
    })


def export_csv(df: pd.DataFrame, output_dir: Path, selat: str) -> pd.Series:
    output_dir.mkdir(parents=True, exist_ok=True)

    all_path = output_dir / f"{selat}_all_stats.csv"
    df.to_csv(all_path, index=False, float_format="%.6f", encoding="utf-8-sig")
    info(f"Tersimpan: {all_path}")

    for col in [
        "speed_mean",
        "speed_spring",
        "speed_neap",
        "power_mean",
        "speed_p90",
        "prob_v075_pct",
        "prob_v10_pct",
    ]:
        out = output_dir / f"{selat}_{col}.csv"
        df[["M", "N", "x", "y", col]].rename(columns={col: "value"}).to_csv(
            out, index=False, float_format="%.6f", encoding="utf-8-sig"
        )
        info(f"Tersimpan: {out}")

    vector_path = output_dir / f"{selat}_vector_stats.csv"
    df[[
        "M", "N", "x", "y", "u_mean", "v_mean", "dir_mean_deg",
        "u_spring", "v_spring", "u_neap", "v_neap",
    ]].to_csv(vector_path, index=False, float_format="%.6f", encoding="utf-8-sig")
    info(f"Tersimpan: {vector_path}")

    idx = df["speed_mean"].idxmax()
    dominant = df.loc[[idx]].copy()
    dominant_path = output_dir / f"{selat}_titik_dominan.csv"
    dominant.to_csv(dominant_path, index=False, float_format="%.6f", encoding="utf-8-sig")
    info(f"Tersimpan: {dominant_path}")

    row = dominant.iloc[0]
    section("TITIK DOMINAN BERDASARKAN SPEED MEAN")
    print(f"  M, N       : {int(row['M'])}, {int(row['N'])}")
    print(f"  X, Y       : {row['x']:.10f}, {row['y']:.10f}")
    print(f"  Speed mean : {row['speed_mean']:.4f} m/s")
    print(f"  Speed max  : {row['speed_max']:.4f} m/s")
    print(f"  Power mean : {row['power_mean']:.2f} W/m2")
    return row


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    selat_key = str(SELAT).strip().lower()
    if selat_key not in CONFIG:
        fail(f"SELAT='{SELAT}' tidak ada dalam CONFIG. Pilihan yang valid: {', '.join(CONFIG.keys())}")

    cfg = CONFIG[selat_key]
    for _key in ("trim_dat", "trim_def", "grid_file", "output_final"):
        if not str(cfg[_key]).strip():
            fail(f"CONFIG['{selat_key}']['{_key}'] masih kosong. Isi path di bagian USER CONFIGURATION.")
    globals()["SELAT"] = selat_key
    trim_dat = cfg["trim_dat"]
    trim_def = cfg["trim_def"]
    grid_file = Path(cfg["grid_file"])
    output_dir = Path(cfg["output_final"])
    temp_speed_path = output_dir / f"_{SELAT}_speed_postspinup.tmp"

    section(f"ONE-RUN DATA PIPELINE — SELAT {SELAT.upper()}")
    print(f"  {WATERMARK_TEXT}")
    print(f"  DAT        : {trim_dat}")
    print(f"  DEF        : {trim_def}")
    print(f"  GRID       : {grid_file}")
    print(f"  OUTPUT     : {output_dir}")

    section("MEMBACA KOORDINAT GRID .GRD")
    grid = read_grd(grid_file)
    x_grid, y_grid, grid_coord_valid = grid_arrays_from_grd(grid)
    crs = determine_crs(grid, cfg.get("cartesian_crs", ""))
    info(f"CRS untuk QGIS: {crs}")

    section("MEMBUKA FILE TRIM NEFIS")
    nef = load_nefis()
    fd = nefis_open(nef, trim_dat, trim_def)

    try:
        section("MEMBACA DIMENSI TRIM")
        mmax = read_int_scalar(nef, fd, "map-const", "MMAX")
        nmax = read_int_scalar(nef, fd, "map-const", "NMAX")
        try:
            kmax = read_int_scalar(nef, fd, "map-const", "KMAX")
        except Exception:
            kmax = 1

        info(f"MMAX, NMAX, KMAX dari trim: {mmax}, {nmax}, {kmax}")
        if kmax != 1:
            fail(
                f"KMAX={kmax}. Script ini hanya tervalidasi untuk model barotropik KMAX=1. "
                "Jangan lanjutkan sebelum metode vertical averaging ditentukan."
            )

        x_grid, y_grid, grid_coord_valid = align_grid_arrays_to_trim(
            grid, x_grid, y_grid, grid_coord_valid, mmax, nmax
        )

        kcs_mask = read_kcs_mask(nef, fd, nmax, mmax)
        nval = nmax * mmax * kmax

        section("MENGHITUNG TIMESTEP")
        n_time = count_timesteps(nef, fd, nval)
        spinup = int(round(cfg["spinup_hari"] * 24 * 60 / cfg["interval_menit"]))
        n_analysis = n_time - spinup
        if n_analysis <= 0:
            fail(f"Spin-up={spinup} timestep tidak boleh >= total timestep={n_time}.")

        spring_t = int(cfg["spring_timestep"])
        neap_t = int(cfg["neap_timestep"])
        for label, t0 in [("spring", spring_t), ("neap", neap_t)]:
            if not (0 <= t0 < n_time):
                fail(f"{label}_timestep={t0} di luar rentang zero-based 0 sampai {n_time - 1}.")

        info(f"Spin-up dibuang : {spinup} timestep")
        info(f"Data dianalisis : {n_analysis} timestep")
        info(f"Spring timestep : {spring_t} = {timestamp_for_timestep(cfg, spring_t)}")
        info(f"Neap timestep   : {neap_t} = {timestamp_for_timestep(cfg, neap_t)}")

        section("MEMBACA U1/V1 DAN MENGHITUNG STATISTIK")
        output_dir.mkdir(parents=True, exist_ok=True)
        shape = (nmax, mmax)
        base_valid = kcs_mask & grid_coord_valid

        sum_speed = np.zeros(shape, dtype=np.float64)
        count_speed = np.zeros(shape, dtype=np.int64)
        speed_max = np.full(shape, np.nan, dtype=np.float32)
        sum_power = np.zeros(shape, dtype=np.float64)
        count_power = np.zeros(shape, dtype=np.int64)
        sum_u = np.zeros(shape, dtype=np.float64)
        sum_v = np.zeros(shape, dtype=np.float64)
        count_uv = np.zeros(shape, dtype=np.int64)
        count_v05 = np.zeros(shape, dtype=np.int64)
        count_v075 = np.zeros(shape, dtype=np.int64)
        count_v10 = np.zeros(shape, dtype=np.int64)

        speed_spring = np.full(shape, np.nan, dtype=np.float32)
        speed_neap = np.full(shape, np.nan, dtype=np.float32)
        u_spring = np.full(shape, np.nan, dtype=np.float32)
        v_spring = np.full(shape, np.nan, dtype=np.float32)
        u_neap = np.full(shape, np.nan, dtype=np.float32)
        v_neap = np.full(shape, np.nan, dtype=np.float32)

        # Disk-backed array agar P90 tidak menghabiskan RAM.
        speed_store = np.memmap(
            temp_speed_path,
            dtype="float32",
            mode="w+",
            shape=(n_analysis, nmax, mmax),
        )
        speed_store[:] = np.nan

        analysis_i = 0
        for t0 in range(n_time):
            u_flat = read_series_time_only(nef, fd, "U1", t0, nval)
            v_flat = read_series_time_only(nef, fd, "V1", t0, nval)
            U = clean_velocity_values(reshape_f(u_flat, nmax, mmax, kmax))
            V = clean_velocity_values(reshape_f(v_flat, nmax, mmax, kmax))
            Uc, Vc, speed = reconstruct_backward_velocity(U, V)

            Uc = np.where(base_valid, Uc, np.nan).astype(np.float32)
            Vc = np.where(base_valid, Vc, np.nan).astype(np.float32)
            speed = np.where(base_valid, speed, np.nan).astype(np.float32)

            if t0 == spring_t:
                speed_spring = speed.copy()
                u_spring = Uc.copy()
                v_spring = Vc.copy()
            if t0 == neap_t:
                speed_neap = speed.copy()
                u_neap = Uc.copy()
                v_neap = Vc.copy()

            if t0 < spinup:
                continue

            speed_store[analysis_i] = speed
            finite_speed = np.isfinite(speed)
            finite_uv = np.isfinite(Uc) & np.isfinite(Vc)
            power = 0.5 * RHO * speed.astype(np.float64) ** 3
            finite_power = np.isfinite(power)

            sum_speed[finite_speed] += speed[finite_speed]
            count_speed[finite_speed] += 1
            if analysis_i == 0:
                speed_max = speed.copy()
            else:
                speed_max = np.fmax(speed_max, speed)

            sum_power[finite_power] += power[finite_power]
            count_power[finite_power] += 1
            sum_u[finite_uv] += Uc[finite_uv]
            sum_v[finite_uv] += Vc[finite_uv]
            count_uv[finite_uv] += 1
            count_v05[finite_speed] += speed[finite_speed] >= 0.50
            count_v075[finite_speed] += speed[finite_speed] >= 0.75
            count_v10[finite_speed] += speed[finite_speed] >= 1.00

            analysis_i += 1
            if analysis_i % 250 == 0 or analysis_i == n_analysis:
                print(f"    {analysis_i}/{n_analysis} timestep post-spin-up selesai")

        speed_store.flush()

        speed_mean = safe_divide(sum_speed, count_speed)
        power_mean = safe_divide(sum_power, count_power)
        u_mean = safe_divide(sum_u, count_uv)
        v_mean = safe_divide(sum_v, count_uv)
        dir_mean = (np.degrees(np.arctan2(u_mean, v_mean)) + 360.0) % 360.0
        prob_v05 = safe_divide(count_v05.astype(float), count_speed) * 100.0
        prob_v075 = safe_divide(count_v075.astype(float), count_speed) * 100.0
        prob_v10 = safe_divide(count_v10.astype(float), count_speed) * 100.0

        info("Menghitung P90 dari array disk-backed ...")
        with np.errstate(all="ignore"):
            speed_p90 = np.nanpercentile(speed_store, 90, axis=0)

        valid_mask = base_valid & np.isfinite(speed_mean)
        if not np.any(valid_mask):
            fail("Tidak ada sel grid valid untuk diekspor.")

        section("MENYUSUN DAN MENYIMPAN CSV FINAL GEOREFERENSI")
        df = build_dataframe(
            x=x_grid,
            y=y_grid,
            valid_mask=valid_mask,
            speed_mean=speed_mean,
            speed_max=speed_max,
            speed_p90=speed_p90,
            speed_spring=speed_spring,
            speed_neap=speed_neap,
            power_mean=power_mean,
            u_mean=u_mean,
            v_mean=v_mean,
            dir_mean=dir_mean,
            u_spring=u_spring,
            v_spring=v_spring,
            u_neap=u_neap,
            v_neap=v_neap,
            prob_v05=prob_v05,
            prob_v075=prob_v075,
            prob_v10=prob_v10,
        )
        dominant_row = export_csv(df, output_dir, SELAT)

        summary_path = output_dir / f"{SELAT}_RINGKASAN_ONE_RUN_GET_DATA.txt"
        summary_path.write_text(
            "\n".join([
                WATERMARK_TEXT,
                f"ONE-RUN DATA PIPELINE DELFT3D-FLOW CLASSIC -> CSV GEOREFERENSI — {SELAT.upper()}",
                "=" * 78,
                f"DAT                 : {trim_dat}",
                f"DEF                 : {trim_def}",
                f"GRID                : {grid_file}",
                f"Coordinate System   : {grid.coordinate_system}",
                f"CRS QGIS            : {crs}",
                f"MMAX, NMAX, KMAX    : {mmax}, {nmax}, {kmax}",
                f"Jumlah timestep     : {n_time}",
                f"Spin-up dibuang     : {spinup} timestep",
                f"Timestep dianalisis : {n_analysis}",
                f"Spring zero-based   : {spring_t} | {timestamp_for_timestep(cfg, spring_t)}",
                f"Neap zero-based     : {neap_t} | {timestamp_for_timestep(cfg, neap_t)}",
                f"Jumlah sel diekspor : {len(df)}",
                f"Titik dominan M/N   : {int(dominant_row['M'])}, {int(dominant_row['N'])}",
                f"Titik dominan x/y   : {dominant_row['x']:.10f}, {dominant_row['y']:.10f}",
                f"Speed mean dominan  : {dominant_row['speed_mean']:.4f} m/s",
                f"Power mean dominan  : {dominant_row['power_mean']:.2f} W/m2",
                "",
                "Output CSV ini sudah memakai koordinat asli dari file .grd.",
                "Untuk QGIS: X field = x, Y field = y.",
            ]),
            encoding="utf-8",
        )
        info(f"Ringkasan tersimpan: {summary_path}")

        section("SELESAI")
        print(f"  Semua CSV final georeferensi berada di: {output_dir}")
        print("  Lanjutkan ke script QGIS one-run final untuk membuat peta.")

    finally:
        nefis_close(nef, fd)
        try:
            if "speed_store" in locals():
                del speed_store
            if temp_speed_path.exists():
                temp_speed_path.unlink()
        except Exception:
            warning(f"File sementara belum terhapus: {temp_speed_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("\n" + "!" * 78)
        print(f"ERROR: {exc}")
        print("!" * 78)
        sys.exit(1)

# © Aimar Rendra P.A.
# Teknik Kelautan — Universitas Hasanuddin