# Delft3D Current Analysis

Python workflow for processing Delft3D-FLOW hydrodynamic model results, extracting current data through NEFIS, performing current and energy analysis, and generating QGIS map outputs.

## Repository Structure

```text
Delft3D-Current-Analysis/
├── README.md
├── delft3d_current_analysis.py
├── delft3d_nefis_to_csv.py
├── qgis_map_export.py
└── .gitignore
```

## Scripts

### `delft3d_current_analysis.py`
Main script for current and kinetic-energy analysis from Delft3D-FLOW output.

Includes current-speed statistics, current direction, dominant-current-point analysis, velocity-threshold probability, spring/neap analysis, kinetic power density, and harmonic analysis when UTide is available.

### `delft3d_nefis_to_csv.py`
Reads Delft3D-FLOW classic NEFIS output (`.dat` + `.def`) and Delft3D grid coordinates (`.grd`), checks/aligns dimensions, processes current data, and exports georeferenced CSV results.

### `qgis_map_export.py`
QGIS workflow for processing spatial CSV results, interpolation/raster products, masks, and final map/layout exports.

## Requirements

### Delft3D / NEFIS

The NEFIS-based scripts require Delft3D 4 and its NEFIS library.

Official Deltares download pages:

- https://download.deltares.nl/delft3d-4-suite
- https://oss.deltares.nl/web/delft3d/downloads

Example NEFIS DLL path:

```text
C:\Program Files\Deltares\Delft3D 4.05.01\x64\share\bin\nefis.dll
```

If Delft3D is installed elsewhere, change the NEFIS path in the script configuration. **Do not upload `nefis.dll` to GitHub.** Install Delft3D separately from Deltares.

### Python

Recommended: Python 3.9+ on Windows.

Install the third-party packages:

```bash
pip install numpy pandas matplotlib scipy utide imageio Pillow
```

Packages/modules used include:

```python
import numpy as np
import pandas as pd
from utide import solve
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, BoundaryNorm
import imageio.v2 as imageio
```

The following are Python standard-library modules and require no separate installation:

```python
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
```

### QGIS

`qgis_map_export.py` is intended to run in a compatible QGIS Python environment / QGIS Python Console.

## Input Data

Typical Delft3D files include:

```text
trim-*.dat
trim-*.def
*.grd
*.dep
*.enc
*.bnd
```

Additional CSV or observation files may be required depending on the configured workflow.

The repository contains the processing scripts, not the original research/model data.

## Configuration

Open the `CONFIG` section near the beginning of each script and replace the example local paths with paths on your own computer. For example:

```python
CONFIG = {
    "trim_dat": r"D:\PROJECT\model\trim.dat",
    "trim_def": r"D:\PROJECT\model\trim.def",
    "grid_file": r"D:\PROJECT\model\model.grd",
    "output_final": r"D:\PROJECT\results",
}
```

Do not change analysis logic unless you understand the effect of the parameter being modified.

## NEFIS Files

For Delft3D-FLOW classic output, the `.dat` and `.def` files belong to the same NEFIS dataset and must correspond to the same simulation output.

Example:

```text
trim-current.dat
trim-current.def
```

The corresponding Delft3D grid (`.grd`) is also required when georeferenced coordinates are needed.

## Installation

1. Install Delft3D 4 from Deltares and confirm that `nefis.dll` is available.
2. Install Python 3.9 or newer.
3. Install packages:

```bash
pip install numpy pandas matplotlib scipy utide imageio Pillow
```

4. Install QGIS if `qgis_map_export.py` is required.
5. Configure local input/output paths in the scripts.
6. Run the script required for the processing stage.

## Workflow

```text
Delft3D-FLOW Output
        │
        ├── NEFIS processing
        │       └── delft3d_nefis_to_csv.py
        │
        ├── Current & energy analysis
        │       └── delft3d_current_analysis.py
        │
        └── GIS processing
                └── qgis_map_export.py
```

These scripts can represent separate processing stages. They are not necessarily a strict dependency chain; use each script according to the input/output files configured for the workflow.

## Data Privacy

Do not upload private or restricted research data to GitHub unless intentionally approved for public release.

Typical files to keep outside the public repository include:

```text
*.dat
*.def
*.grd
*.dep
*.enc
*.bnd
*.obs
*.csv
*.xlsx
```

Also replace personal computer paths such as `C:\Users\...` and private research directories before publishing.

## Purpose

This repository demonstrates a Python-based workflow for Delft3D-FLOW post-processing, NEFIS data extraction, current velocity analysis, kinetic power-density analysis, georeferenced CSV generation, GIS processing, and QGIS map automation for marine current-energy assessment.

## Author

**Aimar Rendra P.A.**  
Teknik Kelautan — Universitas Hasanuddin
