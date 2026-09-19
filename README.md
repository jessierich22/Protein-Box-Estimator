# Protein Box Estimator

Analyze molecular dimensions, estimate simulation box sizes, and screen box adequacy throughout molecular dynamics trajectories.

## Overview

`estimate_boxes.py` uses MDAnalysis to evaluate either a single molecular structure or a complete trajectory.

Two analysis modes are available:

1. **Static structure analysis** using `-s` only;
2. **Trajectory analysis** using `-s` together with `-f`.

The script can calculate the maximum atom-to-atom distance (`Dmax`), Cartesian molecular spans, recommended cubic and orthorhombic box dimensions, approximate solvent content, molecular extension throughout a trajectory, and Cartesian clearance relative to the recorded simulation box.

## Features

### Static structure mode

- Calculates the maximum pairwise atomic distance (`Dmax`);
- Identifies the two atoms defining `Dmax`;
- Calculates the molecular spans along the X, Y, and Z axes;
- Estimates conservative cubic and axis-based orthorhombic boxes;
- Estimates the number of SPC216 water molecules and the approximate total number of atoms.

### Trajectory mode

- Calculates `Dmax` for every analyzed frame;
- Measures changes in `Dmax` relative to the first analyzed frame;
- Identifies the frame and time of maximum molecular extension;
- Tracks the maximum molecular span along each Cartesian axis;
- Estimates the required cubic box edge throughout the trajectory;
- Compares molecular spans with the actual box dimensions;
- Applies a Cartesian clearance screening criterion;
- Writes the per-frame results to a CSV file.

## Requirements

- Python 3.11
- NumPy 1.26.4
- MDAnalysis 2.9.0

## Installation

Clone the repository:

```bash
git clone https://github.com/jessierich22/Protein-Box-Estimator.git
cd Protein-Box-Estimator
```

### Using Conda

```bash
conda env create -f environment.yml
conda activate protein-box-estimator
```

### Using pip

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Usage

### Analyze a static structure

The structure can be provided in any format supported by MDAnalysis, such as GRO or PDB:

```bash
python estimate_boxes.py \
    -s protein.gro \
    --selection "protein" \
    --cutoff 1.2 \
    --extra 1.0
```

Example using a PDB file:

```bash
python estimate_boxes.py -s protein.pdb --selection "protein"
```

In this mode, the results are printed to the terminal and no CSV file is generated.

### Analyze a trajectory

Provide a structure or topology file with `-s` and a trajectory with `-f`:

```bash
python estimate_boxes.py \
    -s system.tpr \
    -f trajectory.xtc \
    --selection "protein" \
    --cutoff 1.2 \
    --extra 1.0 \
    --stride 1 \
    --block-size 2000 \
    -o box_analysis.csv
```

The topology and trajectory files must be mutually compatible.

## Command-line options

| Argument | Default | Description |
|---|---:|---|
| `-s`, `--structure` | Required | Structure or topology file supported by MDAnalysis, such as PDB, GRO, or TPR |
| `-f`, `--trajectory` | None | Optional trajectory file, such as XTC or TRR |
| `--selection` | `protein` | MDAnalysis atom-selection expression |
| `--cutoff` | `1.2` nm | Long-range cutoff used in the box-size estimate |
| `--extra` | `1.0` nm | Additional safety distance |
| `--block-size` | `2000` | Number of atoms in each distance-calculation block |
| `--stride` | `1` | Analyze every Nth trajectory frame |
| `-o`, `--output` | `box_analysis.csv` | Output CSV file used in trajectory mode |

Display the command-line help with:

```bash
python estimate_boxes.py --help
```

## Methodology

### Maximum molecular dimension

The maximum molecular dimension, $D_{\max}$, is the largest distance between any two selected atoms:

$$
D_{\max} = \max_{i,j}\left\|\mathbf{r}_i-\mathbf{r}_j\right\|
$$

Pairwise distances are evaluated in blocks to reduce peak memory usage. The calculation remains quadratic with respect to the number of selected atoms.

### Cartesian molecular spans

For each Cartesian axis $k$, the molecular span is calculated as:

$$
S_k = r_{k,\max}-r_{k,\min}
$$

where $k$ corresponds to X, Y, or Z.

### Estimated box dimensions

The total size added to each molecular dimension is:

$$
d_{\mathrm{added}} = r_{\mathrm{cutoff}} + r_{\mathrm{extra}}
$$

The conservative cubic box edge is calculated from `Dmax`:

$$
L_{\mathrm{cubic}} = D_{\max}+d_{\mathrm{added}}
$$

The axis-based orthorhombic box dimensions are:

$$
L_k = S_k+d_{\mathrm{added}}
$$

These equations describe the total amount added to each full box dimension by the current implementation.

### Water and system-size estimates

Water density is derived from the GROMACS SPC216 reference box:

$$
\rho_{\mathrm{water}} = \frac{216}{(1.86206\ \mathrm{nm})^3}
$$

For a box with dimensions $L_x$, $L_y$, and $L_z$:

$$
V_{\mathrm{box}}=L_xL_yL_z
$$

$$
N_{\mathrm{water}}\approx V_{\mathrm{box}}\rho_{\mathrm{water}}
$$

The approximate total number of atoms is:

$$
N_{\mathrm{total}}=N_{\mathrm{solute}}+3N_{\mathrm{water}}
$$

The solvent estimate is an upper bound because the volume occupied by the solute is not subtracted.

## Trajectory extension analysis

The first analyzed frame is used as the reference. For each frame, the script reports:

$$
\Delta D_{\max}(t)=D_{\max}(t)-D_{\max}(t_0)
$$

and the relative change:

$$
\Delta D_{\max,\%}(t)=
\frac{D_{\max}(t)-D_{\max}(t_0)}{D_{\max}(t_0)}\times100
$$

When `--stride` is greater than 1, the first and maximum values refer only to the frames included in the analysis.

## Cartesian clearance screening

When trajectory box dimensions are available, the script calculates for each axis:

$$
C_k=L_{\mathrm{box},k}-S_k
$$

The minimum clearance observed across all axes and analyzed frames is classified as:

| Status | Criterion | Interpretation |
|---|---|---|
| `CRITICAL` | Clearance < 0 | The molecular span exceeded a box dimension |
| `FLAG` | 0 ≤ clearance < cutoff | The frame should be inspected with an explicit periodic-image distance analysis |
| `PASS` | Clearance ≥ cutoff | No frame was flagged by the Cartesian span criterion |

This clearance is a geometric screening metric. It is **not** the minimum atom-to-atom distance between the selected molecule and its periodic image.

## Trajectory CSV output

In trajectory mode, the CSV file contains:

| Column | Description |
|---|---|
| `frame` | Trajectory frame index |
| `time_ps` | Frame time in picoseconds |
| `time_ns` | Frame time in nanoseconds |
| `dmax_nm` | Maximum pairwise distance |
| `delta_dmax_nm` | Change in `Dmax` relative to the first analyzed frame |
| `delta_dmax_percent` | Relative change in `Dmax` |
| `required_cubic_box_nm` | Estimated required cubic box edge |
| `span_x_nm`, `span_y_nm`, `span_z_nm` | Cartesian molecular spans |
| `box_x_nm`, `box_y_nm`, `box_z_nm` | Box dimensions recorded in the trajectory |
| `clearance_x_nm`, `clearance_y_nm`, `clearance_z_nm` | Cartesian clearance values |

## Periodic boundary conditions

The selected molecule must be whole in the structure and in every analyzed trajectory frame. A molecule fragmented across periodic boundaries can produce artificial spans and incorrect `Dmax` values.

Prepare the trajectory with an appropriate periodic-boundary treatment before running the analysis. For GROMACS trajectories, this commonly involves making molecules whole and centering the selected solute with `gmx trjconv`. The exact procedure depends on the system and topology.

## Assumptions and limitations

- Only atoms matching the MDAnalysis selection are analyzed.
- Static analysis describes only the supplied conformation.
- Trajectory analysis describes only the frames selected by `--stride`.
- Pairwise-distance calculations can be computationally demanding for large selections and long trajectories.
- The `--block-size` option controls peak memory use but does not eliminate the quadratic computational cost.
- Water estimates assume a three-site SPC-like water model.
- Water estimates do not subtract the volume occupied by the solute.
- Additional ions required for neutralization or ionic strength are not predicted.
- Cartesian clearance is a screening metric and does not replace an explicit protein–periodic-image distance calculation.
- Triclinic box geometry and box angles are not considered by the Cartesian clearance screening.
- Results should be treated as planning and diagnostic estimates rather than exact system-construction values.

## Author

Jéssica Ierich

## Repository access

This repository is currently intended for internal research use by authorized collaborators.
