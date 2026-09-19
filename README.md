# Protein Box Estimator

Estimate conservative cubic simulation box dimensions and approximate system sizes for protein molecular dynamics simulations.

## Overview

`estimate_boxes.py` processes a database containing multiple protein structures and estimates:

- Maximum pairwise atomic distance (`Dmax`);
- Conservative cubic box edge length;
- Box volume;
- Approximate number of water molecules;
- Approximate total number of atoms in the solvated system.

The script calculates pairwise distances in blocks to reduce memory usage compared with allocating the complete distance matrix.

## Requirements

- Python 3.9 or later
- NumPy
- MDAnalysis

## Installation

Clone the repository and create a virtual environment:

```bash
git clone https://github.com/jessierich22/Protein-Box-Estimator.git
cd Protein-Box-Estimator

python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Alternatively, install the dependencies directly:

```bash
pip install numpy MDAnalysis
```

## Input directory structure

The base directory must contain one folder for each protein system. Each system folder must contain a file named `protein.gro`.

```text
protein_database/
├── 1abc/
│   └── protein.gro
├── 2def/
│   └── protein.gro
└── 3ghi/
    └── protein.gro
```

The folder name is used as the system identifier in the output CSV file. Directories without a `protein.gro` file are skipped.

## Usage

Basic usage:

```bash
python estimate_boxes.py --base /path/to/protein_database
```

Example with all available options:

```bash
python estimate_boxes.py \
    --base /path/to/protein_database \
    --cutoff 1.2 \
    --extra 1.0 \
    --water-density 33.4 \
    --block-size 5000 \
    --output protein_box_estimates.csv
```

## Command-line options

| Argument | Default | Description |
|---|---:|---|
| `--base` | `.` | Base directory containing the protein system folders |
| `--cutoff` | `1.2` nm | Nonbonded interaction cutoff |
| `--extra` | `1.0` nm | Additional safety margin |
| `--water-density` | `33.4` nm⁻³ | Approximate number density of water |
| `--block-size` | `5000` | Block size used for the `Dmax` calculation |
| `-o`, `--output` | `protein_box_estimates.csv` | Name of the output CSV file |

## Output

The output CSV contains the following columns:

| Column | Description |
|---|---|
| `PDB` | System identifier obtained from the folder name |
| `Protein_atoms` | Number of atoms selected as protein |
| `Dmax_nm` | Maximum pairwise distance between protein atoms |
| `Conservative_cubic_box_nm` | Estimated cubic box edge length |
| `Box_volume_nm3` | Estimated cubic box volume |
| `Estimated_water_molecules` | Approximate number of water molecules |
| `Estimated_total_atoms` | Protein atoms plus three atoms per water molecule |

## Methodology

### Maximum molecular dimension

The maximum molecular dimension, $D_{\max}$, is defined as the largest distance between any two atoms in the protein:

$$
D_{\max} = \max_{i,j} \left\|\mathbf{r}_i-\mathbf{r}_j\right\|
$$

Pairwise distances are evaluated in blocks to avoid allocating the complete $N \times N$ distance matrix.

### Conservative cubic box

The cubic box edge length is calculated as:

$$
L = D_{\max} + 2(r_{\mathrm{cutoff}} + r_{\mathrm{extra}})
$$

where $r_{\mathrm{cutoff}}$ is the nonbonded interaction cutoff and $r_{\mathrm{extra}}$ is an additional safety margin.

The box volume is:

$$
V_{\mathrm{box}} = L^3
$$

### Protein volume

The approximate protein volume is estimated from its molecular mass:

$$
V_{\mathrm{protein}} \approx M_{\mathrm{protein}} \times 0.00121
$$

where $M_{\mathrm{protein}}$ is expressed in daltons and the resulting volume is expressed in nm³.

### Number of water molecules

The available solvent volume is estimated as:

$$
V_{\mathrm{solvent}} = \max(V_{\mathrm{box}}-V_{\mathrm{protein}},0)
$$

The approximate number of water molecules is then:

$$
N_{\mathrm{water}} \approx V_{\mathrm{solvent}}\rho_{\mathrm{water}}
$$

where the default water number density is 33.4 molecules nm⁻³.

The estimated total number of atoms assumes three atoms per water molecule, as in the TIP3P water model:

$$
N_{\mathrm{total}} = N_{\mathrm{protein}} + 3N_{\mathrm{water}}
$$

## Assumptions and limitations

- The estimated box is cubic and intentionally conservative.
- Calculations are based on the conformation provided in `protein.gro`.
- Only atoms recognized by MDAnalysis through the `protein` selection are included.
- The protein volume is approximated from its molecular mass.
- The total atom estimate assumes a three-site water model such as TIP3P.
- Ions, ligands, cofactors, lipids, nanoparticles, and other non-protein components are not included.
- The final number of atoms generated during system solvation may differ from the estimate.
- The input structure must be whole. Molecules fragmented across periodic boundaries may produce incorrect dimensions.
- The method is intended for system planning and computational-resource estimation rather than exact system construction.

## Dependencies

The `requirements.txt` file should contain:

```text
numpy
MDAnalysis
```

## Author

Jéssica Ierich

## Repository access

This repository is currently intended for internal research use by authorized collaborators.
