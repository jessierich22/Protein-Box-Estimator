#!/usr/bin/env python3

import os
import csv
import argparse
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import self_distance_array


def calculate_dmax(atoms, block_size=5000):
    """
    Calculate the maximum pairwise atomic distance using blocks
    to avoid allocating the complete NxN distance matrix.
    """

    positions = atoms.positions.astype(np.float64)
    n_atoms = len(positions)

    dmax = 0.0

    for i in range(0, n_atoms, block_size):

        block_i = positions[i:i + block_size]

        # Distances inside the same block
        if len(block_i) > 1:
            distances = self_distance_array(block_i)

            if len(distances) > 0:
                dmax = max(dmax, distances.max())

        # Distances against subsequent blocks
        for j in range(i + block_size, n_atoms, block_size):

            block_j = positions[j:j + block_size]

            diff = block_i[:, None, :] - block_j[None, :, :]
            distances = np.sqrt(np.sum(diff * diff, axis=2))

            dmax = max(dmax, distances.max())

    # MDAnalysis coordinates are in Angstrom
    return dmax / 10.0


def main():

    parser = argparse.ArgumentParser(
        description="Estimate recommended cubic boxes for a protein database."
    )

    parser.add_argument(
        "--base",
        default=".",
        help="Base directory containing the PDB folders.",
    )

    parser.add_argument(
        "--cutoff",
        type=float,
        default=1.2,
        help="Nonbonded cutoff in nm.",
    )

    parser.add_argument(
        "--extra",
        type=float,
        default=1.0,
        help="Additional safety margin in nm.",
    )

    parser.add_argument(
        "--water-density",
        type=float,
        default=33.4,
        help="Approximate number of water molecules per nm^3.",
    )

    parser.add_argument(
        "--block-size",
        type=int,
        default=5000,
        help="Block size for Dmax calculation.",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="protein_box_estimates.csv",
        help="Output CSV file.",
    )

    args = parser.parse_args()

    results = []

    folders = sorted(
        d for d in os.listdir(args.base)
        if os.path.isdir(os.path.join(args.base, d))
    )

    print("=" * 80)
    print(" PROTEIN DATABASE - recommended BOX ESTIMATION")
    print("=" * 80)
    print(f"Base directory : {args.base}")
    print(f"Cutoff         : {args.cutoff:.2f} nm")
    print(f"Extra margin   : {args.extra:.2f} nm")
    print(f"Output         : {args.output}")
    print("=" * 80)

    for pdb in folders:

        protein_file = os.path.join(
            args.base,
            pdb,
            "protein.gro",
        )

        if not os.path.isfile(protein_file):
            continue

        print()
        print(f"[{pdb}] Reading {protein_file}")

        try:
            u = mda.Universe(protein_file)

            protein = u.select_atoms("protein")

            if protein.n_atoms == 0:
                print("  WARNING: no protein atoms detected.")
                continue

            n_protein = protein.n_atoms

            # -----------------------------------------------
            # Maximum protein dimension
            # -----------------------------------------------

            dmax = calculate_dmax(
                protein,
                block_size=args.block_size,
            )

            # -----------------------------------------------
            # recommended cubic box
            # -----------------------------------------------

            margin = args.cutoff + args.extra

            box_edge = dmax + 2.0 * margin

            box_volume = box_edge ** 3

            # -----------------------------------------------
            # Approximate protein volume
            #
            # Typical partial specific volume:
            # ~0.73 cm3/g
            #
            # This corresponds roughly to:
            # ~0.00121 nm3/Da
            #
            # We can estimate protein volume from masses.
            # -----------------------------------------------

            protein_mass_da = np.sum(protein.masses)

            protein_volume = (
                protein_mass_da * 0.00121
            )

            solvent_volume = max(
                box_volume - protein_volume,
                0.0,
            )

            # -----------------------------------------------
            # Estimate water molecules
            # -----------------------------------------------

            water_molecules = int(
                round(
                    solvent_volume *
                    args.water_density
                )
            )

            # TIP3P = 3 atoms / water
            water_atoms = water_molecules * 3

            estimated_total_atoms = (
                n_protein + water_atoms
            )

            # -----------------------------------------------
            # Store
            # -----------------------------------------------

            results.append(
                {
                    "PDB": pdb,
                    "Protein_atoms": n_protein,
                    "Dmax_nm": round(dmax, 3),
                    "recommended_cubic_box_nm": round(
                        box_edge, 3
                    ),
                    "Box_volume_nm3": round(
                        box_volume, 1
                    ),
                    "Estimated_water_molecules": water_molecules,
                    "Estimated_total_atoms": estimated_total_atoms,
                }
            )

            print(
                f"  Protein atoms : {n_protein:,}"
            )

            print(
                f"  Dmax          : {dmax:.3f} nm"
            )

            print(
                f"  Cubic box     : {box_edge:.3f} nm"
            )

            print(
                f"  Est. atoms    : {estimated_total_atoms:,}"
            )

        except Exception as error:

            print(
                f"  ERROR: {error}"
            )

    # ============================================================
    # Write CSV
    # ============================================================

    output_path = os.path.join(
        args.base,
        args.output,
    )

    fieldnames = [
        "PDB",
        "Protein_atoms",
        "Dmax_nm",
        "recommended_cubic_box_nm",
        "Box_volume_nm3",
        "Estimated_water_molecules",
        "Estimated_total_atoms",
    ]

    with open(
        output_path,
        "w",
        newline="",
    ) as csvfile:

        writer = csv.DictWriter(
            csvfile,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(results)

    print()
    print("=" * 80)
    print(" FINISHED")
    print("=" * 80)
    print(f"Proteins analyzed : {len(results)}")
    print(f"Output            : {output_path}")
    print("=" * 80)


if __name__ == "__main__":
    main()
