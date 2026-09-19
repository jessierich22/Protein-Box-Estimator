#!/usr/bin/env python3

import argparse
import csv
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array


# ============================================================
# SPC216 reference
# ============================================================

SPC216_WATERS = 216
SPC216_BOX_EDGE = 1.86206  # nm
SPC216_VOLUME = SPC216_BOX_EDGE ** 3

WATER_DENSITY = SPC216_WATERS / SPC216_VOLUME


# ============================================================
# Utilities
# ============================================================

def atom_description(atom):
    return (
        f"index={atom.index + 1}, "
        f"name={atom.name}, "
        f"resname={atom.resname}, "
        f"resid={atom.resid}, "
        f"segid={atom.segid}"
    )


def maximum_pairwise_distance(atoms, block_size=2000):
    """
    Calculate maximum atom-atom distance using blocks.

    MDAnalysis coordinates are in Angstrom.
    """

    positions = atoms.positions
    n_atoms = len(atoms)

    max_distance = -1.0
    max_i = None
    max_j = None

    for i_start in range(0, n_atoms, block_size):

        i_end = min(i_start + block_size, n_atoms)
        block_i = positions[i_start:i_end]

        for j_start in range(i_start, n_atoms, block_size):

            j_end = min(j_start + block_size, n_atoms)
            block_j = positions[j_start:j_end]

            distances = distance_array(block_i, block_j)

            # Ignore self-distances and duplicate pairs
            # when comparing the same block.
            if i_start == j_start:
                distances[np.tril_indices_from(distances)] = -1.0

            local_flat_index = np.argmax(distances)
            local_max = distances.flat[local_flat_index]

            if local_max > max_distance:

                local_i, local_j = np.unravel_index(
                    local_flat_index,
                    distances.shape
                )

                max_distance = local_max
                max_i = i_start + local_i
                max_j = j_start + local_j

    return max_distance, atoms[max_i], atoms[max_j]


def molecular_spans(atoms):
    """
    Return molecular X/Y/Z spans in nm.
    """

    positions_nm = atoms.positions / 10.0

    xyz_min = positions_nm.min(axis=0)
    xyz_max = positions_nm.max(axis=0)

    return xyz_max - xyz_min


def estimate_system(box_x, box_y, box_z, n_solute_atoms):

    volume = box_x * box_y * box_z

    waters = int(round(volume * WATER_DENSITY))
    water_atoms = waters * 3
    total_atoms = n_solute_atoms + water_atoms

    return volume, waters, water_atoms, total_atoms


# ============================================================
# Single structure
# ============================================================

def analyze_structure(u, atoms, args):

    print()
    print("=" * 70)
    print(" STRUCTURE / BOX SIZE ANALYSIS")
    print("=" * 70)

    print(f"\nStructure : {args.structure}")
    print(f"Selection : {args.selection}")
    print(f"Atoms     : {len(atoms):,}")

    print("\nCalculating maximum atom-atom distance...")

    max_dist_A, atom1, atom2 = maximum_pairwise_distance(
        atoms,
        block_size=args.block_size
    )

    max_dist_nm = max_dist_A / 10.0

    span_x, span_y, span_z = molecular_spans(atoms)

    total_added_size = args.cutoff + args.extra

    cubic_edge = max_dist_nm + total_added_size

    ortho_x = span_x + total_added_size
    ortho_y = span_y + total_added_size
    ortho_z = span_z + total_added_size

    cubic_stats = estimate_system(
        cubic_edge,
        cubic_edge,
        cubic_edge,
        len(atoms)
    )

    ortho_stats = estimate_system(
        ortho_x,
        ortho_y,
        ortho_z,
        len(atoms)
    )

    print("\n" + "-" * 70)
    print("MAXIMUM ATOM-ATOM DISTANCE")
    print("-" * 70)

    print(f"\nAtom 1 : {atom_description(atom1)}")
    print(f"Atom 2 : {atom_description(atom2)}")

    print(f"\nMaximum distance : {max_dist_nm:.3f} nm")

    print("\n" + "-" * 70)
    print("MOLECULAR SPAN")
    print("-" * 70)

    print(f"\nX : {span_x:.3f} nm")
    print(f"Y : {span_y:.3f} nm")
    print(f"Z : {span_z:.3f} nm")

    print("\n" + "-" * 70)
    print("BOX PARAMETERS")
    print("-" * 70)

    print(f"\nLong-range cutoff : {args.cutoff:.3f} nm")
    print(f"Extra distance    : {args.extra:.3f} nm")
    print(
        f"Total added size  : "
        f"{args.cutoff + args.extra:.3f} nm"
    )

    print("\nFormula:")
    print(
        f"  L = molecular dimension + "
        f"{args.cutoff:.3f} + {args.extra:.3f} nm"
    )

    print("\n" + "-" * 70)
    print("RECOMMENDED BOX DIMENSIONS")
    print("-" * 70)

    print(
        f"\nConservative cubic box (Dmax):\n"
        f"  {cubic_edge:.3f} x "
        f"{cubic_edge:.3f} x "
        f"{cubic_edge:.3f} nm"
    )

    print(
        f"\nAxis-based orthorhombic box:\n"
        f"  {ortho_x:.3f} x "
        f"{ortho_y:.3f} x "
        f"{ortho_z:.3f} nm"
    )

    print("\n" + "-" * 70)
    print("SPC216 WATER ESTIMATE")
    print("-" * 70)

    print(
        f"\nSPC216 reference density : "
        f"{WATER_DENSITY:.3f} waters/nm^3"
    )

    volume, waters, water_atoms, total_atoms = cubic_stats

    print("\nCubic box:")
    print(f"  Volume                  : {volume:.3f} nm^3")
    print(f"  Estimated waters        : {waters:,}")
    print(f"  Estimated water atoms   : {water_atoms:,}")
    print(f"  Solute atoms            : {len(atoms):,}")
    print(f"  Estimated total atoms   : {total_atoms:,}")

    volume, waters, water_atoms, total_atoms = ortho_stats

    print("\nOrthorhombic box:")
    print(f"  Volume                  : {volume:.3f} nm^3")
    print(f"  Estimated waters        : {waters:,}")
    print(f"  Estimated water atoms   : {water_atoms:,}")
    print(f"  Solute atoms            : {len(atoms):,}")
    print(f"  Estimated total atoms   : {total_atoms:,}")

    print(
        "\nNOTE: Water estimates are upper bounds because the volume "
        "occupied by the solute is not subtracted."
    )

    print("\n" + "=" * 70)


# ============================================================
# Trajectory analysis
# ============================================================

def analyze_trajectory(u, atoms, args):

    print()
    print("=" * 70)
    print(" TRAJECTORY MOLECULAR EXTENSION / BOX ANALYSIS")
    print("=" * 70)

    print(f"\nStructure  : {args.structure}")
    print(f"Trajectory : {args.trajectory}")
    print(f"Selection  : {args.selection}")
    print(f"Atoms      : {len(atoms):,}")
    print(f"Stride     : {args.stride}")

    print(
        "\nWARNING: The selected molecule must be whole in every frame."
        "\nA molecule broken across PBC can produce artificial dimensions."
    )

    results = []

    # ------------------------------------------------------------
    # Global Dmax
    # ------------------------------------------------------------

    global_max_dmax = -1.0
    global_max_frame = None
    global_max_time = None
    global_atom1 = None
    global_atom2 = None

    # ------------------------------------------------------------
    # Initial structure
    # ------------------------------------------------------------

    initial_dmax = None
    initial_span_x = None
    initial_span_y = None
    initial_span_z = None

    # ------------------------------------------------------------
    # Maximum axis spans
    # ------------------------------------------------------------

    max_span_x = -1.0
    max_span_y = -1.0
    max_span_z = -1.0

    max_span_x_frame = None
    max_span_y_frame = None
    max_span_z_frame = None

    max_span_x_time = None
    max_span_y_time = None
    max_span_z_time = None

    # ------------------------------------------------------------
    # Clearance screening
    # ------------------------------------------------------------

    minimum_clearance = np.inf
    minimum_clearance_frame = None
    minimum_clearance_time = None
    minimum_clearance_axis = None

    analyzed_frames = 0

    total_added_size = args.cutoff + args.extra

    print("\nAnalyzing trajectory...")

    for ts in u.trajectory[::args.stride]:

        analyzed_frames += 1

        # --------------------------------------------------------
        # Molecular spans
        # --------------------------------------------------------

        span_x, span_y, span_z = molecular_spans(atoms)

        # Track maximum spans

        if span_x > max_span_x:
            max_span_x = span_x
            max_span_x_frame = ts.frame
            max_span_x_time = ts.time

        if span_y > max_span_y:
            max_span_y = span_y
            max_span_y_frame = ts.frame
            max_span_y_time = ts.time

        if span_z > max_span_z:
            max_span_z = span_z
            max_span_z_frame = ts.frame
            max_span_z_time = ts.time

        # --------------------------------------------------------
        # Dmax
        # --------------------------------------------------------

        max_dist_A, atom1, atom2 = maximum_pairwise_distance(
            atoms,
            block_size=args.block_size
        )

        dmax_nm = max_dist_A / 10.0

        # First analyzed frame = reference

        if initial_dmax is None:

            initial_dmax = dmax_nm

            initial_span_x = span_x
            initial_span_y = span_y
            initial_span_z = span_z

        # Changes relative to first frame

        delta_dmax = dmax_nm - initial_dmax

        delta_dmax_percent = (
            delta_dmax / initial_dmax
        ) * 100.0

        required_box = dmax_nm + total_added_size

        # --------------------------------------------------------
        # Global maximum Dmax
        # --------------------------------------------------------

        if dmax_nm > global_max_dmax:

            global_max_dmax = dmax_nm
            global_max_frame = ts.frame
            global_max_time = ts.time

            global_atom1 = (
                atom1.index,
                atom1.name,
                atom1.resname,
                atom1.resid,
                atom1.segid
            )

            global_atom2 = (
                atom2.index,
                atom2.name,
                atom2.resname,
                atom2.resid,
                atom2.segid
            )

        # --------------------------------------------------------
        # Box dimensions
        # --------------------------------------------------------

        if ts.dimensions is not None:

            box_x = ts.dimensions[0] / 10.0
            box_y = ts.dimensions[1] / 10.0
            box_z = ts.dimensions[2] / 10.0

            # Screening metric only:
            #
            # box dimension - molecular Cartesian span

            clearance_x = box_x - span_x
            clearance_y = box_y - span_y
            clearance_z = box_z - span_z

            clearances = np.array([
                clearance_x,
                clearance_y,
                clearance_z
            ])

            local_min_index = np.argmin(clearances)
            local_min_clearance = clearances[local_min_index]

            axes = ["X", "Y", "Z"]

            if local_min_clearance < minimum_clearance:

                minimum_clearance = local_min_clearance
                minimum_clearance_frame = ts.frame
                minimum_clearance_time = ts.time
                minimum_clearance_axis = axes[local_min_index]

        else:

            box_x = np.nan
            box_y = np.nan
            box_z = np.nan

            clearance_x = np.nan
            clearance_y = np.nan
            clearance_z = np.nan

        # --------------------------------------------------------
        # Save frame
        # --------------------------------------------------------

        results.append([
            ts.frame,
            ts.time,
            ts.time / 1000.0,
            dmax_nm,
            delta_dmax,
            delta_dmax_percent,
            required_box,
            span_x,
            span_y,
            span_z,
            box_x,
            box_y,
            box_z,
            clearance_x,
            clearance_y,
            clearance_z
        ])

        if analyzed_frames % 10 == 0:

            print(
                f"\rFrames analyzed: {analyzed_frames}",
                end="",
                flush=True
            )

    print(f"\rFrames analyzed: {analyzed_frames}")

    # ============================================================
    # CSV
    # ============================================================

    with open(args.output, "w", newline="") as csvfile:

        writer = csv.writer(csvfile)

        writer.writerow([
            "frame",
            "time_ps",
            "time_ns",
            "dmax_nm",
            "delta_dmax_nm",
            "delta_dmax_percent",
            "required_cubic_box_nm",
            "span_x_nm",
            "span_y_nm",
            "span_z_nm",
            "box_x_nm",
            "box_y_nm",
            "box_z_nm",
            "clearance_x_nm",
            "clearance_y_nm",
            "clearance_z_nm"
        ])

        writer.writerows(results)

    # ============================================================
    # Extension analysis
    # ============================================================

    max_increase = global_max_dmax - initial_dmax

    max_increase_percent = (
        max_increase / initial_dmax
    ) * 100.0

    required_initial_box = (
        initial_dmax + total_added_size
    )

    required_worst_box = (
        global_max_dmax + total_added_size
    )

    box_growth = (
        required_worst_box - required_initial_box
    )

    # ============================================================
    # Final report
    # ============================================================

    print("\n" + "-" * 70)
    print("MOLECULAR EXTENSION ANALYSIS")
    print("-" * 70)

    print(f"\nInitial Dmax       : {initial_dmax:.3f} nm")
    print(f"Maximum Dmax       : {global_max_dmax:.3f} nm")

    print(
        f"Maximum increase   : "
        f"{max_increase:.3f} nm"
    )

    print(
        f"Relative increase  : "
        f"{max_increase_percent:.2f} %"
    )

    print(
        f"\nMaximum extension frame : "
        f"{global_max_frame}"
    )

    print(
        f"Maximum extension time  : "
        f"{global_max_time / 1000.0:.3f} ns"
    )

    print(f"\nAtom 1 : "
          f"index={global_atom1[0] + 1}, "
          f"name={global_atom1[1]}, "
          f"resname={global_atom1[2]}, "
          f"resid={global_atom1[3]}, "
          f"segid={global_atom1[4]}")

    print(f"Atom 2 : "
          f"index={global_atom2[0] + 1}, "
          f"name={global_atom2[1]}, "
          f"resname={global_atom2[2]}, "
          f"resid={global_atom2[3]}, "
          f"segid={global_atom2[4]}")

    # ------------------------------------------------------------

    print("\n" + "-" * 70)
    print("AXIS EXTENSION ANALYSIS")
    print("-" * 70)

    print("\nInitial spans:")

    print(f"  X : {initial_span_x:.3f} nm")
    print(f"  Y : {initial_span_y:.3f} nm")
    print(f"  Z : {initial_span_z:.3f} nm")

    print("\nMaximum observed spans:")

    print(
        f"  X : {max_span_x:.3f} nm "
        f"(frame {max_span_x_frame}, "
        f"{max_span_x_time / 1000.0:.3f} ns)"
    )

    print(
        f"  Y : {max_span_y:.3f} nm "
        f"(frame {max_span_y_frame}, "
        f"{max_span_y_time / 1000.0:.3f} ns)"
    )

    print(
        f"  Z : {max_span_z:.3f} nm "
        f"(frame {max_span_z_frame}, "
        f"{max_span_z_time / 1000.0:.3f} ns)"
    )

    # ------------------------------------------------------------

    print("\n" + "-" * 70)
    print("BOX REQUIREMENT THROUGH TRAJECTORY")
    print("-" * 70)

    print(
        f"\nCutoff              : "
        f"{args.cutoff:.3f} nm"
    )

    print(
        f"Extra distance      : "
        f"{args.extra:.3f} nm"
    )

    print(
        f"Total added size    : "
        f"{total_added_size:.3f} nm"
    )

    print(
        f"\nInitial required cubic edge : "
        f"{required_initial_box:.3f} nm"
    )

    print(
        f"Worst required cubic edge   : "
        f"{required_worst_box:.3f} nm"
    )

    print(
        f"Increase in required edge   : "
        f"{box_growth:.3f} nm"
    )

    # ------------------------------------------------------------

    print("\n" + "-" * 70)
    print("CARTESIAN BOX CLEARANCE SCREENING")
    print("-" * 70)

    if np.isfinite(minimum_clearance):

        print(
            f"\nMinimum clearance : "
            f"{minimum_clearance:.3f} nm"
        )

        print(
            f"Axis              : "
            f"{minimum_clearance_axis}"
        )

        print(
            f"Frame             : "
            f"{minimum_clearance_frame}"
        )

        print(
            f"Time              : "
            f"{minimum_clearance_time / 1000.0:.3f} ns"
        )

        print(
            "\nDefinition:"
            "\n  clearance = box dimension - molecular Cartesian span"
        )

        print(
            "\nIMPORTANT:"
            "\n  This is a geometric screening metric only."
            "\n  It is NOT the minimum atom-atom distance between"
            "\n  the protein and its periodic image."
        )

        if minimum_clearance < 0:

            screening_status = "CRITICAL"

        elif minimum_clearance < args.cutoff:

            screening_status = "FLAG"

        else:

            screening_status = "PASS"

        print(
            f"\nScreening status : "
            f"{screening_status}"
        )

        if screening_status == "FLAG":

            print(
                "At least one frame should be inspected with an "
                "explicit protein-periodic-image distance analysis."
            )

        elif screening_status == "CRITICAL":

            print(
                "The Cartesian molecular span exceeded a box "
                "dimension. Check PBC treatment immediately."
            )

        else:

            print(
                "No frame was flagged by the Cartesian span "
                "screening criterion."
            )

    else:

        print(
            "\nNo valid box dimensions were found."
        )

    # ------------------------------------------------------------

    print("\n" + "-" * 70)
    print("OUTPUT")
    print("-" * 70)

    print(
        f"\nPer-frame data written to:"
        f"\n  {args.output}"
    )

    print(
        "\nCSV includes:"
        "\n  Dmax"
        "\n  Delta Dmax"
        "\n  Relative Dmax change"
        "\n  Required cubic box"
        "\n  X/Y/Z spans"
        "\n  Actual box dimensions"
        "\n  Cartesian clearance screening"
    )

    print("\n" + "=" * 70)


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Estimate simulation box dimensions from a structure "
            "or analyze molecular extension throughout a trajectory."
        )
    )

    parser.add_argument(
        "-s",
        "--structure",
        required=True,
        help="Structure/topology file: PDB, GRO, TPR, etc."
    )

    parser.add_argument(
        "-f",
        "--trajectory",
        default=None,
        help="Optional trajectory file, e.g. XTC"
    )

    parser.add_argument(
        "--selection",
        default="protein",
        help='MDAnalysis selection (default: "protein")'
    )

    parser.add_argument(
        "--cutoff",
        type=float,
        default=1.2,
        help="Long-range cutoff in nm (default: 1.2)"
    )

    parser.add_argument(
        "--extra",
        type=float,
        default=1.0,
        help="Additional safety distance in nm (default: 1.0)"
    )

    parser.add_argument(
        "--block-size",
        type=int,
        default=2000,
        help="Atoms per distance-calculation block (default: 2000)"
    )

    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Analyze every Nth trajectory frame (default: 1)"
    )

    parser.add_argument(
        "-o",
        "--output",
        default="box_analysis.csv",
        help="Trajectory CSV output"
    )

    args = parser.parse_args()

    # ------------------------------------------------------------
    # Load
    # ------------------------------------------------------------

    if args.trajectory:

        u = mda.Universe(
            args.structure,
            args.trajectory
        )

    else:

        u = mda.Universe(
            args.structure
        )

    atoms = u.select_atoms(args.selection)

    if len(atoms) < 2:

        raise ValueError(
            f"Selection '{args.selection}' contains fewer than 2 atoms."
        )

    # ------------------------------------------------------------
    # Run
    # ------------------------------------------------------------

    if args.trajectory:

        analyze_trajectory(
            u,
            atoms,
            args
        )

    else:

        analyze_structure(
            u,
            atoms,
            args
        )


if __name__ == "__main__":
    main()
