#!/usr/bin/env python
"""
Piotr Zarzycki
email: ppzarzycki@lbl.gov
Lawrence Berkeley National Laboratory 

1 read a Gaussian formatted checkpoint (fchk)
2 extract Hessian
3 read atomic masses from a separate text file
4 compute harmonic vibrational frequencies in cm^-1.

Usage:
    python harmfreq_from_fchk.p  <fchk_file> <masses_file>

to get formatted checkpoint use Gaussian formchk file 
    formchk <chk_file> 
it will generate <fchk_file> 
for example 
    formchk CO3.chk 
producedd CO3.fchk

the masses file one mass per line 
for example for CO3.fchk  
12.0000
15.99491  
15.99491  
15.99491 

2026-03
"""

import sys
import numpy as np

# Constants
HARTREE_TO_JOULE = 4.3597447222071e-18       # 1 Hartree in Joules
BOHR_TO_METER = 5.29177210903e-11            # 1 Bohr in meters
AMU_TO_KG = 1.66053906660e-27               # 1 amu in kg
SPEED_OF_LIGHT_CM = 2.99792458e10            # speed of light in cm/s
PI = np.pi


def read_fchk(filename):
    natoms = None
    atomic_numbers = []
    coordinates = []
    hessian = []

    with open(filename, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]

        # --- Number of atoms ---
        if line.startswith('Number of atoms'):
            natoms = int(line.split()[-1])
            i += 1
            continue

        # --- Atomic numbers ---
        if line.startswith('Atomic numbers'):
            n_vals = int(line.split()[-1])
            i += 1
            atomic_numbers = _read_integer_array(lines, i, n_vals)
            i += _num_lines_for_values(n_vals, per_line=6)
            continue

        # --- Current cartesian coordinates ---
        if line.startswith('Current cartesian coordinates'):
            n_vals = int(line.split()[-1])
            i += 1
            coordinates = _read_real_array(lines, i, n_vals)
            i += _num_lines_for_values(n_vals, per_line=5)
            continue

        # --- Cartesian Force Constants (Hessian) ---
        if line.startswith('Cartesian Force Constants'):
            n_vals = int(line.split()[-1])
            i += 1
            hessian = _read_real_array(lines, i, n_vals)
            i += _num_lines_for_values(n_vals, per_line=5)
            continue

        i += 1

    if natoms is None:
        raise ValueError("Could not find 'Number of atoms' in the fchk file.")
    if len(hessian) == 0:
        raise ValueError("Could not find 'Cartesian Force Constants' in the fchk file.")

    atomic_numbers = np.array(atomic_numbers, dtype=int)
    coordinates = np.array(coordinates).reshape(natoms, 3)

    return natoms, atomic_numbers, coordinates, np.array(hessian)


def _read_integer_array(lines, start, n_vals):
    """Read n_vals integers starting from line index `start`."""
    values = []
    i = start
    while len(values) < n_vals:
        values.extend(int(x) for x in lines[i].split())
        i += 1
    return values[:n_vals]


def _read_real_array(lines, start, n_vals):
    """Read n_vals real numbers starting from line index `start`."""
    values = []
    i = start
    while len(values) < n_vals:
        values.extend(float(x) for x in lines[i].split())
        i += 1
    return values[:n_vals]


def _num_lines_for_values(n_vals, per_line):
    """How many lines are needed to store n_vals values at per_line values per line."""
    return (n_vals + per_line - 1) // per_line


def read_masses(filename):
    masses = []
    with open(filename, 'r') as f:
        for line in f:
            stripped = line.strip()
            if stripped == '' or stripped.startswith('#'):
                continue
            masses.append(float(stripped.split()[0]))
    return np.array(masses)


def lower_triangle_to_full(lt, n):
    mat = np.zeros((n, n))
    idx = 0
    for i in range(n):
        for j in range(i + 1):
            mat[i, j] = lt[idx]
            mat[j, i] = lt[idx]
            idx += 1
    return mat


def compute_mass_weighted_hessian(hessian_full, masses):
    natoms = len(masses)
    n3 = 3 * natoms

    inv_sqrt_m = np.zeros(n3)
    for a in range(natoms):
        for k in range(3):
            inv_sqrt_m[3 * a + k] = 1.0 / np.sqrt(masses[a])

    mw_hessian = hessian_full * np.outer(inv_sqrt_m, inv_sqrt_m)
    return mw_hessian


def eigenvalues_to_frequencies(eigenvalues):
    conv = HARTREE_TO_JOULE / (BOHR_TO_METER ** 2 * AMU_TO_KG)

    frequencies = np.zeros_like(eigenvalues)
    for i, ev in enumerate(eigenvalues):
        if ev >= 0:
            freq_si = np.sqrt(ev * conv)
        else:
            freq_si = -np.sqrt(-ev * conv) 

        frequencies[i] = freq_si / (2.0 * PI * SPEED_OF_LIGHT_CM)

    return frequencies


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print("Error: Please provide an fchk file and a masses file.")
        sys.exit(1)

    fchk_file = sys.argv[1]
    masses_file = sys.argv[2]

    # ---- Read data ----
    print(f"Reading formatted checkpoint file: {fchk_file}")
    natoms, atomic_numbers, coordinates, hessian_lt = read_fchk(fchk_file)
    print(f"  Number of atoms: {natoms}")

    print(f"Reading masses from: {masses_file}")
    masses = read_masses(masses_file)
    print(f"  Number of masses read: {len(masses)}")

    if len(masses) != natoms:
        raise ValueError(
            f"Number of masses ({len(masses)}) does not match "
            f"number of atoms in fchk ({natoms})."
        )

    n3 = 3 * natoms
    expected_lt = n3 * (n3 + 1) // 2
    if len(hessian_lt) != expected_lt:
        raise ValueError(
            f"Expected {expected_lt} lower-triangular Hessian elements, "
            f"but found {len(hessian_lt)}."
        )

    print("Building full Hessian from lower-triangular storage...")
    hessian_full = lower_triangle_to_full(hessian_lt, n3)

    print("Mass-weighting the Hessian...")
    mw_hessian = compute_mass_weighted_hessian(hessian_full, masses)

    print("Diagonalizing mass-weighted Hessian...")
    eigenvalues, eigenvectors = np.linalg.eigh(mw_hessian)

    frequencies = eigenvalues_to_frequencies(eigenvalues)

    print("\n" + "=" * 60)
    print("  Harmonic Vibrational Frequencies")
    print("=" * 60)
    print(f"  {'Mode':>4s}  {'Eigenvalue (au)':>18s}  {'Frequency (cm-1)':>18s}")
    print("-" * 60)

    is_linear = _check_linear(coordinates, natoms)
    n_tr = 5 if is_linear else 6
    n_vib = n3 - n_tr

    for idx in range(n3):
        freq = frequencies[idx]
        ev = eigenvalues[idx]
        label = ""
        if idx < n_tr:
            label = "  (Translation)"
        elif freq < 0:
            label = "  (imaginary)"
        print(f"  {idx + 1:4d}  {ev:18.8e}  {freq:20.8f}{label}")

    print("-" * 60)
    print(f"\n  Molecule is {'linear' if is_linear else 'non-linear'}.")
    print(f"  Expected {n_tr} translational/rotational modes "
          f"and {n_vib} vibrational modes.\n")

    # Print just the vibrational frequencies
    vib_freqs = frequencies[n_tr:]
    print("  Vibrational frequencies (cm-1):")
    for idx, freq in enumerate(vib_freqs):
        sign = "" if freq >= 0 else "i"
        print(f"    {idx + 1:4d}:  {abs(freq):20.10f} {sign}")
    print()

    return frequencies


def _check_linear(coords, natoms):
    if natoms <= 2:
        return True

    ref = coords[0]
    direction = None
    for i in range(1, natoms):
        diff = coords[i] - ref
        norm = np.linalg.norm(diff)
        if norm > 1e-10:
            direction = diff / norm
            break

    if direction is None:
        return True  

    for i in range(2, natoms):
        diff = coords[i] - ref
        norm = np.linalg.norm(diff)
        if norm < 1e-10:
            continue
        diff_normalized = diff / norm
        cross = np.linalg.norm(np.cross(direction, diff_normalized))
        if cross > 1e-6:
            return False
    return True


if __name__ == '__main__':
    main()
