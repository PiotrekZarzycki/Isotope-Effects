# Isotope-Effects:

1. Recalculating harmonic vibrational frequencies
2. Diatomics 
3. Polyatomics
4. Path-Integral


### 1. Recalculating harmonic vibrational frequencies
Recalculating harmonic vibrational frequencies for different isotopologues using a previously computed numerical Hessian 
(second derivatives of energy with respect to atomic coordinates) by reading the Hessian from a file and reconstructing the
mass-weighted Hessian matrix with updated atomic masses and re-diagonalizing it.

#### Gaussian-formatted checkpoint
use 
```formchk <chk_file>``` 
to generate formated checkpoint, and then use harmfreq_from_fchk.py to calculate frequencies for given masses read from an additional input file (provide mass for each atom at separate line), as many as atoms in the system)
Usage: 
``` harmfreq_from_fchk.py <fchk_file> <file_with_masses>```
example for CO32- molecule in the Examples/Gaussian

```python harmfreq_from_fchk.py CO3.fchk masses.txt ```

with masses file (`masses.txt`):
```
12.00000
15.99491  
15.99491  
15.99491  
```
At the end of the output you will see: 
```
============================================================
  Harmonic Vibrational Frequencies
============================================================
  Mode     Eigenvalue (au)    Frequency (cm-1)
------------------------------------------------------------
     1     -4.15276475e-11           -0.03312630  (Translation)
     2      3.33072318e-12            0.00938153  (Translation)
     3      1.31796050e-11            0.01866188  (Translation)
     4      1.59113419e-09            0.20504901  (Translation)
     5      1.20965381e-08            0.56537269  (Translation)
     6      2.13358704e-08            0.75086085  (Translation)
     7      1.56723703e-02          643.53418152
     8      1.56723776e-02          643.53433233
     9      2.79701146e-02          859.70885592
    10      3.86311833e-02         1010.35334358
    11      6.40693051e-02         1301.15574622
    12      6.40701774e-02         1301.16460420
------------------------------------------------------------

  Molecule is non-linear.
  Expected 6 translational/rotational modes and 6 vibrational modes.

  Vibrational frequencies (cm-1):
       1:        643.5341815196 
       2:        643.5343323318 
       3:        859.7088559183 
       4:       1010.3533435837 
       5:       1301.1557462217 
       6:       1301.1646042030
```


## 2. Diatomics 
### 2.1. Dunham Analysis of Diatomic Potential Energy Surfaces
Fits a Dunham polynomial to ab initio Potential Energy Surface (PES) data for diatomic molecule, extracts Dunham Y-parameters,
converts to spectroscopic constants, computes energy levels for multiple isotopologues

Input files are  PES data file in format R/Å and E/Hartree (two columns) and mass file with lines `label  mass1_amu  mass2_amu`




