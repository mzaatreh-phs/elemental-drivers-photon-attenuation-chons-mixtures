# Building XCOM

Vendored from NIST's public-domain XCOM v3.1 distribution (M.J. Berger &
J.H. Hubbell, 1999; `~/Downloads/XCOM.tar.gz` on this machine).

```
gfortran -w -std=legacy -o xcom XCOM.f
```

`-std=legacy` is required: this is fixed-form Fortran 77 with old-style
continuation markers gfortran otherwise rejects. `-w` silences legacy-syntax
warnings (none of them are correctness issues -- the binary reproduces
textbook NIST values exactly, see epixs/xcom_engine.py's docstring).

The binary must be run with its working directory set to this folder
(`vendor/xcom/`), since it opens its per-element `MDATX3.NNN` data files by
bare relative name -- `epixs/xcom_engine.py` handles this via `cwd=` on the
subprocess call.

## Precision patch (FORMAT 820/840)

The stock program prints its cross-section table with `1PE10.3` (4
significant figures). That's fine for the coefficients themselves, but
`epixs/zeff.py`'s Zeq determination needs the RATIO of two of those columns
(incoherent / total-without-coherent), and at Compton-dominated energies
(roughly 0.2-1.5 MeV) photoelectric absorption -- the only thing that still
differentiates that ratio between elements once coherent scattering is
excluded -- is many orders of magnitude smaller than incoherent scattering.
At 4 significant figures the two columns print identically, the ratio
rounds to exactly 1.0, and Zeq becomes undefined/unstable (verified: it was
jumping from 7.3 to 4.0 to 6.6 across adjacent energies for a real test
mixture). Patched both output FORMAT statements from `1PE10.3` to
`1PE20.12` (full double precision) to fix this; `epixs/xcom_engine.py`'s
parser regex was updated to match. Re-run `gfortran -w -std=legacy -o xcom
XCOM.f` after modifying XCOM.f to rebuild with this patch.

`run` (the top-level launcher) compiles this automatically on first use if
`xcom` isn't already built.
