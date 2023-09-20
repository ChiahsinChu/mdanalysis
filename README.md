# README

**This package is modified by Jia-Xin Zhu on the basis of [MDAnalysis v2.3.0](https://github.com/MDAnalysis/mdanalysis/tree/package-2.3.0).**

## selection based on relative positions

Explanation...
Example...

## LAMMPSDUMP reader

Explanation...

## XYZ writer

If writing XYZ file with `mda.Writer.write` and the dimension is not None, the cell parameters and the pbc info will be written in the second row. This allows the cell parameters can be read by `ASE` from the output XYZ file.

Example:

```python
import MDAnalysis as mda

u = mda.Universe("topo.psf",
                 "test.lammpstrj",
                 topology_format="PSF",
                 format="LAMMPSDUMP")
with mda.Writer("output.xyz") as W:
    for ts in u.trajectory:
        W.write(u)
```

## XYZ topology

If there are atomic charges info included in extended xyz file, these charges will be read into the `Universe` object.
