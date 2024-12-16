# README

**Functions implemented in this branch will be PR to the official version gradually.**

**This package is modified by Jia-Xin Zhu on the basis of [MDAnalysis v2.3.0](https://github.com/MDAnalysis/mdanalysis/tree/package-2.3.0).**

## selection based on relative positions

Explanation...
Example...

## XYZ writer

Probably we can use other formats, such as `PDB`. 

Support [ExtXYZ](https://github.com/libAtoms/extxyz) format.

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

# Deprecated

## LAMMPSDUMP reader

> Use `additional_columns` in [`MDAnalysis.coordinates.LAMMPS.DumpReader`](https://docs.mdanalysis.org/2.8.0/documentation_pages/coordinates/LAMMPS.html#MDAnalysis.coordinates.LAMMPS.DumpReader)

Read forces, velocities and charges from lammpsdump file.

## XYZ topology

> Use [`Universe.add_TopologyAttr`](https://userguide.mdanalysis.org/stable/topology_system.html#adding-topologyattrs)

If there are atomic charges info included in extended xyz file, these charges will be read into the `Universe` object.
