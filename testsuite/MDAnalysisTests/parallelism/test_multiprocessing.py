"""Test that MDAnalysis plays nicely with multiprocessing

"""
import multiprocessing
import numpy as np
import pytest
import pickle

import MDAnalysis as mda
from MDAnalysis.coordinates.core import get_reader_for
from MDAnalysisTests.datafiles import (
    AUX_XVG,
    CRD,
    PSF, DCD,
    DMS,
    DLP_CONFIG,
    DLP_HISTORY,
    INPCRD,
    GMS_ASYMOPT,
    GRO,
    LAMMPSdata_mini,
    mol2_molecules,
    MMTF,
    NCDF,
    PDB_small, PDB_multiframe,
    PDBQT_input,
    PQR,
    TRR,
    TRJ,
    TRZ,
    XTC,
    XPDB_small,
    XYZ_mini, XYZ,
)

from numpy.testing import assert_equal


@pytest.fixture(params=[
    (PSF, DCD),
    (GRO, XTC),
    (PDB_multiframe,),
    (XYZ,),
])
def u(request):
    if len(request.param) == 1:
        f = request.param[0]
        return mda.Universe(f)
    else:
        top, trj = request.param
        return mda.Universe(top, trj)

# Define target functions here
# inside test functions doesn't work
def cog(u, ag, frame_id):
    u.trajectory[frame_id]

    return ag.center_of_geometry()


def getnames(u, ix):
    # Check topology stuff works
    return u.atoms[ix].name


def test_multiprocess_COG(u):
    ag = u.atoms[10:20]

    ref = np.array([cog(u, ag, i)
                    for i in range(4)])

    p = multiprocessing.Pool(2)
    res = np.array([p.apply(cog, args=(u, ag, i))
                    for i in range(4)])
    p.close()
    assert_equal(ref, res)


def test_multiprocess_names(u):
    ref = [getnames(u, i)
           for i in range(10)]

    p = multiprocessing.Pool(2)
    res = [p.apply(getnames, args=(u, i))
                   for i in range(10)]
    p.close()

    assert_equal(ref, res)

@pytest.fixture(params=[
    # formatname, filename
    ('CRD', CRD, dict()),
    ('DATA', LAMMPSdata_mini, dict(n_atoms=1)),
    ('DCD', DCD, dict()),
    ('DMS', DMS, dict()),
    ('CONFIG', DLP_CONFIG, dict()),
    ('HISTORY', DLP_HISTORY, dict()),
    ('INPCRD', INPCRD, dict()),
    ('GMS', GMS_ASYMOPT, dict()),
    ('GRO', GRO, dict()),
    ('MMTF', MMTF, dict()),
    ('MOL2', mol2_molecules, dict()),
    ('PDB', PDB_small, dict()),
    ('PQR', PQR, dict()),
    ('PDBQT', PDBQT_input, dict()),
    ('TRR', TRR, dict()),
    ('TRZ', TRZ, dict(n_atoms=8184)),
    ('TRJ', TRJ, dict(n_atoms=252)),
    ('XTC', XTC, dict()),
    ('XPDB', XPDB_small, dict()),
    ('XYZ', XYZ_mini, dict()),
    ('NCDF', NCDF, dict()),
    ('memory', np.arange(60).reshape(2, 10, 3).astype(np.float64), dict()),
    ('CHAIN', [GRO, GRO, GRO], dict()),
])
def ref_reader(request):
    fmt_name, filename, extras = request.param

    r = get_reader_for(filename, format=fmt_name)(filename, **extras)
    try:
        yield r
    finally:
        # make sure file handle is closed afterwards
        r.close()

def test_readers_pickle(ref_reader):
    ps = pickle.dumps(ref_reader)

    reanimated = pickle.loads(ps)

    assert len(ref_reader) == len(reanimated)