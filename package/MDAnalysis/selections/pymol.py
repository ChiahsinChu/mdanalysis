# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
PyMOL selections
=================

Write :class:`MDAnalysis.core.groups.AtomGroup` selection to a
script `pml`_ file that defines PyMOL_ atomselect macros. To be used
in PyMOL like this::

  @macros.pml

The selections should appear in the user interface.

.. _PyMOL: http://www.pymol.org
.. _pml: http://pymol.sourceforge.net/newman/user/S0210start_cmds.html#6_5_1

.. autoclass:: SelectionWriter
   :inherited-members:
"""
from . import base


class SelectionWriter(base.SelectionWriterBase):
    format = ["PyMol", "pml"]
    ext = "pml"
    continuation = '\\'  # quoted backslash!
    commentfmt = "# %s"
    default_numterms = 6

    def _translate(self, atoms, **kwargs):
        # PyMol index is 1-based
        def _index(atom):
            return "index {0:d}".format((atom.index + 1))

        return base.join(atoms, ' |', _index)

    def _write_head(self, out, **kwargs):
        out.write(self.comment("MDAnalysis PyMol selection"))
        out.write("select {name!s}, ".format(**kwargs) + self.continuation + '\n')
