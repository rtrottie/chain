#!/usr/bin/env python

from pylada.vasp.functional import Vasp

class NEBVasp(Vasp):

    def __init__(self, copy=None, species=None, kpoints=None, **kwargs):
        super().__init__(copy, species, kpoints, kwargs)
        self.restart = kwargs.get('restart', None)
        self.initial = kwargs.get('initial', None)
        self.final   = kwargs.get('final'  , None)

    def bringup(self, structure, outdir, **kwargs):
        if self.restart:
            Vasp.