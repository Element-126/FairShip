# -*- coding: utf-8 -*-

from __future__ import absolute_import, division
from future.utils import viewitems

from collections import OrderedDict
import numpy as np

from scalar_portal import *
from scalar_portal.data.particles import *


_placeholder_particles = {
     0: ('r0', 'void', +9900080),
    +1: ('r+', 'r-', +9900081),
    -1: ('r-', 'r+', -9900081)
}

def _pythia_dummy_particle_strings():
    return OrderedDict((name, format_pythia_particle_string(
        pdg_id=pdg_id, name=name, antiname=antiname, spin_type=1,
        charge_type=3*charge, mass=0., lifetime_si=0., new=True,
        may_decay=False, is_visible=False))
            for charge, (name, antiname, pdg_id) in viewitems(_placeholder_particles)
            if charge >= 0)

def _pythia_disable_parent_sm_channels(particle):
    string = '{}:onMode = off'.format(get_pdg_id(particle))
    return OrderedDict([(particle, string)])

def _pythia_complementary_string(parent, branching_ratio):
    '''
    Dummy channels are used to prevent PYTHIA from rescaling the total
    branching ratio of each heavy meson to 1, hence preserving the relative
    production widths for the scalar.

    They cause the heavy mesons to decay to some placeholder state with a
    branching ratio complementary (1-BR) to the process it wraps. These
    placeholder states can then be discarded during the analysis.
    '''
    charge = get_charge(parent)
    rzero = _placeholder_particles[0]
    _, _, placeholder_id = _placeholder_particles[charge]
    return format_pythia_string(get_pdg_id(parent), [placeholder_id, rzero[2]], branching_ratio)

def _root_add_dummy_particles(pdg):
    'Add the placeholder particles to the ROOT database.'
    rzero = _placeholder_particles[0]
    rplus = _placeholder_particles[1]
    pdg.AddParticle(rzero[0], 'Reject (neutral)', 0.0, True, 0.0,  0.0, 'reject', rzero[2])
    pdg.AddParticle(rplus[0], 'Reject (charged)', 0.0, True, 0.0, +1.0, 'reject', rplus[2])
    pdg.AddAntiParticle(rplus[1], -rplus[2])

class RescaledProductionBranchingRatios(ProductionBranchingRatios):
    '''
    Specialization of `scalar_portal.ProductionBranchingRatios` with extra
    quirks to rescale the branching ratios for FairShip usage.
    '''
    def __init__(self, *args, **kwargs):
        super(RescaledProductionBranchingRatios, self).__init__(*args, **kwargs)
        self._parent_particles = set(ch._parent for ch in self._channels.values())
        self._channels_by_parent = OrderedDict(
            (parent, [str(ch) for ch in self._channels.values() if ch._parent == parent])
            for parent in self._parent_particles)
        self._sum_br = OrderedDict(
            (parent, sum(self._br[ch] for ch in self._channels_by_parent[parent]))
            for parent in self._parent_particles)
        # Ref: https://stackoverflow.com/a/39279912
        self._max_sum_br = np.fmax.reduce(self._sum_br.values())
        self._rescaled_br = OrderedDict(
            (st, br / self._max_sum_br) for st, br in viewitems(self._br))
        self._complementary_channels = OrderedDict(
            (parent, 1 - sum_br/self._max_sum_br)
            for parent, sum_br in viewitems(self._sum_br))

    @property
    def branching_ratios(self):
        return self._rescaled_br

    @property
    def maximum_total_branching_ratio(self):
        return self._max_sum_br

    @property
    def scaling_factor(self):
        return 1 / self.maximum_total_branching_ratio

    def pythia_strings(self):
        all_strs = OrderedDict()
        # Add dummy particle definitions
        dummy_particle_strs = _pythia_dummy_particle_strings()
        all_strs.update(dummy_particle_strs)
        # Reset parents to remove their SM channels
        for parent in self._parent_particles:
            new_parent_str = _pythia_disable_parent_sm_channels(parent)
            all_strs.update(new_parent_str)
        # Add BSM channels
        base_strs = super(RescaledProductionBranchingRatios, self).pythia_strings()
        all_strs.update(base_strs)
        # Add dummy channels
        complementary_strs = OrderedDict(
            ('{} (rejected)'.format(parent), _pythia_complementary_string(parent, br))
            for parent, br in viewitems(self._complementary_channels)
            if br > 0)
        all_strs.update(complementary_strs)
        return all_strs

class FairShipBranchingRatiosResult(BranchingRatiosResult):
    '''
    Specialization of `scalar_portal.BranchingRatiosResult` with an extra
    method to setup ROOT.
    '''
    def __init__(self, *args, **kwargs):
        super(FairShipBranchingRatiosResult, self).__init__(*args, **kwargs)

    def root_add_particles(self, pdg):
        if self.decays._mS.ndim > 0:
            raise(ValueError('Can only set up ROOT for a single scalar mass.'))
        pdg.AddParticle(
            'S', 'Scalar', self.decays._mS, False, self.decays.total_width, 0.0,
            'hiddensector', self.decays._scalar_id)
        _root_add_dummy_particles(pdg)

class FairShipScalarModel(Model):
    '''
    Specialization of `scalar_portal.Model` with extra FairShip-specific
    features.
    '''
    def __init__(self, scalar_id=9900015):
        super(FairShipScalarModel, self).__init__(scalar_id)

    def compute_branching_ratios(self, mass, couplings=None, ignore_invalid=False, **kwargs):
        '''
        Compute the decay and (rescaled) production branching ratios of the
        scalar particle, and return a `BranchingRatiosResult` object containing
        the result.
        '''
        if couplings is None:
            couplings = kwargs
        prod_channels  = self.production.get_active_processes()
        decay_channels = self.decays.get_active_processes()
        prod_br  = RescaledProductionBranchingRatios(
            prod_channels , mass, couplings, ignore_invalid, scalar_id=self.scalar_pdg_id)
        decay_br = DecayBranchingRatios(
            decay_channels, mass, couplings, ignore_invalid, scalar_id=self.scalar_pdg_id)
        res = FairShipBranchingRatiosResult(prod_br, decay_br)
        return res
