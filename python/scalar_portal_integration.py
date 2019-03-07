# -*- coding: utf-8 -*-

from __future__ import absolute_import, division
from future.utils import viewitems

import numpy as np

from scalar_portal import *


class RescaledProductionBranchingRatios(ProductionBranchingRatios):
    '''
    Specialization of `scalar_portal.ProductionBranchingRatios` with extra
    quirks to rescale the branching ratios for FairShip usage.

    TODO: add dummy particles.
    '''
    def __init__(self, *args, **kwargs):
        super(RescaledProductionBranchingRatios, self).__init__(*args, **kwargs)
        # Ref: https://stackoverflow.com/a/39279912
        self._max_br = np.fmax.reduce(self._br.values())
        self._rescaled_br = {
            st: br / self._max_br for st, br in viewitems(self._br)}

    @property
    def branching_ratios(self):
        return self._rescaled_br

    @property
    def maximum_branching_ratio(self):
        return self._max_br

    @property
    def scaling_factor(self):
        return 1 / self.maximum_branching_ratio


class FairShipScalarModel(Model):
    '''
    Specialization of `scalar_portal.Model` with extra FairShip-specific
    features.

    TODO: add dummy particles and channels.
    '''
    def __init__(self, scalar_id=9900025):
        super(FairShipScalarModel, self).__init__(scalar_id)

    def compute_branching_ratios(self, mass, coupling, ignore_invalid=False):
        '''
        Compute the decay and (rescaled) production branching ratios of the
        scalar particle, and return a `BranchingRatiosResult` object containing
        the result.
        '''
        prod_channels  = self.production.get_active_processes()
        decay_channels = self.decays.get_active_processes()
        prod_br  = RescaledProductionBranchingRatios(
            prod_channels , mass, coupling, ignore_invalid, scalar_id=self.scalar_pdg_id)
        decay_br = DecayBranchingRatios(
            decay_channels, mass, coupling, ignore_invalid, scalar_id=self.scalar_pdg_id)
        res = BranchingRatiosResult(prod_br, decay_br)
        return res
