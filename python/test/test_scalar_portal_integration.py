# -*- coding: utf-8 -*-

from __future__ import absolute_import, division

from nose.tools import assert_equals, assert_raises
import numpy as np

from scalar_portal.production.two_body_hadronic import TwoBodyHadronic
from scalar_portal import ProductionBranchingRatios
from scalar_portal_integration import (FairShipScalarModel,
                                       RescaledProductionBranchingRatios)

def test_rescaled_production_br():
    channels = [TwoBodyHadronic('B+', 'pi+'), TwoBodyHadronic('B+', 'K*_2(1430)+')]
    mS = np.array([1, 4])
    br  = ProductionBranchingRatios(        channels, mS, 1)
    sbr = RescaledProductionBranchingRatios(channels, mS, 1)
    br_pi  = br.branching_ratios[ 'B+ -> S pi+'        ]
    br_K2  = br.branching_ratios[ 'B+ -> S K*_2(1430)+']
    sbr_pi = sbr.branching_ratios['B+ -> S pi+'        ]
    sbr_K2 = sbr.branching_ratios['B+ -> S K*_2(1430)+']
    assert_equals(list(br_K2 > br_pi), [True , False])
    assert_equals(list(sbr_K2 == 1),   [True , False])
    assert_equals(list(sbr_pi == 1),   [False, True ])
    assert(np.all(np.fmax.reduce(sbr._rescaled_br.values()) == 1))
    sf = sbr.scaling_factor
    epsilon = 1e-14
    assert(np.all(sf * sbr.maximum_branching_ratio - 1) < epsilon)
    assert(abs(sf[0]*br_K2[0] - 1) < epsilon)
    assert(abs(sf[1]*br_pi[1] - 1) < epsilon)

def test_fairship_scalar_model():
    m = FairShipScalarModel()
    m.production.enable('B -> S K' )
    m.production.enable('B -> S K*')
    m.decays.enable('S -> e+ e-'  )
    m.decays.enable('S -> mu+ mu-')
    res = m.compute_branching_ratios(0.5, 1)
    # FIXME: include dummy particles and channels
    assert_equals(res.pythia_full_string(), '''\
9900025:new = S S 1 0 0 0.5 0.0 0.0 0.0 7.22250988672e-05
9900025:isResonance = false
9900025:mayDecay = true
9900025:isVisible = false
521:addChannel = 1 0.830119639659 0 9900025 321
511:addChannel = 1 0.252990892052 0 9900025 30313
511:addChannel = 1 0.400325751322 0 9900025 100313
521:addChannel = 1 0.272630974461 0 9900025 30323
511:addChannel = 1 0.927960927961 0 9900025 313
511:addChannel = 1 0.770318591137 0 9900025 311
521:addChannel = 1 1.0 0 9900025 323
521:addChannel = 1 0.43140367149 0 9900025 100323
9900025:addChannel = 1 3.14194722291e-05 0 -11 11
9900025:addChannel = 1 0.999968580528 0 -13 13'''
    )
