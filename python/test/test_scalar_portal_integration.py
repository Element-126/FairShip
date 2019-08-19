# -*- coding: utf-8 -*-

from __future__ import absolute_import, division

from nose.tools import assert_equals, assert_raises
import numpy as np

import ROOT

from scalar_portal.production.two_body_hadronic import TwoBodyHadronic
from scalar_portal import ProductionBranchingRatios
from scalar_portal_integration import *

def test_rescaled_production_br():
    channels = [TwoBodyHadronic('B+', 'pi+'), TwoBodyHadronic('B+', 'K*_2(1430)+'),
                TwoBodyHadronic('B0', 'pi0'), TwoBodyHadronic('B0', 'K*_2(1430)0')]
    mS = np.array([1, 4])
    br  = ProductionBranchingRatios(        channels, mS, {'theta': 1})
    sbr = RescaledProductionBranchingRatios(channels, mS, {'theta': 1})
    br_pi1  = br.branching_ratios[ 'B+ -> S pi+'        ]
    br_K21  = br.branching_ratios[ 'B+ -> S K*_2(1430)+']
    sbr_pi1 = sbr.branching_ratios['B+ -> S pi+'        ]
    sbr_K21 = sbr.branching_ratios['B+ -> S K*_2(1430)+']
    br_pi0  = br.branching_ratios[ 'B0 -> S pi0'        ]
    br_K20  = br.branching_ratios[ 'B0 -> S K*_2(1430)0']
    sbr_pi0 = sbr.branching_ratios['B0 -> S pi0'        ]
    sbr_K20 = sbr.branching_ratios['B0 -> S K*_2(1430)0']
    assert_equals(list(br_K21 > br_pi1), [True , False])
    epsilon = 1e-14
    assert(np.all(np.abs((br_pi1 + br_K21) / sbr.maximum_total_branching_ratio - 1) < epsilon))
    assert(np.all(np.abs(sbr_pi1 + sbr_K21 - 1) < epsilon))
    assert(np.all(np.abs(( br_K21/ br_pi1) - ( br_K20/ br_pi0)) <= epsilon * ( br_K20/ br_pi0)))
    assert(np.all(np.abs((sbr_K21/sbr_pi1) - (sbr_K20/sbr_pi0)) <= epsilon * (sbr_K20/sbr_pi0)))
    sf = sbr.scaling_factor
    assert(np.all(np.abs(sf * sbr.maximum_total_branching_ratio - 1)) < epsilon)
    assert(np.all(np.abs(sf * br_pi1 - sbr_pi1) <= epsilon * sbr_pi1))
    assert(np.all(np.abs(sf * br_pi0 - sbr_pi0) <= epsilon * sbr_pi0))
    assert(np.all(np.abs(sf * br_K21 - sbr_K21) <= epsilon * sbr_K21))
    assert(np.all(np.abs(sf * br_pi0 - sbr_pi0) <= epsilon * sbr_pi0))

def test_fairship_scalar_model():
    m = FairShipScalarModel(scalar_id=9900025)
    m.production.enable('B -> S K' )
    m.production.enable('B -> S K*')
    m.decays.enable('S -> e+ e-'  )
    m.decays.enable('S -> mu+ mu-')
    res = m.compute_branching_ratios(0.5, theta=1)
    assert_equals(res.pythia_full_string(), '''\
9900025:new = S void 1 0 0 0.5 0.0 0.0 0.0 7.22250988672e-05
9900025:isResonance = false
9900025:mayDecay = true
9900025:isVisible = false
9900080:new = r0 void 1 0 0 0.0 0.0 0.0 0.0 0.0
9900080:isResonance = false
9900080:mayDecay = false
9900080:isVisible = false
9900081:new = r+ r- 1 3 0 0.0 0.0 0.0 0.0 0.0
9900081:isResonance = false
9900081:mayDecay = false
9900081:isVisible = false
511:onMode = off
521:onMode = off
511:addChannel = 1 0.303974622031 0 9900025 311
511:addChannel = 1 0.366181701418 0 9900025 313
511:addChannel = 1 0.0998324740874 0 9900025 30313
511:addChannel = 1 0.157972130424 0 9900025 100313
521:addChannel = 1 0.327572651899 0 9900025 321
521:addChannel = 1 0.394608965081 0 9900025 323
521:addChannel = 1 0.107582626681 0 9900025 30323
521:addChannel = 1 0.170235756339 0 9900025 100323
511:addChannel = 1 0.0720390720391 0 9900080 9900080
9900025:addChannel = 1 0.999968580528 0 -13 13
9900025:addChannel = 1 3.14194722291e-05 0 -11 11'''
    )

def test_root_integration():
    m = FairShipScalarModel(scalar_id=9900025)
    m.production.enable('B -> S K' )
    m.production.enable('B -> S K*')
    m.decays.enable('S -> e+ e-'  )
    m.decays.enable('S -> mu+ mu-')
    res = m.compute_branching_ratios(0.5, theta=1)
    pdg = ROOT.TDatabasePDG()
    res.root_add_particles(pdg)
    S  = pdg.GetParticle('S' )
    r0 = pdg.GetParticle('r0')
    rp = pdg.GetParticle('r+')
    rm = pdg.GetParticle('r-')
    assert_equals( S.PdgCode(), +9900025)
    assert_equals(r0.PdgCode(), +9900080)
    assert_equals(rp.PdgCode(), +9900081)
    assert_equals(rm.PdgCode(), -9900081)
    assert_equals( S.Charge(),  0.0)
    assert_equals(r0.Charge(),  0.0)
    assert_equals(rp.Charge(), +1.0)
    assert_equals(rm.Charge(), -1.0)
    assert_equals( S.Stable(), False)
    assert_equals(r0.Stable(), True )
    assert_equals(rp.Stable(), True )
    assert_equals(rm.Stable(), True )
    assert_equals(S.Width(), res.total_width)

def test_root_integration_error():
    m = FairShipScalarModel()
    m.production.enable('B -> S K' )
    m.production.enable('B -> S K*')
    m.decays.enable('S -> e+ e-'  )
    m.decays.enable('S -> mu+ mu-')
    res = m.compute_branching_ratios([0.5, 1.5, 3.0], theta=1)
    pdg = ROOT.TDatabasePDG()
    assert_raises(ValueError, lambda: res.root_add_particles(pdg))
