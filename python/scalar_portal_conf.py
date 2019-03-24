# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

try:
    from scalar_portal_integration import FairShipScalarModel
    _scalar_portal_available = True
    _import_exception = None
except ImportError as e:
    _scalar_portal_available = False
    _import_exception = e

import ROOT
from method_logger import MethodLogger
from pythia8_conf_utils import make_particles_stable, exit_if_zero_br


_ship_beam_momentum = 400. # GeV
_ship_min_stable_ctau = 1. # mm

# All (positive) PDG codes for the strange hadrons present in the input file.
# TODO: add IDs once cascade data for strange hadrons is available.
_ship_strange_hadron_ids = []
# All (positive) PDG codes for the charm hadrons present in the input file.
_ship_charm_hadron_ids   = [411, 421, 431, 4122, 4132, 4232, 4332]
# All (positive) PDG codes for the beauty hadrons present in the input file.
_ship_beauty_hadron_ids  = [511, 521, 531, 5122, 5132, 5232, 5332]

_pythia_scalar_id = 9900025
_pythia_number_count = 0

def configure_scalar_portal(P8gen, mass, coupling, production_from,
                            deep_copy=False, debug=True):
    '''
    This function configures `HNLPythia8Generator` for the scalar particle.
    '''

    if not _scalar_portal_available:
        print_instructions()
        assert(_import_exception is not None)
        raise(_import_exception)

    if debug:
        pythia_log=open('pythia8_conf.txt','w')
        P8gen = MethodLogger(P8gen, sink=pythia_log)

    # Generic PYTHIA configuration
    # ============================

    P8gen.UseRandom3()
    P8gen.SetMom(_ship_beam_momentum)
    if deep_copy:
        P8gen.UseDeepCopy()
    # Frequency at which PYTHIA displays progress.
    P8gen.SetParameters('Next:numberCount = {}'.format(_pythia_number_count))
    # Make long-lived particles decay in Geant 4.
    make_particles_stable(P8gen, above_lifetime=_ship_min_stable_ctau)

    # Add missing hadronic resonances
    # -------------------------------

    _add_missing_meson_resonances(P8gen)

    # Configure PYTHIA for the scalar portal
    # ======================================

    m = FairShipScalarModel(scalar_id=_pythia_scalar_id)

    # Production
    # ----------

    if production_from == 'K':
        _disable_existing_channels(P8gen, _ship_strange_hadron_ids)
        m.production.enable('K -> S pi')
    elif production_from == 'B':
        _disable_existing_channels(P8gen, _ship_beauty_hadron_ids)
        m.production.enable('B -> S K?')
    else:
        raise(ValueError('Unsupported selection for scalar production: {}.'.format(production_from)))

    # Decay
    # -----

    if mass <= 1.0: # GeV
        m.decays.enable('LightScalar')
    elif mass >= 2.0: # GeV
        m.decays.enable('HeavyScalar')
    else:
        raise(ValueError('Extrapolation between 1 and 2 GeV not supported yet.'))

    # Finalize setup
    # --------------

    br = m.compute_branching_ratios(mass, coupling)
    exit_if_zero_br(br.production.maximum_total_branching_ratio, production_from, mass,
                    particle='Scalar')
    scaling_factor = br.production.scaling_factor
    print('All production branching ratios rescaled by {:.4}'.format(scaling_factor))
    pythia_config = br.pythia_full_string()
    for par in pythia_config.split('\n'):
        P8gen.SetParameters(par)
    P8gen.SetHNLId(_pythia_scalar_id)
    P8gen.List(_pythia_scalar_id)

    if debug:
        pythia_log.close()

    # Configure ROOT
    # ==============

    pdg = ROOT.TDatabasePDG.Instance()
    br.root_add_particles(pdg)

def print_instructions():
    # FIXME: write more complete instructions.
    instructions = '''\
    Please install the `scalar_portal` submodule and its dependencies.

    Dependencies:
    * numpy         >= 1.11.0
    * future        >= 0.17.1
    * scipy         >= 0.15.1
    * pandas        >= 0.19.0
    * particletools >= 1.0.1
    * rundec        >= 0.5.1
    '''
    print(instructions)

# Add particle definitions for missing hadronic resonances.
# ---------------------------------------------------------

# I could not find a detailed explanation of mMin and mMax, or how to choose them.
# Looking at existing resonances in ParticleData.xml, it seems that:
#         (mMin, mMax) - m0 ~ ±5 × mWidth
# Source for masses and widths: average values from the 2018 Review of Particle Physics

# FairShip seems to be working well without adding these resonances to the ROOT
# database, so for now I just add them to PYTHIA.

_missing_meson_resonances = {
    # K*(1410) (S=1)
    'K*(1410)0': '100313:new = K*(1410)0 K*(1410)bar0 3 0 0 1.421 0.236 0.23 2.60 0.0',
    'K*(1410)+': '100323:new = K*(1410)+ K*(1410)- 3 3 0 1.421 0.236 0.23 2.60 0.0',
    # K*(1680) (S=1)
    'K*(1680)0': '30313:new = K*(1680)0 K*(1680)bar0 3 0 0 1.718 0.322 0.11 3.30 0.0',
    'K*(1680)+': '30323:new = K*(1680)+ K*(1680)- 3 3 0 1.718 0.322 0.11 3.30 0.0',
    # K*_0(700) (S=0)
    'K*_0(700)0': '9000311:new = K*_0(700)0 K*_0(700)bar0 1 0 0 0.824 0.478 0.0 3.2 0.0',
    'K*_0(700)+': '9000321:new = K*_0(700)+ K*_0(700)- 1 3 0 0.824 0.478 0.0 3.2 0.0',
}

def _add_missing_meson_resonances(P8gen):
    for string in _missing_meson_resonances.values():
        P8gen.SetParameters(string)

# Disable all SM decay channels from initial heavy hadrons
# --------------------------------------------------------

def _disable_existing_channels(P8gen, particle_ids):
    for pid in particle_ids:
        P8gen.SetParameters('{}:onMode = off'.format(pid))
