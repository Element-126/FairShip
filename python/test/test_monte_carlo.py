# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from future.utils import viewitems

from collections import Counter

from scalar_portal_conf import configure_scalar_portal
from scalar_portal_integration import FairShipScalarModel
from scalar_portal.data.particles import *
from scalar_portal.api.channel import _from_channel_str, _to_channel_str

import ROOT

_beauty_hadron_cascade = '/eos/experiment/ship/data/Beauty/Cascade-run0-19-parp16-MSTP82-1-MSEL5-5338Bpot.root'

_scalar_id = 9900015

def _safe_get_pdg_id(name):
    if name == 'S':
        return _scalar_id
    else:
        return get_pdg_id(name)

def _safe_get_name(pdg_id):
    if pdg_id == _scalar_id:
        return 'S'
    else:
        return get_name(pdg_id)

def _get_antiparticle(pdg_id):
    if pdg_id == _scalar_id:
        return _scalar_id
    else:
        return -pdg_id

def _count_monte_carlo_events(nevents, mass, couplings=[1, 0],
                              production_from='B', print_events=False):
    f = ROOT.TFile.Open(_beauty_hadron_cascade)
    tree = f.Get('pythia6')
    P8gen = ROOT.HNLPythia8Generator()
    configure_scalar_portal(P8gen, mass, couplings, couplings, production_from, debug=True)
    P8gen.SetParameters('ProcessLevel:all = off')
    p8 = P8gen.getPythiaInstance()
    p8.init()
    hadrons = Counter()
    prod_evts = Counter()
    decay_evts = Counter()
    # Generate events
    for i, evt in enumerate(tree):
        if i >= nevents:
            break
        p8.event.reset()
        p8.event.append(int(evt.id), 1, 0, 0, evt.px, evt.py, evt.pz, evt.E, evt.M)
        p8.next()
        if print_events:
            p8.event.list()
        # Analyse event
        hadron_id = p8.event[1].id()
        hadrons[hadron_id] += 1
        for i in range(p8.event.size()):
            pt = p8.event[i]
            assert(pt.id() != -_scalar_id)
            if pt.id() == _scalar_id:
                assert(pt.mother1() == 1)
                sibling_ixs = p8.event[1].daughterList()
                sibling_ids = tuple(p8.event[ix].id() for ix in sibling_ixs)
                prod_evts[(hadron_id, sibling_ids)] += 1
                children_ixs = pt.daughterList()
                children_ids = tuple(p8.event[ix].id() for ix in children_ixs)
                decay_evts[children_ids] += 1
                break # Avoids double-counting if there are two scalars
    return hadrons, prod_evts, decay_evts

_header = \
    '           Channel           | Expected BR | Observed BR (error) | Norm. err. | Pass '
_production_subheader = \
    '------------------------------------ Production -------------------------------------'
_decay_subheader = \
    '--------------------------------------- Decay ---------------------------------------'
_pass_str = '\033[92mYes\033[0m' # Green 'Yes'
_fail_str = '\033[91mNo\033[0m'  # Red 'No'

def _print_br_comparison(ch, expected_br, observed_br, err_br, threshold=5):
    try:
        exp_br = expected_br[ch]
        obs_br = observed_br[ch]
        err = err_br[ch]
        norm_err = 0 if (exp_br == 0 and obs_br == 0) else (obs_br - exp_br) / err
        ok = abs(norm_err) < threshold
        norm_err_str = '{:10.6f}'.format(norm_err) if norm_err != 0 else '        --'
        print(' {: <28}| {:11.6f} | {:7.5f} +/- {:7.5f} | {: <10} | {: <4} '.format(
            ch, exp_br, obs_br, err, norm_err_str, _pass_str if ok else _fail_str))
        return ok
    except KeyError:
        return True

def _compare_branching_ratios(nevents, mass, couplings, production_from,
                              production_proc, decay_proc):
    hadron_ct, prod_ct, decay_ct = _count_monte_carlo_events(
        nevents, mass, couplings, production_from=production_from)
    total_scalars = sum(prod_ct.values())
    assert(sum(decay_ct.values()) == total_scalars)
    m = FairShipScalarModel()
    for proc in production_proc:
        m.production.enable(proc)
    for proc in decay_proc:
        m.decays.enable(proc)
    res = m.compute_branching_ratios(mass, {'theta': couplings[0], 'alpha': couplings[1]})
    prod_br = res.production.branching_ratios
    decay_br = res.decays.branching_ratios
    observed_prod_br = dict()
    err_prod_br = dict()
    observed_decay_br = dict()
    err_decay_br = dict()
    for ch in list(prod_br.keys()):
        parent, children = _from_channel_str(ch)
        parent_id = _safe_get_pdg_id(parent)
        children_ids = tuple(_safe_get_pdg_id(c) for c in children)
        cc_ch = _to_channel_str(_safe_get_name(_get_antiparticle(parent_id)),
                                [_safe_get_name(_get_antiparticle(i)) for i in children_ids])
        try:
            nb_had = hadron_ct[parent_id]
            nb_dec = prod_ct[(parent_id, children_ids)]
            observed_prod_br[ch] = nb_dec / nb_had
            exp_nb_dec = prod_br[ch] * nb_had
            err_prod_br[ch] = exp_nb_dec**(1/2) / nb_had + exp_nb_dec / nb_had**(3/2)
        except:
            print('Channel {} not observed in MC output.'.format(ch))
        try:
            # The expected BR is the same for the charge-conjugate channel.
            prod_br[cc_ch] = prod_br[ch]
            cc_nb_had = hadron_ct[_get_antiparticle(parent_id)]
            cc_nb_dec = prod_ct[(_get_antiparticle(parent_id), tuple(map(_get_antiparticle, children_ids)))]
            observed_prod_br[cc_ch] = cc_nb_dec / cc_nb_had
            cc_exp_nb_dec = prod_br[ch] * cc_nb_had
            err_prod_br[cc_ch] = cc_exp_nb_dec**(1/2) / cc_nb_had + cc_exp_nb_dec / cc_nb_had**(3/2)
        except:
            print('Channel {} not observed in MC output.'.format(cc_ch))
    # Not all channels are simulated (e.g. S -> multi-meson), so we must account for that.
    total_expected_br = 0.0
    for ch, br in viewitems(decay_br):
        try:
            parent, children = _from_channel_str(ch)
            children_ids = tuple(get_pdg_id(c) for c in children)
            nb_dec = decay_ct[children_ids]
            observed_decay_br[ch] = nb_dec / total_scalars
            exp_nb_dec = br * total_scalars
            err_decay_br[ch] = exp_nb_dec**(1/2) / total_scalars
            total_expected_br += br
        except:
            print('Channel {} not observed in MC output.'.format(ch))
    print('Total BR for simulated channels = {}'.format(total_expected_br))
    corrected_decay_br = {ch: br / total_expected_br for ch, br in viewitems(decay_br)}
    corrected_err_decay_br = {ch: br / total_expected_br for ch, br in viewitems(err_decay_br)}
    all_pass = True
    print(_header)
    print(_production_subheader)
    for ch in prod_br.keys():
        all_pass &= _print_br_comparison(ch, prod_br, observed_prod_br, err_prod_br)
    print(_decay_subheader)
    for ch in decay_br.keys():
        all_pass &= _print_br_comparison(ch, corrected_decay_br, observed_decay_br, err_decay_br)
    assert(all_pass)

def test_mc_bmesons_light_scalar():
    _compare_branching_ratios(100000, 1.2, [1e-2, 0], 'B', ['B -> S K?'], ['LightScalar'])

def test_mc_bmesons_heavy_scalar():
    _compare_branching_ratios(100000, 4.0, [1e-4, 0], 'B', ['B -> S K?'], ['HeavyScalar'])

def test_mc_bmesons_light_scalar_quartic():
    _compare_branching_ratios(100000, 1.0, [1e-8, 0.1], 'B', ['B -> S S K?', 'B -> S S'], ['LightScalar'])

def test_mc_bmesons_heavy_scalar_quartic():
    _compare_branching_ratios(100000, 2.0, [1e-8, 0.1], 'B', ['B -> S S K?', 'B -> S S'], ['HeavyScalar'])

if __name__ == '__main__':
    test_mc_bmesons_light_scalar()
    test_mc_bmesons_heavy_scalar()
