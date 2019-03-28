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

_scalar_id = 9900025

def _count_monte_carlo_events(nevents, mass, coupling=1, production_from='B', print_events=False):
    f = ROOT.TFile.Open(_beauty_hadron_cascade)
    tree = f.Get('pythia6')
    P8gen = ROOT.HNLPythia8Generator()
    configure_scalar_portal(P8gen, mass, coupling, production_from, debug=True)
    P8gen.SetParameters('ProcessLevel:all = off')
    p8 = P8gen.getPythiaInstance()
    p8.init()
    hadrons = Counter()
    prod_evts = Counter()
    decay_evts = Counter()
    for i, evt in enumerate(tree):
        if i >= nevents:
            break
        p8.event.reset()
        p8.event.append(int(evt.id), 1, 0, 0, evt.px, evt.py, evt.pz, evt.E, evt.M)
        p8.next()
        if print_events:
            p8.event.list()
        hadron_id = p8.event[1].id()
        hadrons[hadron_id] += 1
        for i in range(p8.event.size()):
            pt = p8.event[i]
            assert(pt.id() != -_scalar_id)
            if pt.id() == _scalar_id:
                assert(pt.mother1() == 1)
                sibling_ix = p8.event[1].daughter2()
                sibling_id = p8.event[sibling_ix].id()
                prod_evts[(hadron_id, sibling_id)] += 1
                children_ixs = pt.daughterList()
                children_ids = tuple(p8.event[ix].id() for ix in children_ixs)
                decay_evts[children_ids] += 1
    return hadrons, prod_evts, decay_evts

_header = \
    '          Channel          | Expected BR | Observed BR (error) | Norm. err. | Pass '
_production_subheader = \
    '----------------------------------- Production ------------------------------------'
_decay_subheader = \
    '-------------------------------------- Decay --------------------------------------'
_pass_str = '\033[92mYes\033[0m' # Green 'Yes'
_fail_str = '\033[91mNo\033[0m'  # Red 'No'

def _print_br_comparison(ch, expected_br, observed_br, threshold=5):
    try:
        exp_br = expected_br[ch]
        obs_br, obs_br_err = observed_br[ch]
        norm_err = (obs_br - exp_br) / obs_br_err if obs_br_err > 0 else 0
        ok = abs(norm_err) < threshold
        print(' {: <26}| {:11.6f} | {:7.5f} +/- {:7.5f} | {:10.6f} | {: <4} '.format(
            ch, exp_br, obs_br, obs_br_err, norm_err, _pass_str if ok else _fail_str))
        return ok
    except KeyError:
        return True

def _compare_branching_ratios(nevents, mass, coupling, production_from,
                              production_proc, decay_proc):
    hadron_ct, prod_ct, decay_ct = _count_monte_carlo_events(
        nevents, mass, coupling, production_from=production_from)
    total_scalars = sum(prod_ct.values())
    assert(sum(decay_ct.values()) == total_scalars)
    m = FairShipScalarModel()
    m.production.enable(production_proc)
    m.decays.enable(decay_proc)
    res = m.compute_branching_ratios(mass, coupling)
    prod_br = res.production.branching_ratios
    decay_br = res.decays.branching_ratios
    observed_prod_br = dict()
    observed_decay_br = dict()
    for ch, br in viewitems(prod_br):
        try:
            parent, children = _from_channel_str(ch)
            parent_id = get_pdg_id(parent)
            sibling_id = get_pdg_id(children[1])
            nb_had = hadron_ct[parent_id]
            nb_dec = prod_ct[(parent_id, sibling_id)]
            obs_br = nb_dec / nb_had
            obs_br_err = nb_dec**(1/2) / nb_had + nb_dec / nb_had**(3/2)
            observed_prod_br[ch] = (obs_br, obs_br_err)
        except:
            print('Channel {} not observed in MC output.'.format(ch))
        cc_ch = _to_channel_str(get_name(-parent_id), ['S', get_name(-sibling_id)])
        try:
            cc_nb_had = hadron_ct[-parent_id]
            cc_nb_dec = prod_ct[(-parent_id, -sibling_id)]
            cc_obs_br = cc_nb_dec / cc_nb_had
            cc_obs_br_err = cc_nb_dec**(1/2) / cc_nb_had + cc_nb_dec / cc_nb_had**(3/2)
            # The expected BR is the same for the charge-conjugate channel.
            prod_br[cc_ch] = prod_br[ch]
            observed_prod_br[cc_ch] = (cc_obs_br, cc_obs_br_err)
        except:
            print('Channel {} not observed in MC output.'.format(cc_ch))
    # Not all channels are simulated (e.g. S -> multi-meson), so we must account for that.
    total_expected_br = 0.0
    for ch, br in viewitems(decay_br):
        try:
            parent, children = _from_channel_str(ch)
            children_ids = tuple(get_pdg_id(c) for c in children)
            nb_dec = decay_ct[children_ids]
            obs_br = nb_dec / total_scalars
            obs_br_err = nb_dec**(1/2) / total_scalars
            observed_decay_br[ch] = (obs_br, obs_br_err)
            total_expected_br += br
        except:
            print('Channel {} not observed in MC output.'.format(ch))
    print('Total BR for simulated channels = {}'.format(total_expected_br))
    all_pass = True
    print(_header)
    print(_production_subheader)
    for ch in prod_br.keys():
        all_pass &= _print_br_comparison(ch, prod_br, observed_prod_br)
    print(_decay_subheader)
    for ch in decay_br.keys():
        all_pass &= _print_br_comparison(ch, decay_br, observed_decay_br)
    assert(all_pass)

def test_mc_bmesons_light_scalar():
    _compare_branching_ratios(10000, 1.2, 1e-2, 'B', 'B -> S K?', 'LightScalar')

def test_mc_bmesons_heavy_scalar():
    _compare_branching_ratios(10000, 4.0, 1e-4, 'B', 'B -> S K?', 'HeavyScalar')

if __name__ == '__main__':
    test_mc_bmesons_light_scalar()
    test_mc_bmesons_heavy_scalar()
