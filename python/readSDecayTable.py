"""
# ==================================================================
#   Python module
#
#   This module provides methods to read a configuration file
#   storing on/off (yes/no) flags for the HNL decay channels
#   and to pass the configured decay table to a Pythia8 generator.
#
#   Created: 30/12/2014 Elena Graverini (elena.graverini@cern.ch)
#   Modified by: Alex Seleznov (alex.seleznov@gmail.com)
#   Integrated: 28/05/2018 Jean-Loup Tastet (jeanloup@nbi.ku.dk)
#       Minimal modifications/fixes were made to ensure
#       compatibility with the existing code base.
#
# ==================================================================
"""
import ROOT, os, csv
from particle_data_group import PDGname

pdg = ROOT.TDatabasePDG.Instance()

def PDGcode(particle):
    """
    Read particle ID from PDG database
    """
    particle = PDGname(particle)
    tPart = pdg.GetParticle(particle)
    if tPart != None:
        return int(tPart.PdgCode())
    else:
        raise ValueError("Particle " + str(particle) + " not found in PDG")


def load(conffile = os.path.expandvars('$FAIRSHIP/python/DecaySelection.conf'), verbose=True):
    f = open(conffile,'r')
    reader = csv.reader(f, delimiter=':')
    configuredDecays = {}
    for row in reader:
        if not row: continue # skip empty lines
	if str(row[0]).strip().startswith('#'):
            continue # skip comment / header lines
        channel = str(row[0]).strip()
        flag = str(row[-1]).partition('#')[0].strip() # skip line comments
        configuredDecays[channel] = flag
    if verbose:
        print 'Activated decay channels (plus charge conjugates): '
        for channel in configuredDecays.keys():
            if configuredDecays[channel] == 'yes':
                print '\t'+channel        
    return configuredDecays 

def addSparticleDecayChannels(P8Gen, sparticle, conffile=os.path.expandvars('$FAIRSHIP/python/DecaySelection_sparticle.conf'), verbose=True):
    """
    Configures the HNL decay table in Pythia8
    Inputs:
    - P8Gen: an instance of ROOT.HNLPythia8Generator()
    - hnl: an instance of hnl.HNL()
    - conffile: a file listing the channels one wishes to activate
    """
    # First fetch the list of kinematically allowed decays
    # allowed = sparticle.allowedDecays()
    allowed = sparticle.allowedChannels()
    # Then fetch the list of desired channels to activate
    wanted = load(conffile=conffile, verbose=verbose)
    # Add decay channels
    for dec in allowed:
        if dec not in wanted:
            print 'addSparticleDecayChannels ERROR: channel not configured!\t', dec
            quit()
        if allowed[dec] == 'yes' and wanted[dec] == 'yes':
            particles = [p for p in dec.replace('->',' ').split()]
            children = particles[1:]
            print(children)
            childrenCodes = [PDGcode(p) for p in children]
            # BR = sparticle.findBranchingRatio(dec)
            BR = sparticle.getDecayBranching(dec)
            codes = ' '.join([str(code) for code in childrenCodes])
            P8Gen.SetParameters("9900055:addChannel =  1 " + str(BR) + " 0 " + codes)
            # print "debug readdecay table",particles,children,BR

if __name__ == '__main__':
    configuredDecays = load()
    print 'Activated decay channels: '
    for channel in configuredDecays.keys():
        if configuredDecays[channel] == 'yes':
            print channel
