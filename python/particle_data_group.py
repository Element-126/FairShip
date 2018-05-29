"""
Written by Alex Seleznov (alex.seleznov@gmail.com)
Added here by Jean-Loup Tastet (jeanloup@nbi.ku.dk) without modifications.
"""

import ROOT

PDGdata = ROOT.TDatabasePDG.Instance()

# SOME USEFUL CONSTANTS
PI = 3.14159265358979323846
VEVhiggs = 246 #(GeV)
GFermi = 1.166379e-05 #(GeV^-2)
alphaStrong = 0.1182#(12)

def PDGname(particle): # PDGname(String particle)
    """ Trim particle name(if needed) for use with the PDG database """
    if 'down' in particle: return 'd'
    elif 'up' in particle: return 'u'
    elif 'strange' in particle: return 's'
    elif 'charm' in particle: return 'c'
    elif 'bottom' in particle: return 'b'
    elif 'beauty' in particle: return 'b'
    elif 'top' in particle: return 't'
    if '1' in particle: particle = particle.replace('1',"'") # why ???!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (not (('-' in particle) or ('+' in particle) or ('0' in particle) or ('nu_' in particle) or ('eta' in particle))
        and (particle not in ['d','u','s','c','b','t','g','d_bar','u_bar','s_bar','c_bar','b_bar','t_bar'])):
        particle += '+'
    return particle

def mass(particle): # PDGname(String particle)
    """ Read particle mass from PDG database """
    particle = PDGname(particle)
    return PDGdata.GetParticle(particle).Mass()

def lifetime(particle): # PDGname(String particle)
    """ Read particle lifetime from PDG database """
    particle = PDGname(particle)
    return PDGdata.GetParticle(particle).Lifetime()

class CKMmatrix():
    """ CKM matrix, from http://pdg.lbl.gov/2013/reviews/rpp2013-rev-ckm-matrix.pdf """
    def __init__(self):
        self.Vud = 0.9742
        self.Vus = 0.2252
        self.Vub = 4.2e-03
        self.Vcd = 0.23
        self.Vcs = 1.
        self.Vcb = 4.1e-02
        self.Vtd = 8.4e-03
        self.Vts = 4.3e-02
        self.Vtb = 0.89
CKM = CKMmatrix()
