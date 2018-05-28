"""
# ==================================================================
#   Python 2 module
#
#   This module provides methods to compute the lifetime and
#   branching ratios of new-physics particles, given their masses
#   and couplings as input parameters.
#
#   Created hnl.py: 30/11/2014 Elena Graverini (elena.graverini@cern.ch)
#   Modified to sparticle.py: 14.02.2018 Alex Seleznov (alex.seleznov@gmail.com)
#
#   Sample of usage:
#     ipython -i sparticle.py
#     In [1]: b = HNL(1.,[1.e-8, 2.e-8, 1.e-9],True)
#     HNLbranchings instance initialized with couplings:
#          U2e   = 1e-08
#	   U2mu  = 2e-08
#	   U2tau = 1e-09
#     and mass:
#	   m = 1.0 GeV
#     In [2]: b.computeNLifetime()
#     Out[2]: 4.777721453160521e-05
#     In [3]: b.findBranchingRatio('N -> pi mu')
#     Out[3]: 0.11826749348890987
#
# ==================================================================
"""

#import sys
#sys.path.append("~/SHIP/ShipSoft/.../new")

import math
import ROOT, shipunit, particle_data_group as pdg # shipunit as shipunits

#
# CUSTOM TOOLS
#
def power(num, p): return math.copysign(pow(abs(num), p), num) # num is a real number

#
# SOME USEFUL(for hnl physics) CONSTANTS
#
class constants():
	def __init__(self):
        	self.decayConstant = {'pi':0.1307*shipunit.GeV,
					'pi0':0.130*shipunit.GeV,
					'rho':0.102*shipunit.GeV, 'rho0':0.102*shipunit.GeV,
					'eta':1.2*0.130*shipunit.GeV,
					'eta1':-0.45*0.130*shipunit.GeV} #GeV^2 ? check units!!
        	self.GFermi = pdg.GFermi / shipunit.GeV**2 # Fermi's constant (GeV^-2)
		self.VEVhiggs = pdg.VEVhiggs * shipunit.GeV
        	self.s2thetaw = 0.23126 # square sine of the Weinberg angle
        	self.toSI_eV = 6.58211928*power(10.,-16) # no units or it messes up!!
        	self.toSI_GeV = self.toSI_eV * power(10.,-9) # no units or it messes up!!
		self.mass_e = pdg.mass('e')*shipunit.GeV
		self.mass_mu = pdg.mass('mu')*shipunit.GeV
		self.mass_tau = pdg.mass('tau')*shipunit.GeV
		self.mass_pi0 = pdg.mass('pi0')*shipunit.GeV
		self.mass_D0 = pdg.mass('D0')*shipunit.GeV
		self.mass_c = pdg.mass('c')*shipunit.GeV
		self.mass_s = pdg.mass('s')*shipunit.GeV
	def alpha_strong(self, x): # quite accurate for 1 < x < 10
		#if x < 1 or 10 < x: return 0.
		return ( 0.9610054796707272 - 0.9862427197771493*x + 0.6771369113160683*x**2 - 0.28886158548276886*x**3 + 0.08121820864497636*x**4 - 0.015582557326431793*x**5 + 0.00207352436895346*x**6 - 0.00019121860059143382*x**7 + 0.000011966730046392183*x**8 - 4.824765324757919e-7*x**9 + 1.1183462961227502e-8*x**10 - 1.1068609414927172e-10*x**11 ) #( 16.645774788988426 - 54.505479945251345*x + 66.77635533674862*x**2 - 29.97395927869271*x**3 - 10.738687566503597*x**4 + 20.71789987640358*x**5 - 12.27014589416471*x**6 + 4.243869291696091*x**7 - 0.9612654432251306*x**8 + 0.14687970811954637*x**9 - 0.014909232015924352*x**10 + 0.0009288709544087798*x**11 - 0.00002567541596100477*x**12 - 6.288931301149232e-7*x**13 + 6.790284314598939e-8*x**14 - 1.483001650401117e-9*x**15 ) / ( 1. + 60.371453676149024*x - 248.36828943918118*x**2 + 428.8060142090867*x**3 - 416.5987126078924*x**4 + 254.47777146665902*x**5 - 103.5049215084586*x**6 + 28.858039894872107*x**7 - 5.548269113704298*x**8 + 0.7185028973006893*x**9 - 0.057347103190308606*x**10 + 0.001848903329423121*x**11 + 0.00012420179699846015*x**12 - 0.000016993171090483*x**13 + 7.7240830563821e-7*x**14 - 1.3368045646953004e-8*x**15 )
	def pipi(self, x): # quite accurate for 0.2780188156668123 < x < 0.9215362374971292
		#if x < 0.2780188156668123 or 0.9215362374971292 < x: return 0.
		return ( -75682.578106848 + 1582850.9028230077*x - 14817302.522803737*x**2 + 81961158.23292133*x**3 - 297704830.6208452*x**4 + 745762615.4045931*x**5 - 1315121176.4816368*x**6 + 1633237497.0807104*x**7 - 1400456995.7400215*x**8 + 790041661.1393455*x**9 - 264034214.31633425*x**10 + 39624763.21078834*x**11 )

CONSTS = constants()

#
# S-particle Id
#
SparticleId = '9900055'

#
# ADDING S-PARTICLE TO ROOT
#
def addSparticleToROOT(pid=9900055 ,m = 0.365, g=1.7678622795749257e-20):
    pdg.PDGdata.AddParticle('S5','S-particle', m, False, g, 0., 'S5', pid) # AddParticle( const char * name, const char * title, Double_t mass, Bool_t stable, Double_t width, Double_t charge, const char * ParticleClass, Int_t PDGcode, Int_t Anti = -1, Int_t TrackingCode = 0 )

#
# PARTICLE INSTANCE CLASS
#
class particleInstance(object): # or newParticle()
	"""
	HNL physics according to the nuMSM
	"""
	def __init__(self, mass, couplings, debug=False):
		"""
		Initialize with mass and couplings of the S-particle

		Inputs:
		mass (GeV)
		couplings (list of [g*])
		"""
		self._mass = mass*shipunit.GeV
		self._couplings = couplings
		self.births = {
			'K0 -> pi0 S': self.birthWidth_fromK0(),
			'B+ -> K+ S': self.birthWidth_B__S_K(),
			'B0 -> K0L S': self.birthWidth_B0__S_K0L(),
			'B0 -> K0S S': self.birthWidth_B0__S_K0S(),
			'D -> S e nu': self.birthWidth_fromD(),
			'B+ -> S pi+': self.birthWidth_B__S_Pi(),
			'B0 -> S pi0': self.birthWidth_B0__S_Pi0()
			#'K+ -> pi+ S': self.birthWidth_fromK(),
			#'K- -> pi- S': self.birthWidth_fromK(),
			#'B+ -> K+ S': self.birthWidth_B__S_K(),
			#'B- -> K- S': self.birthWidth_B__S_K(),
			#'B- -> S pi-': self.birthWidth_B__S_Pi(),
			}
		self.decays = {
			'S -> e+ e-': self.decayWidth_to2Leptons('e+', 'e-'),
			'S -> mu+ mu-': self.decayWidth_to2Leptons('mu+', 'mu-'),
			'S -> tau+ tau-': self.decayWidth_to2Leptons('tau+', 'tau-'),
			'S -> pi0 pi0': self.decayWidth_to2Pions('pi0', 'pi0'),
			'S -> g g': self.decayWidth_to2Gluons(),
			'S -> c c_bar': self.decayWidth_to2CharmQuarks(),
			'S -> s s_bar': self.decayWidth_to2StrangeQuarks()
			#'S -> pi+ pi-': self.decayWidth_to2Pions('pi+', 'pi-'),
			#'S -> d- d': self.decayWidth_to2Quarks('d-', 'd'),
			#'S -> u- u': self.decayWidth_to2Quarks('u-', 'u'),
			#'S -> b- b': self.decayWidth_to2Quarks('b-', 'b'),
			#'S -> t- t': self.decayWidth_to2Quarks('t-', 't')
			}
		if debug:
			print "HNLbranchings instance initialized with coupling:"
            		print "\tg*   = %s"%couplings[0]
            		print "and mass:"
            		print "\tm = %s GeV"%(mass)

	def __repr__(self):
		pass

	@property
	def mass(self):
		return self._mass

	@mass.setter
	def mass(self, new_mass): # new_mass comes in GeV
		self._mass = new_mass * shipunit.GeV # and ends up in SHiP units
		self.births = {
			'K0 -> pi0 S': self.birthWidth_fromK0(),
			'B+ -> K+ S': self.birthWidth_B__S_K(),
			'B0 -> K0L S': self.birthWidth_B0__S_K0L(),
			'B0 -> K0S S': self.birthWidth_B0__S_K0S(),
			'D -> S e nu': self.birthWidth_fromD(),
			'B+ -> S pi+': self.birthWidth_B__S_Pi(),
			'B0 -> S pi0': self.birthWidth_B0__S_Pi0()
			#'K+ -> pi+ S': self.birthWidth_fromK(),
			#'K- -> pi- S': self.birthWidth_fromK(),
			#'B+ -> K+ S': self.birthWidth_B__S_K(),
			#'B- -> K- S': self.birthWidth_B__S_K(),
			#'B- -> S pi-': self.birthWidth_B__S_Pi(),
			}
		self.decays = {
			'S -> e+ e-': self.decayWidth_to2Leptons('e+', 'e-'),
			'S -> mu+ mu-': self.decayWidth_to2Leptons('mu+', 'mu-'),
			'S -> tau+ tau-': self.decayWidth_to2Leptons('tau+', 'tau-'),
			'S -> pi0 pi0': self.decayWidth_to2Pions('pi0', 'pi0'),
			'S -> g g': self.decayWidth_to2Gluons(),
			'S -> c c_bar': self.decayWidth_to2CharmQuarks(),
			'S -> s s_bar': self.decayWidth_to2StrangeQuarks()
			#'S -> pi+ pi-': self.decayWidth_to2Pions('pi+', 'pi-'),
			#'S -> d- d': self.decayWidth_to2Quarks('d-', 'd'),
			#'S -> u- u': self.decayWidth_to2Quarks('u-', 'u'),
			#'S -> b- b': self.decayWidth_to2Quarks('b-', 'b'),
			#'S -> t- t': self.decayWidth_to2Quarks('t-', 't')
			}
	@property
	def couplings(self): return self._couplings

	@couplings.setter
	def couplings(self, new_couplings):
		self._couplings = new_couplings
		self.births = {
			'K0 -> pi0 S': self.birthWidth_fromK0(),
			'B+ -> K+ S': self.birthWidth_B__S_K(),
			'B0 -> K0L S': self.birthWidth_B0__S_K0L(),
			'B0 -> K0S S': self.birthWidth_B0__S_K0S(),
			'D -> S e nu': self.birthWidth_fromD(),
			'B+ -> S pi+': self.birthWidth_B__S_Pi(),
			'B0 -> S pi0': self.birthWidth_B0__S_Pi0()
			#'K+ -> pi+ S': self.birthWidth_fromK(),
			#'K- -> pi- S': self.birthWidth_fromK(),
			#'B+ -> K+ S': self.birthWidth_B__S_K(),
			#'B- -> K- S': self.birthWidth_B__S_K(),
			#'B- -> S pi-': self.birthWidth_B__S_Pi(),
			}
		self.decays = {
			'S -> e+ e-': self.decayWidth_to2Leptons('e+', 'e-'),
			'S -> mu+ mu-': self.decayWidth_to2Leptons('mu+', 'mu-'),
			'S -> tau+ tau-': self.decayWidth_to2Leptons('tau+', 'tau-'),
			'S -> pi0 pi0': self.decayWidth_to2Pions('pi0', 'pi0'),
			'S -> g g': self.decayWidth_to2Gluons(),
			'S -> c c_bar': self.decayWidth_to2CharmQuarks(),
			'S -> s s_bar': self.decayWidth_to2StrangeQuarks()
			#'S -> pi+ pi-': self.decayWidth_to2Pions('pi+', 'pi-'),
			#'S -> d- d': self.decayWidth_to2Quarks('d-', 'd'),
			#'S -> u- u': self.decayWidth_to2Quarks('u-', 'u'),
			#'S -> b- b': self.decayWidth_to2Quarks('b-', 'b'),
			#'S -> t- t': self.decayWidth_to2Quarks('t-', 't')
			}

	@property
	def lifetime(self):
		# if self.mass < 0.2780188156668123: return 0.# None
		if self.mass <= 0.768310546875: return CONSTS.toSI_GeV / self.getDecayWidth1()
		elif self.mass <= 10: return CONSTS.toSI_GeV / self.getDecayWidth2()
		else: return 0.# None
		
	def computeNLifetime(self, system="SI"): # system: choose between default (i.e. SI, result in s) or FairShip (result in ns)
		if system == "FairShip": return 1.e9 * CONSTS.toSI_GeV * self.lifetime
		else: return CONSTS.toSI_GeV * self.lifetime

	# BIRTH
	def pulseS(self, massH1, massH2):
		return power( (massH1**2 - (self.mass + massH2)**2)*(massH1**2 - (self.mass - massH2)**2) , 0.5) / 2 / massH1

	def birthWidth_fromK0(self):
		massH1 = pdg.mass('K0')
		massH2 = pdg.mass('pi0')
		#BR_part = 
		#Pdecay = 
		BRpart_Pdecay = 1
		def my_width(): return BRpart_Pdecay * 2 * self.couplings[0]**2 * self.pulseS(massH1, massH2) / massH1
		return my_width

	def birthWidth_B__S_K(self):
		massH1 = pdg.mass('B+')
		massH2 = pdg.mass('K+')
		#BR_part = 
		#Pdecay = 
		BRpart_Pdecay = 1
		def my_width(): return BRpart_Pdecay * 2 * self.couplings[0]**2 * self.pulseS(massH1, massH2) / massH1
		return my_width

	def birthWidth_B0__S_K0L(self):
		massH1 = pdg.mass('B0')
		massH2 = pdg.mass('K0')
		#BR_part = 
		#Pdecay = 
		# BRpart_Pdecay = 1 / 
		def my_width(): return BRpart_Pdecay * self.couplings[0]**2 * self.pulseS(massH1, massH2) / massH1 # * 2 / 2
		return my_width

	def birthWidth_B0__S_K0S(self):
		massH1 = pdg.mass('B0')
		massH2 = pdg.mass('K0')
		#BR_part = 
		#Pdecay = 
		BRpart_Pdecay = 1 / 1.08
		def my_width(): return BRpart_Pdecay * self.couplings[0]**2 * self.pulseS(massH1, massH2) / massH1 # * 2 / 2
		return my_width

	def birthWidth_fromD(self):
		#BR_part = 
		#Pdecay = 
		BRpart_Pdecay = 1
		def my_width(): return BRpart_Pdecay * self.couplings[0]**2
		return my_width

	def birthWidth_B__S_Pi(self):
		massH1 = pdg.mass('B+')
		massH2 = pdg.mass('pi+')
		#BR_part = 
		#Pdecay = 
		BRpart_Pdecay = 1
		def my_width(): return BRpart_Pdecay * 2 * self.couplings[0]**2 * self.pulseS(massH1, massH2) / massH1
		return my_width

	def birthWidth_B0__S_Pi0(self):
		massH1 = pdg.mass('B0')
		massH2 = pdg.mass('pi0')
		#BR_part = 
		#Pdecay = 
		BRpart_Pdecay = 1 / 1.08
		def my_width(): return BRpart_Pdecay * 2 * self.couplings[0]**2 * self.pulseS(massH1, massH2) / massH1
		return my_width

	def getTotalBirthWidth(self):
		return sum([self.births[channel]() for channel in self.births.values()])

	def getBirthBranching(self, birthChannel): # birthChannel is a string describing the decay, in the form 'stuff0 -> stuff1 ... stuffN + S'
		totalWidth = self.getTotalBirthWidth()
		if birthChannel in self.births.keys():
			if totalWidth == 0.: return 0.
			else: return self.births[birthChannel]() / totalWidth
		else:
		    print 'sparticle.particleInstance.getBirthBranching ERROR: unknown channel %s'%birthChannel
		    quit()

	def getSumBirthBranching(self, birthChannels):
		totalWidth = self.getTotalBirthWidth()
		sumWidth = 0.
		for channel in birthChannels:
			if channel in self.births.keys():
				if totalWidth == 0.: return 0.
				else: sumWidth += self.births[channel]()
			else:
			    print 'sparticle.particleInstance.getBirthBranching ERROR: unknown channel %s'%channel
			    quit()
		return sumWidth / totalWidth

	def getMaxSumBirthBranching(self, birthChannels): # birthChannel is a string describing the decay, in the form 'stuff0 -> stuff1 ... stuffN + S'
		branchings = {}
		for channel in birthChannels:
			meson = channel.replace(' -> ',' ').split(' ')[0]
			if not branchings.has_key(meson): branchings[meson] = 0.
			else: branchings[meson] += self.getBirthBranching(channel)
		return max(branchings.values())

	# DECAYING
	def decayWidth_to2Leptons(self, lepton1, lepton2):
		lepton1_mass = pdg.mass(lepton1)
		lepton2_mass = pdg.mass(lepton2)
		def zero_func():
			return 0.
		if self.mass < (lepton1_mass + lepton2_mass)*shipunit.GeV: return zero_func
		def my_width():
			return ( self.couplings[0]**2 * lepton1_mass**2 * self.mass / 8. / pdg.PI / CONSTS.VEVhiggs**2 ) * power(1 - 4 * lepton1_mass**2 / self.mass**2, 1.5)
		return my_width

	def decayWidth_to2Pions(self, pion1, pion2):
		pion1_mass = pdg.mass(pion1)
		pion2_mass = pdg.mass(pion2)
		def zero_func():
			return 0.
		if self.mass < (pion1_mass + pion2_mass)*shipunit.GeV or self.mass < 0.2780188156668123: return zero_func
		def my_width():
			return CONSTS.pipi(self.mass) * self.decayWidth_to2Leptons('mu+', 'mu-')() * power(1 - (2*pion1_mass/self.mass)**2, 0.5)
		return my_width

	def decayWidth_to2CharmQuarks(self):
		def zero_func():
			return 0.
		if self.mass < 2*pdg.mass('D0')*shipunit.GeV: return zero_func
		quark_mass = pdg.mass('c')
		N_f = 3.
		def c(x):
			if pdg.mass('s') < self.mass < pdg.mass('c'): return power(9.*x/2,4./9)*( 1. + 0.895*x + 1.371*x**2 + 1.952*x**3 )
			elif pdg.mass('c') < self.mass < pdg.mass('b'): return power(25.*x/6,12./25)*( 1. + 1.014*x + 1.389*x**2 + 1.091*x**3 )
			elif pdg.mass('b') < self.mass < pdg.mass('t'): return power(23.*x/6,12./23)*( 1. + 1.175*x + 1.501*x**2 + 0.1725*x**3 )
		m_q = quark_mass * c( CONSTS.alpha_strong(self.mass)/pdg.PI ) / c( CONSTS.alpha_strong(quark_mass)/pdg.PI )
		beta = power(1. - (2*m_q/self.mass)**2, 0.5)
		delta_QCD = 5.67*CONSTS.alpha_strong(self.mass)/pdg.PI + (35.94 - 1.36*N_f)*(CONSTS.alpha_strong(self.mass)/pdg.PI)**2 + (164.14 - 25.77*N_f + 0.259*N_f**2)*(CONSTS.alpha_strong(self.mass)/pdg.PI)**3
		delta_t = (1.57 - 4.*math.log(self.mass/pdg.mass('t')) / 3. + 4.*(math.log(m_q/self.mass))**2 / 9.)*(CONSTS.alpha_strong(self.mass)/pdg.PI)**2
		def my_width():
			return power(1. - (2*pdg.mass('D0')/self.mass)**2, 0.5) * ( 3. * self.couplings[0]**2 * self.mass * m_q**2 / 8. / pdg.PI / CONSTS.VEVhiggs**2 ) * beta**3 * (1. + delta_QCD + delta_t)
		return my_width

	def decayWidth_to2StrangeQuarks(self):
		quark_mass = pdg.mass('s')
		def zero_func():
			return 0.
		if self.mass < 2*quark_mass*shipunit.GeV: return zero_func
		N_f = 3.
		def c(x):
			if pdg.mass('s') < self.mass < pdg.mass('c'): return power(9*x/2,4./9)*( 1. + 0.895*x + 1.371*x**2 + 1.952*x**3 )
			elif pdg.mass('c') < self.mass < pdg.mass('b'): return power(25*x/6,12./25)*( 1. + 1.014*x + 1.389*x**2 + 1.091*x**3 )
			elif pdg.mass('b') < self.mass < pdg.mass('t'): return power(23*x/6,12./23)*( 1. + 1.175*x + 1.501*x**2 + 0.1725*x**3 )
		m_q = quark_mass * c( CONSTS.alpha_strong(self.mass)/pdg.PI ) / c( CONSTS.alpha_strong(2)/pdg.PI )
		beta = power(1. - (2*m_q/self.mass)**2, 0.5)
		delta_QCD = 5.67*CONSTS.alpha_strong(self.mass)/pdg.PI + (35.94 - 1.36*N_f)*(CONSTS.alpha_strong(self.mass)/pdg.PI)**2 + (164.14 - 25.77*N_f + 0.259*N_f**2)*(CONSTS.alpha_strong(self.mass)/pdg.PI)**3
		delta_t = (1.57 - 4.*math.log(self.mass/pdg.mass('t')) / 3. + 4*(math.log(m_q/self.mass))**2 / 9.)*(CONSTS.alpha_strong(self.mass)/pdg.PI)**2
		def my_width():
			return ( 3. * self.couplings[0]**2 * self.mass * m_q**2 / 8. / pdg.PI / CONSTS.VEVhiggs**2 ) * beta**3 * (1. + delta_QCD + delta_t)
		return my_width

	def decayWidth_to2Quarks(self, quark1, quark2):
		quark1_mass = pdg.mass(quark1)
		quark2_mass = pdg.mass(quark2)
		N_f = 3.
		def c(x):
			if pdg.mass('s') < self.mass < pdg.mass('c'): return power(9*x/2,4./9)*( 1. + 0.895*x + 1.371*x**2 + 1.952*x**3 )
			elif pdg.mass('c') < self.mass < pdg.mass('b'): return power(25*x/6,12./25)*( 1. + 1.014*x + 1.389*x**2 + 1.091*x**3 )
			elif pdg.mass('b') < self.mass < pdg.mass('t'): return power(23*x/6,12./23)*( 1. + 1.175*x + 1.501*x**2 + 0.1725*x**3 )
		if quark2 == 'c': m_q = quark1_mass * c( CONSTS.alpha_strong(self.mass)/pdg.PI ) / c( CONSTS.alpha_strong(quark1_mass)/pdg.PI )
		elif quark2 == 's': m_q = quark1_mass * c( CONSTS.alpha_strong(self.mass)/pdg.PI ) / c( CONSTS.alpha_strong(2)/pdg.PI )
		beta = power(1. - (2*m_q/self.mass)**2, 0.5)
		delta_QCD = 5.67*CONSTS.alpha_strong(self.mass)/pdg.PI + (35.94 - 1.36*N_f)*(CONSTS.alpha_strong(self.mass)/pdg.PI)**2 + (164.14 - 25.77*N_f + 0.259*N_f**2)*(CONSTS.alpha_strong(self.mass)/pdg.PI)**3
		delta_t = (1.57 - 4.*math.log(self.mass/pdg.mass('t')) / 3. + 4.*(math.log(m_q/self.mass))**2 / 9.)*(CONSTS.alpha_strong(self.mass)/pdg.PI)**2
		def my_width():
			if self.mass < (quark1_mass + quark2_mass)*shipunit.GeV or self.mass < 2.*shipunit.GeV: return 0.
			else: return ( 3. * self.couplings[0]**2 * self.mass * m_q**2 / 8. / pdg.PI / CONSTS.VEVhiggs**2 ) * beta**3 * (1. + delta_QCD + delta_t)
		return my_width

	def decayWidth_to2Gluons(self):
		def zero_func():
			return 0.
		if self.mass < 2.*pdg.mass('s')*shipunit.GeV: return zero_func
		def y(q): return (2. * pdg.mass(q) / self.mass)**2
		def x(q):
			if y(q) > 1.: return math.atan(1. / math.sqrt(y(q)-1))
			else: return complex(pdg.PI, math.log( (1. + math.sqrt(1-y(q)))/(1. - math.sqrt(1-y(q))) )) / 2
		def F_q(q): return -2. * y(q) * (1. + (1-y(q)) * x(q)**2)
		my_const = ( pdg.alphaStrong / 4. / pdg.PI )**2 * (1. + (pdg.mass('t') / CONSTS.VEVhiggs / pdg.PI)**2 / 8.)
		def my_width():
			return abs(sum([F_q(q) for q in ['u','d','s','c','b','t']]))**2 * (self.couplings[0]**2 * self.mass**3 / 8. / pdg.PI / CONSTS.VEVhiggs**2) * ( CONSTS.alpha_strong(self.mass) / 4. / pdg.PI )**2
		return my_width

	def getDecayWidth1(self):
		return sum([self.decays['S -> e+ e-'](), self.decays['S -> mu+ mu-'](), self.decays['S -> pi0 pi0']()]) #, self.decays['S -> pi+ pi-']()]
	def getDecayWidth2(self):
		return sum([self.decays['S -> tau+ tau-'](), self.decays['S -> mu+ mu-'](), self.decays['S -> s s_bar'](), self.decays['S -> c c_bar'](), self.decays['S -> g g']()])

	def getDecayWidth(self):
		# if self.mass < 0.2780188156668123: return 0.# None
		if self.mass <= 0.768310546875: return self.getDecayWidth1()
		elif self.mass <= 10: return self.getDecayWidth2()
		else: return 0.# None
		# return sum([channel() for channel in self.decays.values()])
		# return sum([self.decays[channel]() for channel in decays])

	def getDecayBranching(self, decayChannel): # decayChannel is a string describing the decay, in the form 'S -> stuff1 ... stuffN'
		totalWidth = self.getDecayWidth()
		if decayChannel in self.decays.keys():
			if totalWidth == 0.: return 0.
			else: return self.decays[decayChannel]() / totalWidth
		else:
		    print 'sparticle.particleInstance.getDecayBranching ERROR: unknown decay %s'%decayChannel
		    quit()

#	def allowedDecays(self): # Returns a dictionary of kinematically allowed decay channels
#		m = self.mass
#		channels = []
#		if m <= 0.768310546875:
#			if m <= 2.*CONSTS.mass_e: channels.append('S -> e+ e-')
#			if m <= 2.*CONSTS.mass_mu: channels.append('S -> mu+ mu-')
#			if m < 0.2780188156668123 or m <= 2.*CONSTS.mass_pi0: channels.append('S -> pi0 pi0')
#		elif m <= 10:
#			if m <= 2.*CONSTS.mass_tau: channels.append('S -> tau+ tau-')
#			if m <= 2.*CONSTS.mass_mu: channels.append('S -> mu+ mu-')
#			if m <= 2.*CONSTS.mass_s: channels.append('S -> s s_bar')
#			if m <= 2.*CONSTS.mass_D0: channels.append('S -> c c_bar')
#			if m <= 2.*CONSTS.mass_s: channels.append('S -> g g')
#		return channels

	def allowedChannels(self): # allowedDecays # Returns a dictionary of kinematically allowed decay channels
		m = self.mass
		channels = {}
		if m <= 0.768310546875:
			if m <= 2.*CONSTS.mass_e: channels.update({'S -> e+ e-':'yes'})
			if m <= 2.*CONSTS.mass_mu: channels.update({'S -> mu+ mu-':'yes'})
			if m < 0.2780188156668123 or m <= 2.*CONSTS.mass_pi0: channels.update({'S -> pi0 pi0':'yes'})
		elif m <= 10:
			if m <= 2.*CONSTS.mass_tau: channels.update({'S -> tau+ tau-':'yes'})
			if m <= 2.*CONSTS.mass_mu: channels.update({'S -> mu+ mu-':'yes'})
			if m <= 2.*CONSTS.mass_s: channels.update({'S -> s s_bar':'yes'})
			if m <= 2.*CONSTS.mass_D0: channels.update({'S -> c c_bar':'yes'})
			if m <= 2.*CONSTS.mass_s: channels.update({'S -> g g':'yes'})
		return channels


# probeParticle = particleInstance(0.365, [1.])
# brFair = []
# m = 0.365
# while m < 10:
# 	brFair.append([m, probeParticle.births['B+ -> K+ S']()])
# 	m += 0.005
# 	probeParticle.mass = m
# print brFair
