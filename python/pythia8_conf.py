import ROOT, os, sys, yaml
import shipunit as u
import hnl,rpvsusy
from pythia8_conf_utils import *
from method_logger import MethodLogger
import readDecayTable

def configurerpvsusy(P8gen, mass, couplings, sfermionmass, benchmark, inclusive, deepCopy=False, debug=True):
    # configure pythia8 for Ship usage
    if debug:
        pythia_log=open('pythia8_conf.txt','w')
        P8gen = MethodLogger(P8gen)
    h=readFromAscii("branchingratiosrpvsusybench%d"%benchmark)
    P8gen.UseRandom3() 
    P8gen.SetMom(400)  # beam momentum in GeV 
    if deepCopy: P8gen.UseDeepCopy()
    pdg = ROOT.TDatabasePDG.Instance()
    # let strange particle decay in Geant4
    p8 = P8gen.getPythiaInstance()
    n=1
    while n!=0:
     n = p8.particleData.nextId(n)
     p = p8.particleData.particleDataEntryPtr(n)
     if p.tau0()>1: 
      command = str(n)+":mayDecay = false"
      p8.readString(command)
      print "Pythia8 configuration: Made %s stable for Pythia, should decay in Geant4"%(p.name())
    if inclusive=="True":
        P8gen.SetParameters("SoftQCD:inelastic = on")
        P8gen.SetParameters("PhotonCollision:gmgm2mumu = on")
        P8gen.SetParameters("PromptPhoton:all = on")
        P8gen.SetParameters("WeakBosonExchange:all = on")

    # generate RPV neutralino from inclusive charm hadrons
    if inclusive=="c":
        P8gen.SetParameters("HardQCD::hardccbar  = on")
        # add RPVSUSY
        rpvsusy_instance = rpvsusy.RPVSUSY(mass, couplings, sfermionmass, benchmark, debug=True)
        ctau = rpvsusy_instance.computeNLifetime(system="FairShip") * u.c_light * u.cm
        print "RPVSUSY ctau ",ctau
        P8gen.SetParameters("9900015:new = N2 N2 2 0 0 "+str(mass)+" 0.0 0.0 0.0 "+str(ctau/u.mm)+"  0   1   0   1   0") 
        P8gen.SetParameters("9900015:isResonance = false")
        P8gen.SetParameters("Next:numberCount    =  0")
        # Configuring decay modes...
        rpvsusy_instance.AddChannelsToPythia(P8gen)

        # Finish HNL setup...
        P8gen.SetParameters("9900015:mayDecay = on")
        P8gen.SetHNLId(9900015)
        # also add to PDG
        gamma = u.hbarc / float(ctau) #197.3269631e-16 / float(ctau) # hbar*c = 197 MeV*fm = 197e-16 GeV*cm
        addHNLtoROOT(pid=9900015,m=mass,g=gamma)
        # 12 14 16 neutrinos replace with N2
        charmhistograms = ['d_mu','ds_mu']
        # no tau decay here to consider
        totaltauBR      = 0.0 
        maxsumBR        = getmaxsumbrrpvsusy(h,charmhistograms,mass,couplings)
        if maxsumBR==0.:
           print "No phase space for RPV neutralino from c at this mass:",mass,". Quitting."
           sys.exit()
        totalBR         = gettotalbrrpvsusy(h,charmhistograms,mass,couplings)


        #overwrite D_s+ decays
        P8gen.SetParameters("431:new  D_s+  D_s-    1   3   0    1.96849"\
                            "    0.00000    0.00000    0.00000  1.49900e-01   0   1   0   1   0")
        sumBR=0.
        if getbr(h,'ds_mu',mass,couplings[1])>0.:
           P8gen.SetParameters("431:addChannel      1  "+str(getbr(h,'ds_mu',mass,couplings[1])/maxsumBR)+\
                               "    0      -13       9900015")
           sumBR+=float(getbr(h,'ds_mu',mass,couplings[1])/maxsumBR) 
        if sumBR<1. and sumBR>0.:
            P8gen.SetParameters("431:addChannel      1   "+str(1.-sumBR)+"    0       22      -11")

        #overwrite D+ decays
        P8gen.SetParameters("411:new  D+ D-    1   3   0    1.86962"\
                            "    0.00000    0.00000    0.00000  3.11800e-01   0   1   0   1   0")
        sumBR=0.
        if getbr(h,'d_mu',mass,couplings[1])>0.:
           P8gen.SetParameters("411:addChannel      1  "+str(getbr(h,'d_mu',mass,couplings[1])/maxsumBR)+\
                               "    0      -13       9900015")
           sumBR+=float(getbr(h,'d_mu',mass,couplings[1])/maxsumBR) 
        if sumBR<1. and sumBR>0.:
           P8gen.SetParameters("411:addChannel      1   "+str(1.-sumBR)+"    0       22      -11")

        P8gen.List(9900015)

    if inclusive=="b":
        P8gen.SetParameters("HardQCD::hardbbbar  = on")
        # add RPVSUSY
        rpvsusy_instance = rpvsusy.RPVSUSY(mass, couplings, sfermionmass, benchmark, debug=True)
        ctau = rpvsusy_instance.computeNLifetime(system="FairShip") * u.c_light * u.cm
        P8gen.SetParameters("9900015:new = N2 N2 2 0 0 "+str(mass)+" 0.0 0.0 0.0 "+str(ctau/u.mm)+"  0   1   0   1   0") 
        P8gen.SetParameters("9900015:isResonance = false")
        # Configuring decay modes...
        rpvsusy_instance.AddChannelsToPythia(P8gen)
        # Finish HNL setup...
        P8gen.SetParameters("9900015:mayDecay = on")
        P8gen.SetHNLId(9900015)
        # also add to PDG
        gamma = u.hbarc / float(ctau) #197.3269631e-16 / float(ctau) # hbar*c = 197 MeV*fm = 197e-16 GeV*cm
        addHNLtoROOT(pid=9900015,m=mass,g=gamma)
        # 12 14 16 neutrinos replace with N2
        beautyhistograms = ['b_mu','b_tau','b0_nu_mu','b0_nu_tau']
        maxsumBR=getmaxsumbrrpvsusy(h,beautyhistograms,mass,couplings)
        if maxsumBR==0.:
           print "No phase space for HNL from b at this mass:",mass,". Quitting."
           sys.exit()
        totalBR=gettotalbrrpvsusy(h,beautyhistograms,mass,couplings)

        #overwrite B+ decays
        P8gen.SetParameters("521:new  B+               B-    1   3   0    5.27925"\
                            "0.00000    0.00000    0.00000  4.91100e-01   0   1   0   1   0")
        sumBR=0.
        if getbr(h,'b_tau',mass,couplings[1])>0.:
           P8gen.SetParameters("521:addChannel      1  "+\
                               str(getbr(h,'b_tau',mass,couplings[1])/maxsumBR)+"    0       9900015      -15")
           sumBR+=float(getbr(h,'b_tau',mass,couplings[1])/maxsumBR) 
        if sumBR<1. and sumBR>0.:
           P8gen.SetParameters("521:addChannel      1   "+\
                               str(1.-sumBR)+"    0       22      22")

        #overwrite B0 decays
        P8gen.SetParameters("511:new  B0  Bbar0    1   0   0    5.27958"\
                            "    0.00000    0.00000    0.00000  4.58700e-01   0   1   0   1   0")
        sumBR=0.
        if getbr(h,'b0_nu_tau',mass,couplings[1])>0.:
           P8gen.SetParameters("511:addChannel      1  "+\
                               str(getbr(h,'b0_nu_tau',mass,couplings[1])/maxsumBR)+\
                               "   22       9900015      16")
        if sumBR<1. and sumBR>0.:
           P8gen.SetParameters("511:addChannel      1   "+\
                               str(1.-sumBR)+"    0       22      22")

        P8gen.List(9900015)

    if debug: pythia_log.close()


def configure(P8gen, mass, production_couplings, decay_couplings, inclusive,
              deepCopy=False, debug=True):
    """
    This function configures a HNLPythia8Generator instance for SHiP usage.
    """

    # Wrap the Pythia8 object into a class logging all of its method calls
    if debug:
        pythia_log=open('pythia8_conf.txt','w')
        P8gen = MethodLogger(P8gen, sink=pythia_log)

    fairship_root = os.environ['FAIRSHIP'] 
    histograms = build_histograms('{0}/shipgen/branchingratios.dat'.format(fairship_root))
    P8gen.UseRandom3() # TRandom1 or TRandom3 ?
    P8gen.SetMom(400)  # beam momentum in GeV 
    if deepCopy: P8gen.UseDeepCopy()
    pdg = ROOT.TDatabasePDG.Instance()
    # let strange particle decay in Geant4
    make_particles_stable(P8gen, above_lifetime=1)

    # Load particle & decay data
    # ==========================

    datafile = '{0}/python/hnl_production.yaml'.format(fairship_root)
    with open(datafile, 'rU') as f:
        data = yaml.load(f)
    all_channels  = data['channels']

    # Inclusive
    # =========

    if inclusive=="True":
        P8gen.SetParameters("SoftQCD:inelastic = on")
        P8gen.SetParameters("PhotonCollision:gmgm2mumu = on")
        P8gen.SetParameters("PromptPhoton:all = on")
        P8gen.SetParameters("WeakBosonExchange:all = on")

    # Charm decays only
    # =================

    if inclusive=="c":

        P8gen.SetParameters("HardQCD::hardccbar  = on")
        add_hnl(P8gen, mass, decay_couplings)

        # Add new charmed particles
        # -------------------------

        # Select all charmed particles
        c_particles = data['selections']['c']
        add_particles(P8gen, c_particles + [15], data)

        # Add HNL production channels from charmed particles
        # --------------------------------------------------

        # Find charm and tau decays
        c_channels = [ch for ch in all_channels
                      if ch['id'] in c_particles and ch['decay'] != 'sm']
        tau_channels = [ch for ch in all_channels if abs(ch['id']) == 15]
        # Tau production from D_s+ decay
        c_to_tau_channel = [ch for ch in all_channels
                            if ch['decay'] == 'sm'
                            and abs(ch['children'][0]) == 15][0]

        # Generate all decay chains to compute the branching ratio scaling factor
        # Most charm particles directly decay to HNLs
        c_decay_chains = [DecayChain([make_channel(ch, histograms, mass, production_couplings)])
                          for ch in c_channels]
        # The D_s+ can indirectly produce a HNL by first producing a tau+
        c_to_tau_decay = make_sm_channel(c_to_tau_channel)
        tau_decays = [make_channel(ch, histograms, mass, production_couplings)
                      for ch in tau_channels]
        # D_s+ -> tau+ -> N chains
        c_tau_decay_chains = [DecayChain([c_to_tau_decay, tau_decay])
                              for tau_decay in tau_decays]
        all_c_chains = c_decay_chains + c_tau_decay_chains

        # Compute maximum total branching ration (to rescale all BRs)
        max_total_br = compute_max_total_br(all_c_chains)

        if max_total_br <= 0:
           print("No phase space for HNL from c at this mass: {0}. Quitting.".format(mass))
           sys.exit()

        # Add charm decays
        for ch in c_channels:
            add_channel(P8gen, ch, histograms, mass, production_couplings, 1/max_total_br)

        # Add tau production from D_s+
        # We can also rescale Br(Ds -> tau) and Br(tau -> N X...) as long as
        # Br(Ds -> tau -> N X...) remains the same.
        # Here, we set Br(tau -> N) to unity.
        total_dstau_br = sum(ch.branching_ratio() for ch in c_tau_decay_chains)
        total_tau_br = sum(ch.branching_ratio() for ch in tau_decays)
        P8gen.SetParameters("431:addChannel      1  "+str(total_dstau_br/max_total_br)+"    0      -15       16")

        # Add secondary HNL production from tau
        tau_channels = [ch for ch in data['channels'] if ch['id'] == 15]
        for ch in tau_channels:
            add_tau_channel(P8gen, ch, histograms, mass, production_couplings, 1/total_tau_br)

        # Add dummy channels in place of SM processes
        fill_missing_channels(P8gen, max_total_br, all_c_chains)

        # List channels to confirm that Pythia has been properly set up
        P8gen.List(9900015)

    # Beauty decays only
    # ==================

    if inclusive=="b":

        P8gen.SetParameters("HardQCD::hardbbbar  = on")
        add_hnl(P8gen, mass, decay_couplings)

        # Add beauty particles
        b_particles = data['selections']['b']
        add_particles(P8gen, b_particles, data)

        # Find all decay channels from beauty particles
        b_channels = [ch for ch in all_channels if ch['id'] in b_particles]
        b_decays = [make_channel(ch, histograms, mass, production_couplings)
                    for ch in b_channels]

        # Compute scaling factor
        max_total_br = compute_max_total_br(b_decays)

        if max_total_br <= 0:
           print("No phase space for HNL from b at this mass: {0}. Quitting.".format(mass))
           sys.exit()

        # Add beauty decays
        for ch in b_channels:
            add_channel(P8gen, ch, histograms, mass, production_couplings, 1/max_total_br)

        # Add dummy channels in place of SM processes
        fill_missing_channels(P8gen, max_total_br, b_decays)

        P8gen.List(9900015)

    # Bc decays only
    # ==============

    if inclusive=="bc":

        P8gen.SetParameters("HardQCD::hardbbbar  = on")
        add_hnl(P8gen, mass, decay_couplings)

        # Add B_c+/- particles
        bc_particles = data['selections']['bc']
        add_particles(P8gen, bc_particles, data)

        # Find all of their decay channels
        bc_channels = [ch for ch in all_channels if ch['id'] in bc_particles]
        bc_decays = [make_channel(ch, histograms, mass, production_couplings)
                     for ch in bc_channels]

        # Compute scaling factor
        max_total_br = compute_max_total_br(bc_decays)

        if max_total_br <= 0:
           print("No phase space for HNL from bc at this mass: {0}. Quitting.".format(mass))
           sys.exit()

        # Add B_c+/- decays
        for ch in bc_channels:
            add_channel(P8gen, ch, histograms, mass, production_couplings, 1/max_total_br)

        # Add dummy channels in place of SM processes
        fill_missing_channels(P8gen, max_total_br, bc_decays)

        P8gen.List(9900015)

    if debug: pythia_log.close()

def add_hnl(P8gen, mass, decay_couplings):
    "Adds the HNL to Pythia and ROOT"
    hnl_instance = hnl.HNL(mass, decay_couplings, debug=True)
    ctau = hnl_instance.computeNLifetime(system="FairShip") * u.c_light * u.cm
    print("HNL ctau {}".format(ctau))
    P8gen.SetParameters("9900015:new = N2 N2 2 0 0 "+str(mass)+" 0.0 0.0 0.0 "+str(ctau/u.mm)+"  0   1   0   1   0")
    P8gen.SetParameters("9900015:isResonance = false")
    P8gen.SetParameters("Next:numberCount    =  0")
    # Configuring decay modes...
    readDecayTable.addHNLdecayChannels(P8gen, hnl_instance, conffile=os.path.expandvars('$FAIRSHIP/python/DecaySelection.conf'), verbose=False)
    # Finish HNL setup...
    P8gen.SetParameters("9900015:mayDecay = on")
    P8gen.SetHNLId(9900015)
    # also add to PDG
    gamma = u.hbarc / float(ctau) #197.3269631e-16 / float(ctau) # hbar*c = 197 MeV*fm = 197e-16 GeV*cm
    addHNLtoROOT(pid=9900015,m=mass,g=gamma)
