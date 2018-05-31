import os, sys, warnings, re
import six
import numpy as np
import scipy.interpolate
import ROOT
import shipunit as u

def addHNLtoROOT(pid=9900015 ,m = 1.0, g=3.654203020370371E-21):
    pdg = ROOT.TDatabasePDG.Instance()
    pdg.AddParticle('N2','HNL', m, False, g, 0., 'N2', pid) 

def readFromAscii(filename="branchingratios"):
    FairShip = os.environ['FAIRSHIP'] 
    ascii = open(FairShip+'/shipgen/%s.dat'%filename)
    h={}
    content = ascii.readlines()
    n = 0
    while n<len(content):
       line = content[n]
       if not line.find('TH1F')<0:
          keys = line.split('|')
          n+=1
          limits = content[n].split(',')
          hname = keys[1]
          if len(keys)<5: keys.append(',') 
          h[ hname ] = ROOT.TH1F(hname,keys[2]+';'+keys[3]+';'+keys[4],int(limits[0]),float(limits[1]),float(limits[2]) )
       else:
          keys = line.split(',')
          h[ hname ].SetBinContent(int(keys[0]),float(keys[1]) )
       n+=1
    return h

def getbr(h,histoname,mass,coupling):
    #0 MeV< mass < 3.200 GeV 
    # FIXME: the above assumption is not valid anymore
    # FIXME: assert that bin width == 1 MeV
    mass = int(mass*1000) 
    try:
        br=h[histoname].GetBinContent(mass)
        br=br*coupling
    except:
        br=0
    return br

def getmaxsumbr(h,histograms,mass,couplings,totaltaubr):
    #0 MeV< mass < 3.200 GeV 
    # FIXME: the above assumption is not valid anymore
    maxsumbr=0.0
    sumbrs={}
    brdstauadded=0
    leptons=['e','mu','tau'] 
    for histoname in histograms:
       item = histoname.split('_') 
       lepton = item[len(item)-1]
       meson = item[0]
       try:
          coupling=couplings[leptons.index(lepton)]
       except:
          coupling=couplings[2] 
       if histoname[:3]=='tau': 
          coupling=couplings[2]
       if not sumbrs.has_key(meson):sumbrs[meson]=0
       sumbrs[meson]+=getbr(h,histoname,mass,coupling)
       if meson=="ds" and brdstauadded==0 and totaltaubr>0.:
          sumbrs[meson]+=totaltaubr
          brdstauadded=1       	  
    maxsumbr=max(sumbrs.values())
    return maxsumbr

def getmaxsumbrrpvsusy(h,histograms,mass,couplings):
    #0 MeV< mass < 3.200 GeV 
    maxsumbr=0.0
    sumbrs={}
    for histoname in histograms: 
       item = histoname.split('_') 
       lepton = item[len(item)-1]
       meson = item[0]
       coupling=couplings[1]
       try:
           sumbrs[meson]+=getbr(h,histoname,mass,coupling)
       except:
           sumbrs[meson]=getbr(h,histoname,mass,coupling)
    print sumbrs.values()
    maxsumbr=max(sumbrs.values())
    return maxsumbr

def gettotalbr(h,histograms,mass,couplings,totaltaubr):
    # FIXME: this function does not handle couplings correctly.
    # It should use the channel dictionary instead of parsing the histogram
    # string, which can lead to a misidentification of the coupling.
    totalbr=0.0 
    leptons=['e','mu','tau'] 
    for histoname in histograms: 
       item = histoname.split('_') 
       lepton = item[len(item)-1]
       # FIXME: Workaround to use the correct couplings for tau decays
       if histoname == 'tau_nu_tau_e':
           coupling = couplings[0]
       elif histoname == 'tau_nu_tau_mu':
           coupling = couplings[1]
       else:
           try:
               coupling=couplings[leptons.index(lepton)]
           except:
               coupling=couplings[2]
               if histoname[:3]=='tau': coupling=couplings[2]
       totalbr+=getbr(h,histoname,mass,coupling)
    return totalbr

def gettotalbrrpvsusy(h,histograms,mass,couplings):
    totalbr=0.0 
    for histoname in histograms: 
       item = histoname.split('_') 
       coupling=couplings[1]
       totalbr+=getbr(h,histoname,mass,coupling)
    return totalbr

def checkChannel(channel):
    """
    Checks consistency between decay channel, lepton ID and coupling.
    """
    lepton_ids = [11, 13, 15]
    lepton_str = ['e', 'mu', 'tau']
    success = True
    if 'decay' in channel:
        particles = channel['decay'].split('_')
        parent = particles[0]
        last_children = particles[-1]
        if parent != 'tau' and 'coupling' in channel:
            if 'idlepton' in channel:
                success = success and abs(channel['idlepton']) == lepton_ids[channel['coupling']]
            if last_children in lepton_str:
                success = success and last_children == lepton_str[channel['coupling']]
    if not success:
        warnings.warn("Consistency checks failed for channel " + str(channel))

def setChannels(P8gen,h,channels,mass,couplings,maxsumBR):
    pdg = P8gen.getPythiaInstance().particleData
    sumBR = 0
    for channel in channels:
        if channel['decay']=="ds_production_tau":
            tauhistograms= ['tau_nu_e_bar_e','tau_nu_mu_bar_mu','tau_nu_tau_e','tau_nu_tau_mu','tau_pi-','tau_K-','tau_rho-']
            BrDs2tauSM = 0.0548
            totaltauBR=BrDs2tauSM * gettotalbr(h,tauhistograms,mass,couplings,0.) # FIXME
            P8gen.SetParameters("431:addChannel      1  "+str(totaltauBR/maxsumBR)+"    0      -15       16")
            sumBR+=totaltauBR/maxsumBR
        else:
            checkChannel(channel)
            br = getbr(h,channel['decay'],mass,couplings[channel['coupling']])
            if br>0:
                if channel['id']=='15':
                    tauhistograms= ['tau_nu_e_bar_e','tau_nu_mu_bar_mu','tau_nu_tau_e','tau_nu_tau_mu','tau_pi-','tau_K-','tau_rho-']
                    totaltauBR=gettotalbr(h,tauhistograms,mass,couplings,0.) # FIXME
                    if len(channel)==4:
                        P8gen.SetParameters(channel['id']+":addChannel      1  "+str(br/totaltauBR)+"    1521       9900015      "+str(channel['idhadron']))
                    else:
                        P8gen.SetParameters(channel['id']+":addChannel      1  "+str(br/totaltauBR)+"    1531       9900015      "+str(channel['idlepton'])+" "+str(channel['idhadron']))
                    sumBR+=br/totaltauBR
                elif len(channel)==4:
                    P8gen.SetParameters(channel['id']+":addChannel      1  "+str(br/maxsumBR)+"    0       9900015      "+str(channel['idlepton']))
                    sumBR+=br/maxsumBR
                else:
                    P8gen.SetParameters(channel['id']+":addChannel      1  "+str(br/maxsumBR)+"   22      "+str(channel['idlepton'])+"       9900015   "+str(channel['idhadron']))
                    sumBR+=br/maxsumBR
    if sumBR<1. and sumBR>0.:
        charge = pdg.charge(int(channel['id']))
        if charge>0:
            P8gen.SetParameters(channel['id']+":addChannel      1   "+str(1.-sumBR)+"    0       22      -11")
        elif charge<0:
            P8gen.SetParameters(channel['id']+":addChannel      1   "+str(1.-sumBR)+"    0       22       11")
        else:
            P8gen.SetParameters(channel['id']+":addChannel      1   "+str(1.-sumBR)+"    0       22      22")

def make_particles_stable(P8gen, above_lifetime):
    p8 = P8gen.getPythiaInstance()
    n=1
    while n!=0:
        n = p8.particleData.nextId(n)
        p = p8.particleData.particleDataEntryPtr(n)
        if p.tau0() > above_lifetime:
            command = str(n)+":mayDecay = false"
            p8.readString(command)
            print "Pythia8 configuration: Made %s stable for Pythia, should decay in Geant4"%(p.name())

def parse_histograms(filepath):
    """
    This function parses a file containing histograms of branching ratios.

    It places them in a dictionary indexed by the decay string (e.g. 'd_K0_e'),
    as a pair (mass, branching ratio), where the mass is expressed in GeV.
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()
    # Define regular expressions matching (sub-)headers and data lines
    th1f_exp      = re.compile(r'^TH1F\|.+')
    header_exp    = re.compile(r'^TH1F\|(.+?)\|BR/U2(.+?)\|HNL mass \(GeV\)\|$')
    subheader_exp = re.compile(r'^\s*?(\d+?),\s*(\d+?\.\d+?),\s*(\d+\.\d+)\s*$')
    data_exp      = re.compile(r'^\s*(\d+?)\s*,\s*(\d+\.\d+)\s*$')
    # Locate beginning of each histogram
    header_line_idx = [i for i in range(len(lines)) if th1f_exp.match(lines[i]) is not None]
    # Iterate over histograms
    histograms = {}
    for offset in header_line_idx:
        # Parse header
        mh = header_exp.match(lines[offset])
        if mh is None or len(mh.groups()) != 2:
            raise ValueError("Malformed header encountered: {0}".format(lines[offset]))
        decay_code = mh.group(1)
        # Parse sub-header (min/max mass and number of points)
        ms = subheader_exp.match(lines[offset+1])
        if ms is None or len(ms.groups()) != 3:
            raise ValueError("Malformed sub-header encountered: {0}".format(lines[offset+1]))
        npoints  = int(ms.group(1))
        min_mass = float(ms.group(2))
        max_mass = float(ms.group(3))
        masses = np.linspace(min_mass, max_mass, npoints, endpoint=False)
        branching_ratios = np.zeros(npoints)
        # Now read the data lines (skipping the two header lines)
        for line in lines[offset+2:offset+npoints+1]:
            md = data_exp.match(line)
            if md is None or len(md.groups()) != 2:
                raise ValueError("Malformed data row encountered: {0}".format(line))
            idx = int(md.group(1))
            br  = float(md.group(2))
            branching_ratios[idx] = br
        histograms[decay_code] = (masses, branching_ratios)
    return histograms

def build_histograms(filepath):
    """
    This function reads a file containing branching ratio histograms, and
    returns a dictionary of linear interpolators of the branching ratios,
    indexed by the decay string.
    """
    histogram_data = parse_histograms(filepath)
    histograms = {}
    for (hist_string, (masses, br)) in six.iteritems(histogram_data):
        histograms[hist_string] = scipy.interpolate.interp1d(
            masses, br, kind='linear', bounds_error=False, fill_value=0, assume_sorted=True)
    return histograms

def get_br(histograms, channel, mass, couplings):
    """
    Utility function used to reliably query the branching ratio for a given
    channel at a given mass, taking into account the correct coupling.
    """
    hist = histograms[channel['decay']]
    coupling = couplings[channel['coupling']]
    normalized_br = hist(mass)
    return normalized_br * coupling

def get_branching_ratios(histograms, channels, mass, couplings, scale=False):
    """
    Query the branching ratios for the passed channels, for a given mass and
    set of couplings.

    If the option `scale` is set to true, the branching ratios are rescaled in
    order to make the event generation more efficient.
    """
    branching_ratios = { ch['decay']: get_br(histograms, ch, mass, couplings)
                         for ch in channels }
    if scale:
        branching_ratios = scale_branching_ratios(branching_ratios, channels)
    return branching_ratios

def scale_branching_ratios(branching_ratios, channels):
    """
    Rescale the branching ratios to optimize the event generation.

    This function rescales the passed branching ratios in order to make the
    event generation as efficient as possible when studying very rare processes,
    while enforcing the invariant that any inclusive branching ratio must
    remain lower that unity.

    This is accomplished by computing, for each particle, the total branching
    ratio to processes of interest, and then dividing all branching ratios by
    the highest of those.
    """
    # Total branching ratios for each particle
    total_brs = {}
    for ch in channels:
        if ch['id'] not in total_brs:
            total_brs[ch['id']] = 0
        total_brs[ch['id']] += branching_ratios[ch['decay']]
    # Find the maximum total branching ratio (over all particles)
    max_total_br = max(total_brs.values())
    # Compute and apply the scaling factor
    scaling_factor = 1 / max_total_br
    return { decay: br * scaling_factor
             for (decay, br) in six.iteritems(branching_ratios) }
