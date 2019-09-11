import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import TLorentzVector
from math import cosh, sqrt
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import LatinoAnalysis.NanoGardener.data.HHjjlnu_vars as hh_vars


class HHjjlnu_kin(Module):
    '''
    This module calculates several HH semileptonic analysis observables. 
    The H_jets and W_jets have been already associated by the HH_JetPairing module. 
    The mode selects the tagging algorithm.   
    '''
    def __init__(self, minptjet = 20, mode="nearest_massW", debug=False):
        self.minptjet = minptjet
        self.H_jets_var = { "H_jets"}
        self.W_jets_var = { "W_jets_"+mode}
        self.debug = debug      

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        #self.initReaders(inputTree)
        self.out = wrappedOutputTree

        # New Branches
        for typ, branches in hh_vars.hh_branches.items():
            for var in branches:
                self.out.branch(var, typ)
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    
    def analyze(self, event):
        # Read branches that may be created by previous step in the chain
        # It's important to read them like this in case they 
        # are created by the step before in a PostProcessor chain. 
        self.rawJet_coll    = Collection(event, 'Jet')
        self.Jet_coll       = Collection(event, 'CleanJet')
  

        # # do this check at every event, as other modules might have read further branches
        # if event._tree._ttreereaderversion > self._ttreereaderversion: 
        #     self.initReaders(event._tree)

        lepton_raw = Object(event, "Lepton", index=0)
        puppiMET    = Object(event, "PuppiMET")
        
        lep = TLorentzVector()
        lep.SetPtEtaPhiE(lepton_raw.pt, lepton_raw.eta,lepton_raw.phi, lepton_raw.pt * cosh(lepton_raw.eta))
        puppimet = TLorentzVector()
        puppimet.SetPtEtaPhiE(puppiMET.pt, 0., puppiMET.phi, puppiMET.pt)
       
        output = {}
        output = hh_vars.getDefault()

        jets, jets_ids = self.get_jets_vectors(self.minptjet)
        bjets = []
        wjets = []
        other_jets = []
        for jet, jetind in zip(jets, jets_ids):
            if jetind in H_jets_index:  
                bjets.append(jet)
            elif jetind in W_jets_var:                      
                wjets.append(jet)
            else:
                other_jets.append(jet)

        if -1 not in self.H_jets_var:
                output = hh_vars.getHHkinematics_b(bjets, lep, puppimet, output, self.debug)

        if -1 not in self.W_jets_var:
                output = hh_vars.getHHkinematics_w(wjets, lep, puppimet, output, self.debug)

        if -1 not in self.H_jets_var or -1 not in self.W_jets_var:
                output = hh_vars.getHHkinematics_0(bjets, wjets, lep, other_jets, output, self.debug)

        
        # Fill the branches
        for var, val in output.items():
            self.out.fillBranch(var, val)        

        """return True (go to next module) or False (fail, go to next event)"""
        return True

    def get_jets_vectors(self, ptmin):
        
        jets = []
        coll_ids = []
        for ijnf in range(len(self.rawJet_coll)):
            jetindex = self.rawJet_coll[ijnf].jetIdx
            # index in the original Jet collection
            rawjetid = self.Jet_coll[jetindex].jetIdx
            pt, eta, phi, mass = self.Jet_coll[jetindex].pt, \
                        self.Jet_coll[jetindex].eta,\
                        self.Jet_coll[jetindex].phi, \
                        self.rawJet_coll[rawjetid].mass

            if pt < ptmin or pt<0: 
                break
            if abs(eta) > 10 : continue
            p = pt * cosh(eta)
            en = sqrt(p**2 + mass**2)
            vec = TLorentzVector()
            vec.SetPtEtaPhiE(pt, eta, phi, en)
            # check if different from the previous one
            if self.debug:
                print "Jet index: ", jetindex, "> pt:", pt ," eta:", eta, " phi:", phi, " mass:", mass

            jets.append(vec)
            coll_ids.append(jetindex)
        return jets, coll_ids


