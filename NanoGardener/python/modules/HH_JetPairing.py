import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import re

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from LatinoAnalysis.NanoGardener.modules.PairingUtils import *

nearest_massWZ = lambda jets: nearest_mass_pair(jets,85.7863)
nearest_massWZ.__name__ = "nearest_massWZ"

# The dictionary define the name of the tagging strategy and functions to 
# use. The order of the list defines the order of the tagging of VBS and V jets
pairing_strategies = {
    "nearest_massW"      : ["W", nearest_massWZ],
    "max_pt_pair"        : ["W", max_pt_pair],
    "min_deltaeta_pair"  : ["W", min_deltaeta_pair],

}

bTaggingWPs = {
    "deepCSV" : {
        "L" : 0.1522,
        "M" : 0.4941,
        "T" : 0.8001
    }
}


class HH_JetPairing(Module):

    def __init__(self, minpt=20, mode="ALL", bWP="M", debug = False):
        '''
        This modules performs the Jet pairing for HH semileptonic analysis. 
        '''
        self.minpt = minpt
        self.mode = mode
        self.bWP = bWP
        self.debug = debug


    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        #self.initReaders(inputTree)
        self.out = wrappedOutputTree

        # New Branches
        self.out.branch("H_jets", "I", n=2)
        for key in pairing_strategies.keys():
            self.out.branch("W_jets_"+key, "I", n=2)

        

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        # Read branches that may be created by previous step in the chain
        # It's important to read them like this in case they 
        # are created by the step before in a PostProcessor chain. 
#        self.nFatJet = event.nCleanFatJet
        self.rawJet_coll    = Collection(event, 'Jet')
        self.Jet_coll       = Collection(event, 'CleanJet')
#        self.JetNotFat_coll = Collection(event,'CleanJetNotFat')

        # do this check at every event, as other modules might have read further branches
        # if event._tree._ttreereaderversion > self._ttreereaderversion: 
        #     self.initReaders(event._tree)

        
        good_jets, good_jets_ids, b_scores = self.get_jets_vectors(self.minpt)

        bjets = [(i, bscore) for i, bscore in enumerate(b_scores)
                    if bscore >= bTaggingWPs['deepCSV'][self.bWP]]

        good_jets_b_ord = list(sorted(bjets, key=itemgetter(1), reverse=True))

        if self.debug:  print "btag jets: ", good_jets_b_ord

        hpair = [-1,-1]

        if len(bjets) == 0: 
            # Cut the event
            return False
  
        elif len(bjets) == 1 :
            # IN good jets index
            hpair[0] = good_jets_b_ord[0][0]
            
        elif len(bjets) >= 2:
            hpair[0] = good_jets_b_ord[0][0]
            hpair[1] = good_jets_b_ord[1][0]

        # get the remaiming jets with index in good_jets collection
        remain_jets = [(i,j) for i,j in enumerate(good_jets) if i not in hpair]
        
        if len(remain_jets) >= 2: 
            if self.mode=="ALL":
                for key in pairing_strategies:
                    tag, algo = pairing_strategies[key]
                    if self.debug: print "Association: ", tag, algo.__name__,
                    W_jets = [remain_jets[k][0] for k in algo([rj[1] for rj in remain_jets])]

                    # Go back to CleanJet index
                    W_cleanjets = [good_jets_ids[ij] for ij in W_jets]
                    self.out.fillBranch("{}_jets_{}".format(tag,key), W_cleanjets)
            else:
                if self.mode in pairing_strategies:
                    tag, algo = pairing_strategies[self.mode]
                    if self.debug: print "Association: ", tag, algo.__name__,
                    W_jets = [remain_jets[k][0] for k in algo([rj[1] for rj in remain_jets])]
                    # Go back to CleanJet index
                    W_cleanjets = [good_jets_ids[ij] for ij in W_jets]
                    self.out.fillBranch("{}_jets_{}".format(tag,self.mode), W_cleanjets)
                else:
                    print("ERROR! Selected pairing mode not found!!")
                    return False
            
            
        else:
            # Not enought jets left for pairing
            return False

        # Now going back to CleanJet indexes 
        H_cleanjets = [good_jets_ids[ij] for ij in hpair]
        self.out.fillBranch("H_jets", H_cleanjets)
        
                
        """return True (go to next module) or False (fail, go to next event)"""
        return True


    def get_jets_vectors(self, ptmin):
        '''
        Returns a list of 4-momenta for jets looking only at jets
        that are cleaned from FatJets.
        A list of indexes in the collection of CleanJet is returned as a reference. 

        Inserted here an eta interval cut to avoid using the jets in a specified region for 
        tagging. 
        '''
        jets = []
        coll_ids = []
        b_scores = []
        for jetindex in range(len(self.Jet_coll)):
            # index in the original Jet collection
            rawjetid = self.Jet_coll[jetindex].jetIdx
            pt, eta, phi, mass, bvalue = self.Jet_coll[jetindex].pt, \
                        self.Jet_coll[jetindex].eta,\
                        self.Jet_coll[jetindex].phi, \
                        self.rawJet_coll[rawjetid].mass, \
                        self.rawJet_coll[rawjetid].btagDeepB

            if pt < ptmin or pt<0:
                break


            if abs(eta) < 10 :
                p = pt * cosh(eta)
                en = sqrt(p**2 + mass**2)
                vec = TLorentzVector()
                vec.SetPtEtaPhiE(pt, eta, phi, en)
                # check if different from the previous one
                if self.debug:
                    print "Jet index: ", jetindex, "> pt:", pt ," eta:", eta, " phi:", phi, " mass:", mass
                jets.append(vec)
                coll_ids.append(jetindex)
                b_scores.append(bvalue)
        return jets, coll_ids, b_scores


