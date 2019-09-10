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
pairing_strategies_resolved = {
    "nearest_massW"      : ["W", nearest_massW],
    "max_pt_pair"        : ["W", max_pt_pair],
    "min_deltaeta_pair"  : ["W", min_deltaeta_pair],

}




class HH_JetPairing(Module):

    def __init__(self, minpt=20, mode="ALL", bWP="M" debug = False):
        '''
        This modules performs the Jet pairing for HH semileptonic analysis. 
        It separates events in three categories: boosted and resolved. 

        In the boosted category, only events with 1 FatJet are saved. Events with more FatJets 
        are vetoed. In the remaining jets the VBS pair is selected using the maximum invariant mass. 

        In the resolved category (>= 4 jets) different algorithms can be used to 
        choose the VBS jets and V jets. 

        An eta interval cut can be specified to avoid using the jets in those regions for tagging

        Modes (for resolved category):
        "maxmjj_massWZ" : before VBS jets with max Mjj, than V-jets with mass nearest to (mW+mZ)/2
        "maxmjj_maxPt" : before VBS jets with max Mjj, then V-jets as the pair with max Pt,
        "maxmjj_maxPtsingle"  : before VBS jets with max Mjj, then V-jets as the two jets with highest Pt
        "maxPt_massWZ"  : before VBS jets as pair with max Pt, then V-jets with mass W,Z
        "maxPtsingle_massWZ" : befor VBS jets as the two jets with highest Pt, then V-jet with mass W,Z

        "massWZ_maxmjj": before V jets with mass nearest to W,Z, then VBS jets with MaxMjj
        "massWZ_maxPt" :before V jets with mass nearest to W,Z, then VBS jets with pair with max Pt
       

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
        for key in pairing_strategies_resolved.keys():
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

        if self.debug:  print bjets

        hpair = [-1,-1]
  
        if len(bjets) ==1 :
            hpair = [j[0] for j in list(sorted(bjets, key=itemgetter(1), reverse=True))[:1]]
            hpair.append(-1)

        elif len(bjets) >= 2:
            hpair = [j[0] for j in list(sorted(bjets, key=itemgetter(1), reverse=True))[:2]]        

        self.out.fillBranch("H_jets", hpair)


            # Cache of association algos
            # (N.B. indexes by good_jets collections)
            cache = { }  # algo: ( associated_jets, remaining jets)

            if self.mode=="ALL":
                for key, algos in pairing_strategies_resolved.items():
                    self.perform_jet_association(key, good_jets, good_jets_ids, cache, hpair)
            else:
                if self.mode in pairing_strategies_resolved:
                    self.perform_jet_association(self.mode, good_jets, good_jets_ids, cache, hpair)
                else:
                    print("ERROR! Selected pairing mode not found!!")
                    return False
        else:
            # Cut the event:
            # or it's boosted but with not enough jets, 
            # or it is not boosted and it has less than 4 jets with minpt
            #print("Event removed")
            return False


                
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
        for ijnf in range(len(self.JetNotFat_coll)):
            jetindex = self.JetNotFat_coll[ijnf].jetIdx
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
        return jets, coll_ids, bvalue


    def perform_jet_association(self, mode, good_jets, good_jets_ids, cache, hpair):
        '''
            This function perform the association of the jets with one of the 
            algorithm in pairing_strategies_resolved map.
        '''
        (tag1, algo1), (tag2, algo2) = pairing_strategies_resolved[mode]
        if self.debug: print "Association: ", tag1, algo1.__name__, tag2, algo2.__name__


        if mode == "nearest_massW":
            W_jets = utils.nearest_mass_pair_notH(good_jets, 80.385, hpair)
      

        elif mode == "max_pt_pair":
            W_jets = utils.max_pt_pair_notH(good_jets, hpair)

        elif mode == "min_deltaeta_pair":
            W_jets = utils.min_deltaeta_pairs_notH(jets, hpair)

        # Now going back to CleanJet indexes 
        W_cleanjets = [good_jets_ids[ij] for ij in W_jets]


        self.out.fillBranch("{}_jets_{}".format(tag2,mode), W_cleanjets)

