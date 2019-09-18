from itertools import chain
from math import cosh, sqrt, cos
from ROOT import TLorentzVector

hh_branches = {
            "F": [ "mjj_b", "mjj_wjet",
                "b_pt_high", "b_pt_low", "b_etaprod",
                "wjet_pt_high", "wjet_pt_low", 
                "b_eta_high", "b_eta_low",
                "wjet_eta_high", "wjet_eta_low",
                "deltaeta_b",  "deltaphi_b", 
                "deltaeta_wjet", "deltaphi_wjet", 
                "deltaphi_lep_b_high", "deltaphi_lep_b_low", 
                "deltaeta_lep_b_high", "deltaeta_lep_b_low", 
                "deltaphi_lep_wjet_high", "deltaphi_lep_wjet_low",
                "deltaeta_lep_wjet_high", "deltaeta_lep_wjet_low",
                "deltaphi_met_b_high", "deltaphi_met_b_low",
                "deltaphi_met_wjet_high", "deltaphi_met_wjet_low",
                "deltaeta_met_b_high", "deltaeta_met_b_low",
                "deltaeta_met_wjet_high", "deltaeta_met_wjet_low",
                "deltaR_lep_b", "deltaR_lep_wjet",                
                "deltaR_b", "deltaR_wjet",
                "Rwjets_high", "Rwjets_low",
                "A_b", "A_wjet",  "Ht"]
            }

def getDefault():
    output = {}
    for br in hh_branches["F"]:
        output[br] = -9999.
   
    return output


def getHHkinematics_b(bjets, wjets, lepton, met, other_jets, output, debug=False):
    # variables extraction for b-jets
    total_b = TLorentzVector(0,0,0,0)
    b_etas = []
    b_phis = []
    b_pts = []
    for i, j in enumerate(bjets):
        total_b+= j
        b_etas.append(j.Eta())
        b_phis.append(j.Phi())
        b_pts.append(j.Pt())
    if debug:
        print "b-jets pts", b_pts
        print "b-jets etas", b_etas
    deltaeta_b = abs(b_etas[0]- b_etas[1])
    mean_eta_b = sum(b_etas) / 2 
    output["b_pt_high"] = b_pts[0]
    output["b_pt_low"] = b_pts[1]
    output["mjj_b"] = total_b.M()
    output["deltaeta_b"] = deltaeta_b
    output["deltaphi_b"] = abs(bjets[0].DeltaPhi(bjets[1]))
    output["deltaR_b"] = bjets[0].DrEtaPhi(bjets[1])
    output["b_etaprod"] = b_etas[0]*b_etas[1]
    output["b_eta_high"] = abs(b_etas[0])
    output["b_eta_low"] = abs(b_etas[1])

    # Delta Phi with lepton
    output["deltaphi_lep_b_high"] = abs(lepton.DeltaPhi(bjets[0]))
    output["deltaphi_lep_b_low"] = abs(lepton.DeltaPhi(bjets[1]))

    # Delta Eta with lepton
    output["deltaeta_lep_b_high"] = abs(lepton.Eta() - b_etas[0])
    output["deltaeta_lep_b_low"]  = abs(lepton.Eta() - b_etas[1])

    #Delta Phi with MET
    output["deltaphi_met_b_high"] = abs(met.DeltaPhi(bjets[0]))
    output["deltaphi_met_b_low"] = abs(met.DeltaPhi(bjets[1]))

    # Delta Eta with MET
    output["deltaeta_met_b_high"] = abs(met.Eta() - b_etas[0])
    output["deltaeta_met_b_low"]  = abs(met.Eta() - b_etas[1])

    # Look for nearest  jet from lepton
    output["deltaR_lep_b"] = min( [ lepton.DrEtaPhi(bjets[0]), lepton.DrEtaPhi(bjets[1])])

    #Asymmetry
    output["A_b"]  = (b_pts[0] - b_pts[1]) / sum(b_pts)

    Ht = 0.
    for oj in other_jets:
        j_pt, j_eta = oj.Pt(), oj.Eta()
        Ht += j_pt
    # Add b-jets and w-jet to Ht
    for jet in chain(bjets, wjets):
        Ht += jet.Pt()
            
    output["Ht"] = Ht

    return output


def getHHkinematics_w(bjets, wjets, lepton, met, other_jets, output, debug=False):
    # variables extraction for w-jets
    total_wjet = TLorentzVector(0,0,0,0)
    wjet_etas = []
    wjet_phis = []
    wjet_pts = []
    for i, j in enumerate(wjets):
        total_wjet += j
        wjet_etas.append(j.Eta())
        wjet_phis.append(j.Phi())
        wjet_pts.append(j.Pt())
    if debug:
        print "Wjet pts", wjet_pts
        print "Wjet etas", wjet_etas
    output["wjet_pt_high"] = wjet_pts[0]
    output["wjet_pt_low"] = wjet_pts[1]
    output["mjj_wjet"] = total_wjet.M()
    output["deltaphi_wjet"] =  abs(wjets[0].DeltaPhi(wjets[1]))
    output["deltaeta_wjet"] = abs(wjet_etas[0] - wjet_etas[1])
    output["deltaR_wjet"] = wjets[0].DrEtaPhi(wjets[1])
    output["wjet_eta_high"] = abs(wjet_etas[0])
    output["wjet_eta_low"] = abs(wjet_etas[1])


    # Delta Phi with lepton
    output["deltaphi_lep_wjet_high"] = abs(lepton.DeltaPhi(wjets[0]))
    output["deltaphi_lep_wjet_low"] = abs(lepton.DeltaPhi(wjets[1]))

    # Delta Eta with lepton
    output["deltaeta_lep_wjet_high"] = abs(lepton.Eta() - wjet_etas[0])
    output["deltaeta_lep_wjet_low"] = abs(lepton.Eta() - wjet_etas[1])
       
    #Delta Phi with MET
    output["deltaphi_met_wjet_high"] = abs(met.DeltaPhi(wjets[0]))
    output["deltaphi_met_wjet_low"] = abs(met.DeltaPhi(wjets[1]))
 
    # Delta Eta with MET
    output["deltaeta_met_wjet_high"] = abs(met.Eta() - wjet_etas[0])
    output["deltaeta_met_wjet_low"] = abs(met.Eta() - wjet_etas[1])

    # Look for nearest  jet from lepton
    output["deltaR_lep_wjet"] = min( [ lepton.DrEtaPhi(wjets[0]), lepton.DrEtaPhi(wjets[1])])

    #Asymmetry
    output["A_wjet"] = (wjet_pts[0] - wjet_pts[1]) / sum(wjet_pts)

    Ht = 0.
    for oj in other_jets:
        j_pt, j_eta = oj.Pt(), oj.Eta()
        Ht += j_pt
    # Add b-jets and w-jet to Ht
    for jet in chain(bjets, wjets):
        Ht += jet.Pt()
            
    output["Ht"] = Ht

    return output
