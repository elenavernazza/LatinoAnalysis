# formulas to be added, in python language
# call to branches have to be in the form event.branchName
# if you want to use logical operators, the have to be the python ones (i.e "and" not " and ")

#formulas = {}

formulas_source = {}
formulas_source['SFweight1l'] = 'event.puW*\
                          event.effTrigW1l*\
                          event.std_vector_lepton_recoW[0]*\
                          event.std_vector_lepton_etaW[0]*event.std_vector_lepton_ptW[0]*\
                          event.veto_EMTFBug \
                          if hasattr(event, \'std_vector_lepton_recoW\') else 1.'

muWP='cut_Tight80x'
eleWP='cut_WP_Tight80X'
formulas_source['LepSF1l__ele_'+eleWP+'__mu_'+muWP] = 'event.std_vector_electron_idisoW_'+eleWP+'[0]*\
                                                  event.std_vector_muon_idisoW_'+muWP+'[0]\
                                                  if hasattr(event, \'std_vector_electron_idisoW_'+eleWP+'\') and hasattr(event, \'std_vector_muon_idisoW_'+muWP+'\') else 1.'
formulas_source['LepCut1l__ele_'+eleWP+'__mu_'+muWP] = '((event.std_vector_electron_isTightLepton_'+eleWP+'[0]>0.5 or event.std_vector_muon_isTightLepton_'+muWP+'[0]>0.5)) \
                                                   if hasattr(event, \'std_vector_electron_isTightLepton_'+eleWP+'\') and hasattr(event, \'std_vector_muon_isTightLepton_'+muWP+'\') else 1.'

formulas_source['GenLepMatch1l'] = 'event.std_vector_lepton_genmatched[0] \
                             if hasattr(event, \'std_vector_lepton_genmatched\') else 1. '

METFilter_Common = '(event.std_vector_trigger_special[0]*\
                     event.std_vector_trigger_special[1]*\
                     event.std_vector_trigger_special[2]*\
                     event.std_vector_trigger_special[3]*\
                     event.std_vector_trigger_special[5]\
                   )'

METFilter_MCver  =  '(event.std_vector_trigger_special[8]==-2.)'
METFilter_MCOld  =  '(event.std_vector_trigger_special[6]*event.std_vector_trigger_special[7])'
METFilter_MCNew  =  '(event.std_vector_trigger_special[8]*event.std_vector_trigger_special[9])'
METFilter_MC     =  METFilter_Common + '*' + '(('+METFilter_MCver+'*'+METFilter_MCOld+') or ((not '+METFilter_MCver+')*'+METFilter_MCNew+'))' 

formulas_source['METFilter_MC'] = METFilter_MC

formulas['VBSjjlnu_combined_weight'] = \
  '(' + formulas_source['GenLepMatch1l'] + ')*' + \
  '(' + formulas_source['METFilter_MC']  + ')*' + \
  '(' + formulas_source['SFweight1l']    + ')*' + \
  '(bPogSF_CMVAL)*' + \
  '(' + formulas_source['LepSF1l__ele_'+eleWP+'__mu_'+muWP] + ')*' + \
  '(' + formulas_source['LepCut1l__ele_'+eleWP+'__mu_'+muWP] + ')*'

