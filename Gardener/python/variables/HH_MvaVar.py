from LatinoAnalysis.Gardener.gardening import TreeCloner


print    ''' 

  _    _ _    _   __  ____      ____      __        
 | |  | | |  | | |  \/  \ \    / /\ \    / /        
 | |__| | |__| | | \  / |\ \  / /  \ \  / /_ _ _ __ 
 |  __  |  __  | | |\/| | \ \/ / /\ \ \/ / _` | '__|
 | |  | | |  | | | |  | |  \  / ____ \  / (_| | |   
 |_|  |_|_|  |_| |_|  |_|   \/_/    \_\/ \__,_|_|   
                                                    
                                                    

''' 


import optparse
import os
import sys
import ROOT
import numpy
import array
import re
import warnings
import os.path
from math import *
import math

class HH_MvaVarFiller(TreeCloner):

    def __init__(self):
        pass

    def createHH_MVA(self):
        self.getHH_MVAV = ROOT.TMVA.Reader();

        self.getHH_MVAV.AddVariable("mjj_b",                  (self.var1))
        self.getHH_MVAV.AddVariable("deltaR_lep_wjet",        (self.var2))
        self.getHH_MVAV.AddVariable("deltaR_b",               (self.var3))
        self.getHH_MVAV.AddVariable("deltaphi_lep_wjet_high", (self.var4))
        self.getHH_MVAV.AddVariable("deltaR_lep_b",           (self.var5))
        self.getHH_MVAV.AddVariable("deltaeta_lep_wjet_high", (self.var6))
        self.getHH_MVAV.AddVariable("deltaphi_lep_b_high",    (self.var7))
        self.getHH_MVAV.AddVariable("deltaphi_met_wjet_high", (self.var8))
        self.getHH_MVAV.AddVariable("deltaphi_lep_b_low",     (self.var9))

        baseCMSSW = os.getenv('CMSSW_BASE')
        self.getHH_MVAV.BookMVA("BDT","/gwpool/users/achiapparini/CMSSW_8_0_26_patch1/src/LatinoAnalysis/Gardener/python/data/TMVAClassification_BDT.weights.xml")

    def help(self):
        return '''Add mva (BDT) variable'''


    def addOptions(self,parser):
#        description = self.help()
#        group = optparse.OptionGroup(parser,self.label, description)
#        group.add_option('-k', '--kind',   dest='kind', help='Which background training to be used', default='1')
#        parser.add_option_group(group)
#        return group
        pass

    def checkOptions(self,opts):
#        if not (hasattr(opts,'kind')):
#            raise RuntimeError('Missing parameter')
#        self.kind   = opts.kind
#        print " kind = ", self.kind
        pass
    def process(self,**kwargs):
        self.getHH_MVAV = None

        self.var1  = array.array('f',[0])
        self.var2  = array.array('f',[0])
        self.var3  = array.array('f',[0])
        self.var4  = array.array('f',[0])
        self.var5  = array.array('f',[0])
        self.var6  = array.array('f',[0])
        self.var7  = array.array('f',[0])
        self.var8  = array.array('f',[0])
        self.var9  = array.array('f',[0])

        tree  = kwargs['tree']
        input = kwargs['input']
        output = kwargs['output']

        self.connect(tree,input)
        newbranches = ['BDT_var']

        self.clone(output,newbranches)

        BDTmva = numpy.ones(1, dtype=numpy.float32)

        self.otree.Branch('BDT_var',  BDTmva,  'BDT_var' + '/F')

        self.createHH_MVA()

        nentries = self.itree.GetEntries()
        print 'Total number of entries: ',nentries 

        # avoid dots to go faster
        itree     = self.itree
        otree     = self.otree


        print '- Starting eventloop'
        step = 5000
        for i in xrange(nentries):
            itree.GetEntry(i)

            ## print event count
            if i > 0 and i%step == 0.:
                print i,'events processed.'

            BDTmva[0] = -9999.


            pt1 = itree.std_vector_lepton_pt[0]
            pt2 = itree.std_vector_lepton_pt[1]
            
            if pt1>0 : 
                self.var1[0]   =  itree.mjj_b
                self.var2[0]   =  itree.deltaR_lep_wjet
                self.var3[0]   =  itree.deltaR_b
                self.var4[0]   =  itree.deltaphi_lep_wjet_high
                self.var5[0]   =  itree.deltaR_lep_b
                self.var6[0]   =  itree.deltaeta_lep_wjet_high
                self.var7[0]   =  itree.deltaphi_lep_b_high
                self.var8[0]   =  itree.deltaphi_met_wjet_high
                self.var9[0]   =  itree.deltaphi_lep_b_low

                BDTmva[0] = self.getHH_MVAV.EvaluateMVA("BDT")
            print "variabile: ", BDTmva[0]

            otree.Fill()
            
        self.disconnect()
        print '- Eventloop completed'

