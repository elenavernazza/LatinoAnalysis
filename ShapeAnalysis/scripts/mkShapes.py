#!/usr/bin/env python

import json
import sys
import ROOT
import optparse
#import hwwinfo
#import hwwsamples
#import hwwtools
import LatinoAnalysis.Gardener.hwwtools as hwwtools
import os.path
import string
import logging
import LatinoAnalysis.Gardener.odict as odict
#from HWWAnalysis.Misc.ROOTAndUtils import TH1AddDirSentry
import traceback
from array import array



# ----------------------------------------------------- ShapeFactory --------------------------------------

class ShapeFactory:
    _logger = logging.getLogger('ShapeFactory')
 
    # _____________________________________________________________________________
    def __init__(self):
      
        variables = {}
        self._variables = variables

        cuts = {}
        self._cuts = cuts

        samples = {}
        self._samples = samples




    # _____________________________________________________________________________
    def __del__(self):
        pass

    # _____________________________________________________________________________
    def getvariable(self,tag,mass,cat):

        if tag in self._variables :
            try:
                theVariable = (self._variables[tag])(mass,cat)
            except KeyError as ke:
                self._logger.error('Variable '+tag+' not available. Possible values: '+', '.join(self._variables.iterkeys()) )
                raise ke
        else :
            theVariable = tag

        return theVariable


    # _____________________________________________________________________________
    def makeNominals(self, inputDir, outputDir, variables, cuts, samples, nuisances, supercuts):

        print "======================"
        print "==== makeNominals ===="
        print "======================"
        
        self._variables = variables
        self._samples   = samples
        self._cuts      = cuts
        
        print " supercut = ", supercut
        
        self._outputFileName = outputDir+'/plots_'+self._tag+".root"
        print " outputFileName = ", self._outputFileName
        os.system ("mkdir " + outputDir + "/")
        self._outFile = ROOT.TFile.Open( self._outputFileName, 'recreate')
        
        ROOT.TH1.SetDefaultSumw2(True)
        
        selections = "1"

        #---- first create structure in output root file
        for cutName in self._cuts :
          print "cut = ", cutName, " :: ", cuts[cutName]
          self._outFile.mkdir (cutName)
          for variableName, variable in self._variables.iteritems():
            self._outFile.mkdir (cutName+"/"+variableName)
            #print "variable[name]  = ", variable ['name']
            #print "variable[range] = ", variable ['range']
            #for nuisance in self._nuisances :
              #print "nuisance = ", nuisance
              # open the root file

        # connect the trees
        list_of_trees_to_connect = {}
        for sampleName, sample in self._samples.iteritems():
          list_of_trees_to_connect[sampleName] = sample['name']
              
        inputs = self._connectInputs( list_of_trees_to_connect, inputDir)
                   

        ## test with LOOP: not very fast ... to be improved ...

        ##---- now plot and save into output root file
        ## - loops
        ##    -> samples
        ##      -> cut 
        ##        -> variables
        ##

        ## all the histograms
        #histos = {}       
        #weightToPlotFormulas = {}
        #variableToPlotFormulas = {}
        
        #for sampleName, sample in self._samples.iteritems():
          
          #print "sampleName = ", sampleName
          ## first prepare all dummy histograms
          #for cutName, cut in self._cuts.iteritems():
            #for variableName, variable in self._variables.iteritems():

              #numTree = 0
              #for tree in inputs[sampleName]:      
                ## new histogram
                #shapeName = 'histo_' + sampleName + "_" + cutName + "_" + variableName + "_" + str(numTree)
                ## prepare a dummy to fill
                #histos[shapeName] = self._makeshape(shapeName,variable['range'])


                #global_weight = sample ['weight']
                #globalCut = "(" + cut + ") * (" + global_weight + ")"  

                ## if weights vector is not given, do not apply file dependent weights
                #if 'weights' in sample.keys() and len(sample ['weights']) != 0 :
                  ## if weight is not given for a given root file, '-', do not apply file dependent weight for that root file
                  #if sample ['weights'][numTree] != '-' :
                    #globalCut = "(" + globalCut + ") * (" +  sample ['weights'][numTree] + ")" 

                #weightFormulaName = 'weight_' + sampleName + '_' + cutName + '_' + variableName + '_' + str(numTree)
                #weightToPlotFormulas[weightFormulaName] = ROOT.TTreeFormula(weightFormulaName, globalCut , tree)

                #numTree += 1


        #for sampleName, sample in self._samples.iteritems():          
          #for variableName, variable in self._variables.iteritems():
            #numTree = 0
            #for tree in inputs[sampleName]:      
              #variableFormulaName = 'weight_' + sampleName + '_' + '_' + variableName + '_' + str(numTree)
              #variableToPlotFormulas[variableFormulaName] = ROOT.TTreeFormula(variableFormulaName, variable['name'] , tree)
              #numTree += 1


          
          ## now really fill the histograms
          #numTree = 0
          #for tree in inputs[sampleName]:
           
            #print '        {0:<20} : {1:^9}'.format(sampleName,tree.GetEntries()),
            ##print '        {0:<20} : {1:^9}'.format(sampleName,tree.GetEntries())
            ## loop over events
            #step = 5000
            #nentries = tree.GetEntries()
            
            ## pre-filter
            ##  to speed up
            #tree.SetEntryList(0)
            ## get the list
            ##myList = ROOT.TEntryList(tree)
            #myList = ROOT.TEntryList('myList'+'_'+str(numTree)+'_'+sampleName,"")
            #nentriesFiltered = myList.GetN()
            #print "      -> nentriesFiltered = ", nentriesFiltered, 
            #tree.Draw('>> myList'+'_'+str(numTree)+'_'+sampleName , supercut, "entrylist");
            ## apply the list
            #tree.SetEntryList(myList)
            #nentriesFiltered = myList.GetN()
            #print "      -> nentriesFiltered = ", nentriesFiltered
 

            ## mild the steps
            #while nentriesFiltered / step > 10 :
              #step *= 10


            ## now really looping over events
            ##for iEvent in xrange(nentries):
            
            #for iEventSel in xrange(nentriesFiltered):
              #iEvent = -1
              #if iEventSel == 0 :
               #iEvent = myList.GetEntry(0)
               #iEventSel+=1
              #else :
               #iEvent = myList.Next()
               #iEventSel+=1
               #tree.GetEntry(iEvent)
              
              ### print event count
              #if iEventSel > 0 and iEventSel%step == 0.:
                #print '   >> ', iEvent, 'events processed ::', nentries, ' [', iEventSel, '(', nentriesFiltered, ')] --> %.2f' % (1. * iEventSel / nentriesFiltered * 100), ' % '

              #for cutName, cut in self._cuts.iteritems():
                ##print "cut = ", cutName, " :: ", cut

                ##global_weight = sample ['weight']
                ##globalCut = "(" + cut + ") * (" + global_weight + ")"  

                ### if weights vector is not given, do not apply file dependent weights
                ##if 'weights' in sample.keys() and len(sample ['weights']) != 0 :
                  ### if weight is not given for a given root file, '-', do not apply file dependent weight for that root file
                  ##if sample ['weights'][numTree] != '-' :
                    ##globalCut = "(" + globalCut + ") * (" +  sample ['weights'][numTree] + ")" 

                ##weightToPlotFormula = ROOT.TTreeFormula('blabla', globalCut , tree)
                
                #for variableName, variable in self._variables.iteritems():
                  ##print "  variable[name]  = ", variable['name']
                  ##print "  variable[range] = ", variable['range']
                  #self._outFile.cd (cutName+"/"+variableName)

                  ##variableToPlotFormula = ROOT.TTreeFormula('blabla', variable['name'] , tree)

                  #shapeName = 'histo_' + sampleName + '_' + cutName + '_' + variableName + '_' + str(numTree)
              
                  #self._logger.debug('---'+sampleName+'---')
                  #self._logger.debug('Formula: '+variable['name']+'>>'+shapeName)
                  #self._logger.debug('Cut:     '+cut)
                  #self._logger.debug('ROOTFiles:'+'\n'.join([f.GetTitle() for f in tree.GetListOfFiles()]))

                  #weightFormulaName = 'weight_' + sampleName + '_' + cutName + '_' + variableName + '_' + str(numTree)
                  #variableFormulaName = 'weight_' + sampleName + '_' + '_' + variableName + '_' + str(numTree)

                  #histos[shapeName].Fill( variableToPlotFormulas[variableFormulaName].EvalInstance(), weightToPlotFormulas[weightFormulaName].EvalInstance() )
 
            #numTree += 1
          
        ##
        ## now post filling processing
        ## loop
        ##  -> cut
        ##    -> variable
        ##      -> sample
        #for cutName, cut in self._cuts.iteritems():
          #for variableName, variable in self._variables.iteritems():      
            #for sampleName, sample in self._samples.iteritems():
              #bigName = 'histo_' + sampleName
              #hTotal = self._makeshape(bigName,variable['range'])
              #numTree = 0
              #for tree in inputs[sampleName] :
                #shapeName = 'histo_' + sampleName + "_" + cutName + "_" + variableName + "_" + str(numTree)
                #if (numTree == 0) :
                  #histos[shapeName].SetTitle(bigName)
                  #histos[shapeName].SetName(bigName)
                  #hTotal = histos[shapeName]
                #else :
                  #hTotal.Add(histos[shapeName])

                #numTree += 1
 

              ## fold if needed
              #doFold = 0
              #if 'fold' in variable.keys() :
                ##print "    variable[fold] = ", variable ['fold']
                #doFold = variable ['fold']

              #if doFold == 1 or doFold == 3 :
                #self._FoldOverflow  (hTotal)
              #if doFold == 2 or doFold == 3 :
                #self._FoldUnderflow (hTotal)
        
        
              ## go 1d
              #self._outFile.cd (cutName+"/"+variableName)
              #outputsHisto = self._h2toh1(hTotal)

              ## eventually write to root file!
              #outputsHisto.Write()              
     
            
              ## prepare nuisance MC/data statistics
              ## - uniform
              ## - uniform method 2
              ## - bin by bin (in selected bins)
              #for nuisanceName, nuisance in nuisances.iteritems():
                #if nuisanceName == 'stat' : # 'stat' has a separate treatment, it's the MC/data statistics
                  ##print "nuisance[type] = ", nuisance ['type']
                  #for sampleNuisName, configurationNuis in nuisance['samples'].iteritems() :
                    #if sampleNuisName == sampleName: # check if it is the sample I'm analyzing!
                      #if configurationNuis['typeStat'] == 'uni' :
                        ##print "     >> uniform"
                        ## take histogram --> outputsHisto
                        #outputsHistoUp = outputsHisto.Clone("histo_"+sampleName+"_statUp")
                        #outputsHistoDo = outputsHisto.Clone("histo_"+sampleName+"_statDown")
                        ## scale up/down
                        #self._scaleHistoStat (outputsHistoUp,  1 )
                        #self._scaleHistoStat (outputsHistoDo, -1 )
                        ## save the new two histograms in final root file
                        #outputsHistoUp.Write()
                        #outputsHistoDo.Write()
                                 
        ## - then disconnect the files
        #self._disconnectInputs(inputs)
        
        
        #################################
        # old method ... but working ...
        
        #---- now plot and save into output root file
        for cutName, cut in self._cuts.iteritems():
          print "cut = ", cutName, " :: ", cut

          # create the list of events -> speed up!          
          # for each tree!!!
          for sampleName, sample in self._samples.iteritems():
            if 'weights' in sample.keys() :
              self._filterTrees( sample ['weight'], sample ['weights'], '(' + cut + ') && (' + supercut + ')' , inputs[sampleName], cutName, sampleName)
            else :
              self._filterTrees( sample ['weight'], []                , '(' + cut + ') && (' + supercut + ')' , inputs[sampleName], cutName, sampleName)

           
          for variableName, variable in self._variables.iteritems():
            print "  variable[name]  = ", variable['name']
            print "  variable[range] = ", variable['range']
            self._outFile.cd (cutName+"/"+variableName)
            
            for sampleName, sample in self._samples.iteritems():
              print "    sample[name]    = ", sample ['name']
              print "    sample[weight]  = ", sample ['weight']
              if 'weights' in sample.keys() :
                print "    sample[weights] = ", sample ['weights']

              doFold = 0
              if 'fold' in variable.keys() :
                print "    variable[fold] = ", variable ['fold']
                doFold = variable ['fold']
              
              # create histogram: already the "hadd" of possible sub-contributions
              if 'weights' in sample.keys() :
                outputsHisto = self._draw( variable['name'], variable['range'], sample ['weight'], sample ['weights'], cut, sampleName, inputs[sampleName], doFold)
              else :
                outputsHisto = self._draw( variable['name'], variable['range'], sample ['weight'], [],                 cut, sampleName, inputs[sampleName], doFold)
               
              outputsHisto.Write()              
              
            
              # prepare nuisance MC/data statistics
              # - uniform
              # - uniform method 2
              # - bin by bin (in selected bins)
              for nuisanceName, nuisance in nuisances.iteritems():
                if nuisanceName == 'stat' : # 'stat' has a separate treatment, it's the MC/data statistics
                  #print "nuisance[type] = ", nuisance ['type']
                  for sampleNuisName, configurationNuis in nuisance['samples'].iteritems() :
                    if sampleNuisName == sampleName: # check if it is the sample I'm analyzing!
                      if configurationNuis['typeStat'] == 'uni' :
                        #print "     >> uniform"
                        # take histogram --> outputsHisto
                        outputsHistoUp = outputsHisto.Clone("histo_"+sampleName+"_statUp")
                        outputsHistoDo = outputsHisto.Clone("histo_"+sampleName+"_statDown")
                        # scale up/down
                        self._scaleHistoStat (outputsHistoUp,  1 )
                        self._scaleHistoStat (outputsHistoDo, -1 )
                        # save the new two histograms in final root file
                        outputsHistoUp.Write()
                        outputsHistoDo.Write()
                        
          
            #for nuisance in self._nuisances :
              #print "nuisance = ", nuisance
              #for sample in self._samples :
                 #print "sample = ", sample
                 
                 # get the weight
                 # open the root file
                 # plot
                 #self._draw(doalias, rng, selections, output, inputs)
                 # save histogram

                 #print 'Output file:',self._outputFile

        # - then disconnect the files
        self._disconnectInputs(inputs)

        


    # _____________________________________________________________________________
    def _filterTrees(self, global_weight, weights, cut, inputs, sampleName, cutName):       
        '''
        global_weight :   the global weight for the samples
        weights       :   the wieghts 'root file' dependent
        cut           :   the selection
        inputs        :   the list of input files for this particular sample
        '''
        self._logger.info('filter Trees to speed up')

        numTree = 0

        for tree in inputs:
          globalCut = "(" + cut + ") * (" + global_weight + ")"  
          # if weights vector is not given, do not apply file dependent weights
          if len(weights) != 0 :
            # if weight is not given for a given root file, '-', do not apply file dependent weight for that root file
            if weights[numTree] != '-' :
              globalCut = "(" + globalCut + ") * (" +  weights[numTree] + ")" 
          
          
          #print " ::: ", tree.GetEntries(),
          # clear list
          tree.SetEntryList(0)
          # get the list
          myList = ROOT.TEntryList('myList'+'_'+str(numTree)+'_'+sampleName+'_'+cutName,"")
          #myList = ROOT.TEntryList(tree)
          tree.Draw('>> myList'+'_'+str(numTree)+'_'+sampleName+'_'+cutName, globalCut, "entrylist");
          #gDirectory = ROOT.gROOT.GetGlobal("gDirectory")
          #myList = gDirectory.Get("myList")
          # apply the list
          tree.SetEntryList(myList)
          
          #print " --> ", tree.GetEntries()
          
          numTree += 1


          
    # _____________________________________________________________________________
    def _draw(self, var, rng, global_weight, weights, cut, sampleName, inputs, doFold):       
        '''
        var           :   the variable to plot
        rng           :   the variable to plot
        global_weight :   the global weight for the samples
        weights       :   the wieghts 'root file' dependent
        cut           :   the selection
        inputs        :   the list of input files for this particular sample
        '''
        
        self._logger.info('Yields by process')
  
        numTree = 0
        bigName = 'histo_' + sampleName
        hTotal = self._makeshape(bigName,rng)
        for tree in inputs:
          print '        {0:<20} : {1:^9}'.format(sampleName,tree.GetEntries()),
          # new histogram
          shapeName = 'histo_' + sampleName + str(numTree)

          # prepare a dummy to fill
          shape = self._makeshape(shapeName,rng)

          self._logger.debug('---'+sampleName+'---')
          self._logger.debug('Formula: '+var+'>>'+shapeName)
          self._logger.debug('Cut:     '+cut)
          self._logger.debug('ROOTFiles:'+'\n'.join([f.GetTitle() for f in tree.GetListOfFiles()]))

          globalCut = "(" + cut + ") * (" + global_weight + ")"  
          # if weights vector is not given, do not apply file dependent weights
          if len(weights) != 0 :
            # if weight is not given for a given root file, '-', do not apply file dependent weight for that root file
            if weights[numTree] != '-' :
              globalCut = "(" + globalCut + ") * (" +  weights[numTree] + ")" 
            
          # in principle now that the trees are filtered
          # I may remove the globalCut here
          # ... but it doesn't hurt leaving it
          entries = tree.Draw( var+'>>'+shapeName, globalCut, 'goff')
          #shape = (ROOT.TH1D*) gDirectory->Get(shapeName)
          print '     >> ',entries,':',shape.Integral()
          
          if (numTree == 0) :
            shape.SetTitle(bigName)
            shape.SetName(bigName)
            hTotal = shape
          else :
            hTotal.Add(shape)

          numTree += 1

        print ' ~~~~ '
        
        # fold if needed
        if doFold == 1 or doFold == 3 :
          self._FoldOverflow  (hTotal)
        if doFold == 2 or doFold == 3 :
          self._FoldUnderflow (hTotal)
        
        
        # go 1d
        hTotalFinal = self._h2toh1(hTotal)
        
        return hTotalFinal



    def _FoldUnderflow(self, h):       

        #print " >> fold underflow"
        if h.GetDimension() == 1:
          nx = h.GetNbinsX()
          # 0 --> 1
          ShapeFactory._moveAddBin(h, (0,),(1,) )
          return
        elif h.GetDimension() == 2:
          nx = h.GetNbinsX()
          ny = h.GetNbinsY()
          for i in xrange(1,nx+1):
            ShapeFactory._moveAddBin(h,(i,0   ),(i, 1 ) )

          for j in xrange(1,ny+1):
            ShapeFactory._moveAddBin(h,(0,    j),(1, j) )

          # 0,0 -> 1,1
          # 0,ny+1 -> 1,ny+1
          # nx+1,0 -> nx+1,1
          
          ShapeFactory._moveAddBin(h, (0,0),(1,1) )
          ShapeFactory._moveAddBin(h, (0,ny+1),(1,ny+1) )
          ShapeFactory._moveAddBin(h, (nx+1,0),(nx+1,1) )
          
        
        
    def _FoldOverflow(self, h):       

        #print " >> fold overflow"
        if h.GetDimension() == 1:
          nx = h.GetNbinsX()
          # n+1 --> n
          ShapeFactory._moveAddBin(h, (nx+1,),(nx,) )
          return
        elif h.GetDimension() == 2:
          nx = h.GetNbinsX()
          ny = h.GetNbinsY()
          for i in xrange(1,nx+1):
            ShapeFactory._moveAddBin(h,(i,ny+1),(i, ny) )

          for j in xrange(1,ny+1):
            ShapeFactory._moveAddBin(h,(nx+1, j),(nx,j) )

            # 0,ny+1 -> 0,ny
            # nx+1,0 -> nx,0
            # nx+1,ny+1 ->nx,ny

            ShapeFactory._moveAddBin(h, (0,ny+1),(0,ny) )
            ShapeFactory._moveAddBin(h, (nx+1,0),(nx,0) )
            ShapeFactory._moveAddBin(h, (nx+1,ny+1),(nx,ny) )


# --- transform 2D into 1D: unrolling
# 
#      3    6    9
#      2    5    8
#      1    4    7
#
    def _h2toh1(self, h):
        import array
        
        if not isinstance(h,ROOT.TH2):
            return h
           
        sentry = TH1AddDirSentry()

#         H1class = getattr(ROOT,h.__class__.__name__.replace('2','1'))
        nx = h.GetNbinsX()
        ny = h.GetNbinsY()

        h_flat = ROOT.TH1D(h.GetName(),h.GetTitle(),nx*ny,0,nx*ny)
 
        sumw2 = h.GetSumw2()
        sumw2_flat = h_flat.GetSumw2()

        for i in xrange(1,nx+1):
            for j in xrange(1,ny+1):
                # i,j must be mapped in 
                b2d = h.GetBin( i,j )
#                 b2d = h.GetBin( j,i )
#                 b1d = ((i-1)+(j-1)*nx)+1
                b1d = ((j-1)+(i-1)*ny)+1

                h_flat.SetAt( h.At(b2d), b1d )
                sumw2_flat.SetAt( sumw2.At(b2d), b1d ) 

        h_flat.SetEntries(h.GetEntries())
        
        stats2d = array.array('d',[0]*7)
        h.GetStats(stats2d)

        stats1d = array.array('d',[0]*4)
        stats1d[0] = stats2d[0]
        stats1d[1] = stats2d[1]
        stats1d[2] = stats2d[2]+stats2d[4]
        stats1d[3] = stats2d[3]+stats2d[5]

        h_flat.PutStats(stats1d)

        xtitle = h.GetXaxis().GetTitle()
        v1,v2 = xtitle.split(':') # we know it's a 2d filled by an expr like y:x
        xtitle = '%s #times %s bin' % (v1,v2)

        h_flat.GetXaxis().SetTitle(xtitle)

        return h_flat
       
 
 
    # _____________________________________________________________________________
    def _scaleHistoStat(self, histo, direction):
        
        for iBin in range(1, histo.GetNbinsX()+1):
          error = histo.GetBinError(iBin)
          value = histo.GetBinContent(iBin)
          newvalue = value + direction * error
          histo.SetBinContent(iBin, newvalue)
  
  
  

    # _____________________________________________________________________________
    @staticmethod
    def _moveAddBin(h, fromBin, toBin ):
        if not isinstance(fromBin,tuple) or not isinstance(toBin,tuple):
            raise ValueError('Arguments must be tuples')

        dims = [h.GetDimension(), len(fromBin), len(toBin) ]

        if dims.count(dims[0]) != len(dims):
            raise ValueError('histogram and the 2 bins don\'t have the same dimension')
        
        #print " h = " , h.GetTitle(), " --> ", h.GetNbinsX()
        
        # get bins
        b1 = h.GetBin( *fromBin )
        b2 = h.GetBin( *toBin )

        # move contents
        c1 = h.At( b1 )
        c2 = h.At( b2 )

        h.SetAt(0, b1)
        h.SetAt(c1+c2, b2)

        # move weights as well
        sumw2 = h.GetSumw2()

        w1 = sumw2.At( b1 )
        w2 = sumw2.At( b2 )

        sumw2.SetAt(0, b1)
        sumw2.SetAt(w1+w2, b2)


    # _____________________________________________________________________________
    @staticmethod
    def _bins2hclass( bins ):
        '''
        Fixed bin width
        bins = (nx,xmin,xmax)
        bins = (nx,xmin,xmax, ny,ymin,ymax)
        Variable bin width
        bins = ([x0,...,xn])
        bins = ([x0,...,xn],[y0,...,ym])  
        '''

        from array import array
        if not bins:
            return name,0
        elif not ( isinstance(bins, tuple) or isinstance(bins,list)):
            raise RuntimeError('bin must be an ntuple or an arrays')

        l = len(bins)
        # 1D variable binning
        if l == 1 and isinstance(bins[0],list):
            ndim=1
            hclass = ROOT.TH1D
            xbins = bins[0]
            hargs = (len(xbins)-1, array('d',xbins))
        elif l == 2 and  isinstance(bins[0],list) and  isinstance(bins[1],list):
            ndim=2
            hclass = ROOT.TH2D
            xbins = bins[0]
            ybins = bins[1]
            hargs = (len(xbins)-1, array('d',xbins),
                    len(ybins)-1, array('d',ybins))
        elif l == 3:
            # nx,xmin,xmax
            ndim=1
            hclass = ROOT.TH1D
            hargs = bins
        elif l == 6:
            # nx,xmin,xmax,ny,ymin,ymax
            ndim=2
            hclass = ROOT.TH2D
            hargs = bins
        else:
            # only 1d or 2 d hist
            raise RuntimeError('What a mess!!! bin malformed!')
        
        return hclass,hargs,ndim

    @staticmethod
    def _bins2dim(bins):
        hclass,hargs,ndim = ShapeFactory._bins2hclass( bins )
        return ndim

    @staticmethod
    def _makeshape( name, bins ):
        hclass,hargs,ndim = ShapeFactory._bins2hclass( bins )
        return hclass(name, name, *hargs)
      
       
 
    # _____________________________________________________________________________
    def _connectInputs(self, samples, inputDir):
        inputs = {}
        treeName = 'latino'
        for process,filenames in samples.iteritems():
          tree = self._buildchain(treeName,[ (inputDir + '/' + f) for f in filenames])
          inputs[process] = tree
          # FIXME: add possibility to add Friend Trees for new variables   
         
        return inputs

    # _____________________________________________________________________________
    def _disconnectInputs(self,inputs):
        for n in inputs.keys():
          #friends = inputs[n].GetListOfFriends()
          #if friends.__nonzero__():
            #for fe in friends:
              #friend = fe.GetTree()
              #inputs[n].RemoveFriend(friend)
              #ROOT.SetOwnership(friend,True)
              #del friend
          # remove the entire list of trees
          del inputs[n]
    
    # _____________________________________________________________________________
    def _buildchain(self,treeName,files):
        listTrees = []
        for path in files:
            self._logger.debug('     '+str(os.path.exists(path))+' '+path)
            if not os.path.exists(path):
                raise RuntimeError('File '+path+' doesn\'t exists')
            tree = ROOT.TChain(treeName)
            tree.Add(path)
            listTrees.append(tree)
        return listTrees



if __name__ == '__main__':
    print '''
--------------------------------------------------------------------------------------------------

   ___|   |                               \  |         |                
 \___ \   __ \    _` |  __ \    _ \      |\/ |   _` |  |  /   _ \   __| 
       |  | | |  (   |  |   |   __/      |   |  (   |    <    __/  |    
 _____/  _| |_| \__,_|  .__/  \___|     _|  _| \__,_| _|\_\ \___| _|    
                       _|                                               

--------------------------------------------------------------------------------------------------
'''    
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)

    parser.add_option('--tag'            , dest='tag'            , help='Tag used for the shape file name'           , default=None)
    parser.add_option('--sigset'         , dest='sigset'         , help='Signal samples [SM]'                        , default='SM')
    parser.add_option('--outputDir'      , dest='outputDir'      , help='output directory'                           , default='./')
    parser.add_option('--inputDir'       , dest='inputDir'       , help='input directory'                            , default='./data/')
    parser.add_option('--nuisancesFile'  , dest='nuisancesFile'  , help='file with nuisances configurations'         , default=None )

          
    # read default parsing options as well
    hwwtools.addOptions(parser)
    hwwtools.loadOptDefaults(parser)
    (opt, args) = parser.parse_args()

    sys.argv.append( '-b' )
    ROOT.gROOT.SetBatch()


    print " configuration file = ", opt.pycfg
    print " lumi =               ", opt.lumi
    
    print " inputDir =           ", opt.inputDir
    print " outputDir =          ", opt.outputDir
 
    

    if not opt.debug:
        pass
    elif opt.debug == 2:
        print 'Logging level set to DEBUG (%d)' % opt.debug
        logging.basicConfig( level=logging.DEBUG )
    elif opt.debug == 1:
        print 'Logging level set to INFO (%d)' % opt.debug
        logging.basicConfig( level=logging.INFO )

      
    factory = ShapeFactory()
    factory._energy    = opt.energy
    factory._lumi      = opt.lumi
    factory._tag       = opt.tag
    
    
    variables = {}
    if os.path.exists(opt.variablesFile) :
      handle = open(opt.variablesFile,'r')
      exec(handle)
      handle.close()
    
    supercut = '1'
    cuts = {}
    if os.path.exists(opt.cutsFile) :
      handle = open(opt.cutsFile,'r')
      exec(handle)
      handle.close()
    
    samples = {}
    if os.path.exists(opt.samplesFile) :
      handle = open(opt.samplesFile,'r')
      exec(handle)
      handle.close()
    

    # ~~~~
    nuisances = {}
    if opt.nuisancesFile == None :
      print " Please provide the nuisances structure if you want to add nuisances "      
    elif os.path.exists(opt.nuisancesFile) :
        handle = open(opt.nuisancesFile,'r')
        exec(handle)
        handle.close()
    
    
    factory.makeNominals( opt.inputDir ,opt.outputDir, variables, cuts, samples, nuisances, supercut)
    
        
        
        
        
        
        
        
        

    #except Exception as e:
        #print '*'*80
        #print 'Fatal exception '+type(e).__name__+': '+str(e)
        #print '*'*80
        #exc_type, exc_value, exc_traceback = sys.exc_info()
        #traceback.print_tb(exc_traceback, file=sys.stdout)
##         traceback.print_tb(exc_traceback, limit=3, file=sys.stdout)
        #print '*'*80
    #finally:
        #print 'Used options'
        #print ', '.join([ '{0} = {1}'.format(a,b) for a,b in opt.__dict__.iteritems()])