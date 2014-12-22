from numpy import *

class phaseIdent:
    """Class for identifying phases from peak positions and intensities"""
    def __init__(self):
        self.epsilon = 0.02
        self.miller = dict()
        self.ratios = dict()
        self.fitLimits = dict()

    #test comment

        ###  pm3m miller indices
        pm3m_ind = []
        (self.miller)['pm3m'] = pm3m_ind
        pm3m_ind.append([1,0,0])
        pm3m_ind.append([1,1,0])
        pm3m_ind.append([1,1,1])
        pm3m_ind.append([2,0,0])
        pm3m_ind.append([2,1,0])
        pm3m_ind.append([2,1,1])
        pm3m_ind.append([2,2,0])
        pm3m_ind.append([2,2,1])
        pm3m_ind.append([3,1,0])
             
        pm3mratios = self.calcRatios(pm3m_ind)
      
        (self.ratios)['pm3m'] = pm3mratios
        (self.fitLimits)['pm3m_NumFitPeaks'] = len(pm3m_ind)
        (self.fitLimits)['pm3m_StartPeakLimit'] = 7

        ###  pn3m miller indices
        pn3m_ind = []
        (self.miller)['pn3m'] = pn3m_ind
        pn3m_ind.append([1,1,0])
        pn3m_ind.append([1,1,1])
        pn3m_ind.append([2,0,0])
        pn3m_ind.append([2,1,1])
        pn3m_ind.append([2,2,0])
        pn3m_ind.append([2,2,1])
        pn3m_ind.append([3,1,0])
              
        pn3mratios = self.calcRatios(pn3m_ind)
      
        (self.ratios)['pn3m'] = pn3mratios
        (self.fitLimits)['pn3m_NumFitPeaks'] = len(pn3m_ind)
        (self.fitLimits)['pn3m_StartPeakLimit'] = 7
      
        ### ia3d
        ia3d_ind = []
        (self.miller)['ia3d'] = ia3d_ind
        ia3d_ind.append([2,1,1])
        ia3d_ind.append([2,2,0])
        ia3d_ind.append([3,2,1])
        ia3d_ind.append([4,0,0])
        ia3d_ind.append([4,2,0])
        ia3d_ind.append([3,3,2])
      
        ia3dratios = self.calcRatios(ia3d_ind)
        (self.ratios)['ia3d'] = ia3dratios  
        (self.fitLimits)['ia3d_NumFitPeaks'] = len(ia3d_ind)
        (self.fitLimits)['ia3d_StartPeakLimit'] = 7
      
        ### im3m
        im3m_ind = []
        (self.miller)['im3m'] = im3m_ind
        im3m_ind.append([1,1,0])
        im3m_ind.append([2,0,0])
        im3m_ind.append([2,1,1])
        im3m_ind.append([2,2,0])
        im3m_ind.append([3,1,0])
        im3m_ind.append([2,2,2])
        im3m_ind.append([3,2,1])
      
        im3mratios = self.calcRatios(im3m_ind)
        (self.ratios)['im3m'] = im3mratios  
        (self.fitLimits)['im3m_NumFitPeaks'] = len(im3m_ind)
        (self.fitLimits)['im3m_StartPeakLimit'] = 7
      
        ### Fd3m
        fd3m_ind = []
        (self.miller)['fd3m'] = fd3m_ind
        fd3m_ind.append([1,1,1])
        fd3m_ind.append([2,2,0])
        fd3m_ind.append([3,1,1])
        fd3m_ind.append([2,2,2])
        fd3m_ind.append([4,0,0])
        fd3m_ind.append([3,3,1])
        fd3m_ind.append([4,2,2])
      
        fd3mratios = self.calcRatios(fd3m_ind)
        (self.ratios)['fd3m'] = fd3mratios  
        (self.fitLimits)['fd3m_NumFitPeaks'] = len(fd3m_ind)
        (self.fitLimits)['fd3m_StartPeakLimit'] = 7
        
      
        ### hii
        hii_ind = []
        (self.miller)['hii'] = hii_ind
        hii_ind.append([1,0])
        hii_ind.append([1,1])
        hii_ind.append([2,0])
        hii_ind.append([2,1])
        hii_ind.append([3,0])
          
        hiiratios = self.calcRatios(hii_ind, HEX_2D=True)
        (self.ratios)['hii'] = hiiratios
        (self.fitLimits)['hii_NumFitPeaks'] = len(hii_ind)
        (self.fitLimits)['hii_StartPeakLimit'] = 7
      
        ### l
        l_ind = []
        (self.miller)['l'] = l_ind
        l_ind.append([1,0])
        l_ind.append([2,0])
        l_ind.append([3,0])
          
        lratios = self.calcRatios(l_ind)
        (self.ratios)['l'] = lratios
        (self.fitLimits)['l_NumFitPeaks'] = len(l_ind)
        (self.fitLimits)['l_StartPeakLimit'] = 7
    
    def calcRatios(self, indices, HEX_2D = False):
        
        # convert to numpy array incase input is a list.
        indices = array(indices)
        
        if HEX_2D:
            firstPeak = (sum(indices[0]**2))**0.5
            if indices.size == 2:
                ratios = [(sum(indices[1]**2))**0.5/firstPeak]
            else:
                ratios = [(sum(index**2))**0.5/firstPeak for index in indices[1:]]
        
        else :
            firstPeak = (sum(indices[0]**2) + (indices[0])[0]*(indices[0])[1])**0.5
            if indices.size == 2:
                ratios = (sum(indices[1]**2) + (indices[1])[0]*(indices[1])[1])/firstPeak
            else:
                ratios = [(sum(index**2) + index[0]*index[1])**0.5/firstPeak for index in indices[1:]]
        
        return ratios
    
    def peakAssign(self, peaks, xMax):
        
        # Ensure in numpy array
        peaks = array( peaks)
        
        if peaks.size <= 1:
            return 'Not Enough Peaks'
        
        structureFits = dict()
        for startPeak in range(peaks.size - 2) :
            expPeakRatios = array((peaks/peaks[startPeak])[startPeak + 1:])
            for key in self.ratios:
                structRatios = self.ratios[key]
                if startPeak >= (self.fitLimits)[key + '_StartPeakLimit']:
                    continue
                numFitPeaks = (self.fitLimits)[key + '_NumFitPeaks']
                if numFitPeaks -1  <= 0:
                    continue
                expPeakPresent = empty(numFitPeaks - 1,dtype=int)
                expPeakPresent[:] = -1
                
                for i, ratio in enumerate(structRatios):
                    if i >= numFitPeaks-1:
                        break
                     
                    try:
                        expPeakPresent[i] = where(abs(expPeakRatios - ratio) < self.epsilon)[0][0]
                    except IndexError:
                        pass
                    
                    if i > 1 :
                        try:
                            if (where(expPeakPresent[0:i-1] == expPeakPresent[i]))[0][0] >= 0:
                                expPeakPresent[i] = -1
                        except IndexError:
                            pass
                
                peakPresent = expPeakPresent >= 0
                try:
                    expPeakPresent = array([0,1+expPeakPresent[where(expPeakPresent >= 0)][0]])
                except IndexError:
                    pass
                
                numAccessPeaks = where(peaks[startPeak]*structRatios < xMax)[0].size
                ++numAccessPeaks
                numAccessPeaks = min(numAccessPeaks,numFitPeaks)
                
                #print, 'Matched :' + StrCompress(String(total(peakPresent) + 1,FORMAT = '(i2)')) + ' peaks out of ' + StrCompress(numAccessPeaks) + ' peaks for structure type ' + key + '.'
                #IF Total(peakPresent) GT 0 THEN 
                
                                
                structureFits[key+'_'+str(startPeak).strip()] = { 'TotalPeaks' : peaks.size,'NumPeaksFound' : sum(peakPresent)+1, 'NumAccessiblePeaks' : numAccessPeaks, 'MatchedPeakIndices' : where(peakPresent == 1)[0], 'ExpectedPeakPositions' : peaks[startPeak]*array([1] + structRatios[0:numAccessPeaks-2]),
                                                                  'ExperimentPeakPositions' : peaks[expPeakPresent+startPeak]}

        return structureFits   

    def likelyPhases(self, structureFits):
        fits = zeros((len(structureFits),4))
        labels = []
        for i,key in enumerate(structureFits):
            fit = structureFits[key]
            labels.append(key)
            fits[i,:] = [i,fit['NumPeaksFound'],fit['NumAccessiblePeaks'],fit['TotalPeaks']]
            #fits[i,0] = key
            #fits[i,1] = fit['NumPeaksFound']
            #fits[i,2] = fit['NumAccessiblePeaks']
            #fits[i,3] = fit['TotalPeaks']
        
        print fits