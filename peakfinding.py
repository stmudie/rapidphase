from numpy import *

class peakfinding:
    """class for finding peaks"""
    def __init__(self,y, peakDescrim = 1.0, fittingRange = [-1,-1]):

        self.y = y
        self.peakDescrim = peakDescrim
        fittingRange[0] = fittingRange[0] if fittingRange[0] != -1 else 0
        fittingRange[1] = fittingRange[1] if fittingRange[1] != -1 else self.y.size
        self.fittingRange = fittingRange
        self.peakIndicies = []
        self.peakIntensities = []
    
    
    def findPeaks(self):
    
        y = self.y[self.fittingRange[0]:self.fittingRange[1]]
        d0 = (y-roll(y,1))[1:-2]
        d1 = (y-roll(y,-1))[1:-2]
        peaks = array([i+1 for i in range(d0.size) if d0[i] > 0 and d1[i] > 0])
        troughs = array([i+1 for i in range(d0.size) if d0[i] < 0 and d1[i] < 0])
        
        for pos in peaks :

            try:
                trLowPos = (where(troughs < pos)[0])[-1]
            except IndexError:
                trLowPos = 0
            trLow = troughs[trLowPos]
            if trLowPos == troughs.size - 1:
                continue
            trHigh = troughs[trLowPos + 1]
            yback = y[trLow] + (y[trHigh]-y[trLow])*(pos-trLow)/(trHigh-trLow) 
            if y[pos] - yback > self.peakDescrim :
                self.peakIndicies.append(pos + self.fittingRange[0])
                self.peakIntensities.append(y[pos] - yback)
            
        
        return self.peakIndicies