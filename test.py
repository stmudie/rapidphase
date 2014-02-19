#!./env/bin/python
from numpy import *
from peakfinding import peakfinding
from phaseIdent import phaseIdent

x,y = genfromtxt("Well_1B_50_0399.dat",usecols=(0,1), unpack=True)

#y = 1000*y
minx = where(x > 0.08)[0][0]
g = peakfinding(y, peakDescrim = 1, fittingRange=[minx,-1])
k = g.findPeaks()
print x[k]

h = phaseIdent()

peakinfo = h.peakAssign(x[k],x[g.fittingRange[1]-1])



h.likelyPhases(peakinfo)
#intensityMostIntenseFittedPeak = max((peakIntensities)[expPeakPresent+startpeak])
#                intensityMostIntenseUnfittedPeak = max((peakIntensities)[tempArr+startpeak])
