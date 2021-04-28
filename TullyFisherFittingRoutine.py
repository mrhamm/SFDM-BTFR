from matplotlib.pylab import *
import SSSStateGenerator as SG
from matplotlib.pyplot import plot, show, errorbar
from numpy import savetxt, loadtxt, array
from scipy import optimize, stats
import DiskFractionSolver as disk
import BoundaryProblemSolver as BP
import pandas as pd
from SPARCTullyFisherRelation import FramedData, TullyFisherRelation


#rescales the theoretical data
def RescaledBTFR(data, alpha):
    logV = array(data[0]) + 0.5*log10(alpha)
    logM = array(data[1]) + 0.5*log10(alpha)
    return(logV,logM)




#The Following is a basic example of how a fitting routine can work.  This takes about 1-2 hours to process, 
   # increase the NumberOfStates for more data points, decrease it for shorter run time
NumberOfState = 25
m22 = 1
BV = .3
Lambda = 1500
#Amp = BV/Lambda**2
fraction = 0.111
halffraction = 0.65
#observed tully fisher relation
data0 = TullyFisherRelation(FramedData)
errorbar(data0[0],data0[1],xerr=data0[2],yerr=data0[3],fmt='ko')
slope1,intercept1,r1,p1,std1=stats.linregress(data0[0].astype(float),data0[1].astype(float))

#initial theoretical data  
data1 = BP.TullyFisherTrend(6,BV,Lambda,fraction,halffraction,m22)
plot(data1[0],data1[1],'ro')
slope2,intercept2,r2,p2,std2=stats.linregress(data1[0],data1[1])

#finds the factor needed to make the intercepts match
difference = intercept1-intercept2
Alpha = 10**(2*difference/(slope2-1))

#final output "best fit" data
data2 = RescaledBTFR(data1,Alpha)
plot(data2[0],data2[1],'go')
show()