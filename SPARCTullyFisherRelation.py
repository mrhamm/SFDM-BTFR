
#enable the below import statement if you wish to render plots in the script
from matplotlib.pylab import plot, show
import pandas as pd

#This is the original data file from SPARC Survey at http://astroweb.cwru.edu/SPARC/BTFR_Lelli2016a.mrt
#Please cite SPARC appropriately if you use this script
Data = pd.read_csv('TullyFisherData.txt')

#this is the format of the text file --- 
#note that the collumn numbers here are INCORRECT
  #1- 11 A11    ---           Name     Galaxy Name
  #12- 17 F6.2   Mpc           Dist     Galaxy Distance
  #18- 23 F5.2   Mpc         e_Dist     Mean error on Distance
  #23- 23 I1     ---         f_Dist     Distance Method (1)
  #24- 28 F5.2   solMass       log(Mb)  Log of the baryonic mass (2)
  #29- 33 F5.2   solMass     e_log(Mb)  Mean error on log(Mb)
  #34- 38 F5.2   km/s          log(Vf)  Log of the rotation velocity (3)
  #39- 43 F5.2   km/s        e_log(Vf)  Mean error on log(Vf)
  #44- 48 F5.2   ---           Fg       Gas fraction Fg=Mg/Mb (2)
  #49- 55 F7.2   solMass/pc2   SBeff    Effective stellar surface density (2)
  #56- 60 F5.2   kpc           Reff     Effective stellar radius

#Builds the dataframe
FramedData = pd.DataFrame(columns = ['Name', 'Distance', 'SigmaD','DMethod','LogM','SigmaLogM','LogV','SigmaLogV'
,'GasFraction','SurfaceDensity','StellarRadius'], index=range(len(Data)))


#populates the dataframe
for i in range(len(Data)):
    FramedData.iloc[i][0] = str(Data.iloc[i][0][0:11])
    FramedData.iloc[i][1] = float(Data.iloc[i][0][12:18])
    FramedData.iloc[i][2] = float(Data.iloc[i][0][19:24])
    FramedData.iloc[i][3] = float(Data.iloc[i][0][25])
    FramedData.iloc[i][4] = float(Data.iloc[i][0][26:32])
    FramedData.iloc[i][5] = float(Data.iloc[i][0][33:38])
    FramedData.iloc[i][6] = float(Data.iloc[i][0][39:44])
    FramedData.iloc[i][7] = float(Data.iloc[i][0][45:50])
    FramedData.iloc[i][8] = float(Data.iloc[i][0][51:56])
    FramedData.iloc[i][9] = float(Data.iloc[i][0][57:64])
    FramedData.iloc[i][10] = float(Data.iloc[i][0][65:70])

#returns the BTFR data for the given filtered data frame
def TullyFisherRelation(Frame):
  logM = Frame['LogM']
  logV = Frame['LogV']
  sigmaM = Frame['SigmaLogM']
  sigmaV = Frame['SigmaLogV']
  return(logV,logM,sigmaV,sigmaM)

#The following test generates the BTFR for galaxies with LogM>10
'''
TestFrame = FramedData.query('LogM > 10')
Test = TullyFisherRelation(TestFrame)
plot(Test[1].iloc[:],Test[0].iloc[:],'o')
show()'''