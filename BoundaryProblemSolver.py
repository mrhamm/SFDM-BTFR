from matplotlib.pylab import *
from matplotlib.pyplot import plot, show
from numpy import savetxt
import DiskFractionSolver as disk

#This file solves a boundary value problem for fixing the amplitude-wavelength value of 
# galaxies with a SSS state as the dark matter component.  

     

def Invariant(data):
    inv = zeros(len(data[0]))
    for j in range(len(inv)):
        inv[j] = abs(data[-2][j])*(data[-3][j])**2
    return inv 


#this will find the intersection point of the halo's invariant with the value of I
# data input is meant to be output in form of disk.ExponentialGalaxy
def BoundaryIntersection(data,I):
    HaloInvariant = Invariant(data)
    Difference = abs(HaloInvariant - I)
    Difference2 = Difference[int(len(Difference)/2):]
    index = where( Difference == min(Difference2) )
    return(index)

#this adds physical units to the results  
def AddUnits(R,Psi,M,Mext,Mtot,Vcirc,m22):
    Upsilon=100
    Mass =(1E-22)*m22#quote mass in eV/c^2
    Wavenumber = Mass / (1.973E-7)  #give wavenumber in inverse meter (M/hbarc where hbarc =197.3evnm)
    Wavenumber = Wavenumber * (9.461E15) #converts to inverse light years
    Conversion = Wavenumber / Upsilon
    Lengthscale = 1.0/Conversion  #length scale in Ly
    Lengthscale2 = Lengthscale / 3262.  #length scale in kpc
    Lengthscale3 = Lengthscale*9.461E15 #length in meter
    Timescale = Lengthscale #timescale in years
    Massscale = Lengthscale*(6.41E12) #In solar masses
    Massscale2 = Massscale*1.988E30 #in kg
    DensityScale = Massscale2/(Lengthscale3)**3
    bigG = 6.674E-11 #Newton constant
    SpeedOfLight = 2.99E8 #speed of light m/2
    return(Lengthscale2*R,sqrt(DensityScale)*Psi,Massscale*M,Massscale*Mext,
    Massscale*Mtot,SpeedOfLight/1000.0 *Vcirc)



#BV = .3
#Lambda = 1500
#Amp = BV/Lambda**2
#fraction = 0.111

# produces a trend by fixing the boundary valuewavelength 
def TullyFisherTrend(N,BV,Wavel,BaryonFraction,DesiredHalfFraction,m22):
    Amp = BV/Wavel**2
    fraction = BaryonFraction
    Fguess =0
    Vguess =0
    Cguess =0
    Kguess =0
    Solution = []
    index = []
    logM = []
    logV = []
    for i in range(N):
        print(i)
        if i ==0:
            Data = disk.ExponentialGalaxy(False,False,0,DesiredHalfFraction,fraction,False,False)
        else: 
            Data = disk.ExponentialGalaxy(Fguess,Vguess,i,DesiredHalfFraction,fraction,Cguess,Kguess)
        Fguess = Data[0][1][0]
        Vguess = Data[0][2][0]
        Cguess = Data[1]
        Kguess = Data[2]
        Solution_i = [Fguess,Vguess,Cguess,Kguess]
        Solution = Solution+[Solution_i]
        intersection = BoundaryIntersection(Data[0],BV)
        index = index + [intersection]
    
        if i>0:
            Halo = disk.HaloData(Data[0],Cguess,Kguess)
            factor = Amp/Data[0][-2][index[i][0][0]]
            radius = array(Data[0][0])*1.0/sqrt(factor)
            MDM = array(Data[0][3])*sqrt(factor)
            Psi = array(Data[0][1])*factor
            Mext = array(Halo[0])*sqrt(factor)
        #print(Mext[-1])
        #print(Mtot[-1])
            Mtot = array(Halo[1])*sqrt(factor)
            Vcirc = array(Halo[2])*sqrt(factor)
            data = AddUnits(radius,Psi,MDM,Mext,Mtot,Vcirc,m22)
            Vcirc = data[-1]
            Mext = data[3]
            logV = logV + [log10(Vcirc[int(len(Vcirc)/2)])]
            logM = logM + [log10(Mext[int(len(Mext)/2)])]
    return(logV,logM)
    


#BV = .3
#Lambda = 1500
#Amp = BV/Lambda**2
#fraction = 0.111

#data =  TullyFisherTrend(12,0.3,1500,0.111,0.65,1.0)

#plot(data[0],data[1],'o')
#show()
  
