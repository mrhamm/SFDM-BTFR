from matplotlib.pylab import *
import SSSStateGenerator as SG
from matplotlib.pyplot import plot, show
from numpy import savetxt



#This file contains a specialized computation. 
# For mathematical details see chapter 3 of https://arxiv.org/pdf/2004.07792.pdf



dr = .02 #should be 0.02 to be consistent with the SSSStateGenerator package



NumAlpha = 5  #this variable is used to adjust the granularity of the continuation parameter alpha, which appears in the "Exponential Galaxy" routine
                # in the future, one should consider a more intelligent way to choose this parameter



#this is the spherical exponential potential
def Vext(C,K,r):
    if r>0:
        return -K/C**3*(2/r*(1-exp(-C*r))-C*exp(-C*r))
    else: 
        return -K/C**2
        
#the corresponding density to the potential
def Rhoext(C,K,r):
    if r>0:
        return K*exp(-C*r)
    else:
        return K 

#integrates Rhoext up to r
def Mext(C,K,r):
    if r == 0:
        return 0
    else:
        Sum = 0 
        radii = linspace(0,r,int(r/dr)+2)
        dR = radii[1]
        for j in range(len(radii)):
            Sum = Sum + dR*Rhoext(C,K,radii[j])*radii[j]*radii[j]*4*pi
        return Sum        

#some relevant physical information
#newdata is a vector in the format of the output of SG.PhysicalState routing
#C0 and K0 are the parameters of the background exponential distribution. 
def HaloData(newdata,C0,K0):
    Mextdata = zeros(len(newdata[0]))
    Mtot = zeros(len(newdata[0]))
    Vrot = zeros(len(newdata[0]))
    fdm = zeros(len(newdata[0]))
    fext = zeros(len(newdata[0]))
    extdensity = zeros(len(newdata[0]))
    totdensity = zeros(len(newdata[0]))
    for j in range(1,len(Mextdata)):
        Mextdata[j] = Mext(C0,K0,newdata[0][j])
        Mtot[j] = newdata[3][j]+Mextdata[j]
        Vrot[j] = sqrt(Mtot[j]/newdata[0][j])
        fdm[j] = newdata[3][j]/Mtot[j]
        fext[j] = Mextdata[j]/Mtot[j]
        extdensity[j] = Rhoext(C0,K0,newdata[0][j])
        totdensity[j] = extdensity[j] + newdata[1][j]*newdata[1][j]*2 
    Rhalf = 2.67/C0
    index = int(Rhalf/dr)
    return(Mextdata,Mtot,Vrot,fdm,fext,extdensity,totdensity,Rhalf,index)

#The difficulty of this computation scales linearly with N,
# be cautious about computing large values of N.  We suggest one chooses
# educated values for the guesses of FG,VG,CG,KG or else risk a large
# increase in computation time.  Successive values of N have solutions 
# that are 'nearby' in parameter space; the output values of F,G,C,K
# for a state of order N will be suitable guesses for FG,VG,CG,KG

#A sufficiently poor guess for FG,VG,CK, or KG will result in non-convergence
# or an error throw.  In the case that you have no guess, set the values
# of FG,VG,CG,KG to false and generic values will be used instead.         
# 
# DesiredHalfFraction is the total DM fraction at the location of the external contribution half-radius.  

#outputs parameter values for both the DM halo and external component.  
def ExponentialGalaxy(FG,VG,N,DesiredHalfFraction,BaryonFraction,CG,KG):

    DMO = SG.DMOExcitedState(N,1,25,False,FG,VG)[0]
    DMOMass = DMO[3][-1]
    
    if CG == False:
        RDM0 = DMO[-6]
        widthfraction = 1
        C0 = 2.67/(RDM0*widthfraction)
        K0 = BaryonFraction*DMOMass*C0**3/(8*pi)
    else:
        C0 = CG
        K0 = KG

    alphas = linspace(0,1,NumAlpha)
    REXT = DMO[0]
    VEXT = []
    for i in range(len(REXT)):
        VEXT = VEXT+[Vext(C0,K0,REXT[i])]

    FG = DMO[1][0]
    #print(FG)
    VG = DMO[2][0]

    
    for j in range(len(alphas)):
        #print('Computing Initial Guess Alpha = ' + str(alphas[j]))
        NewSolution = SG.PhysicalState(FG,VG,N,0.8,alphas[j],VEXT,REXT,99.9)[0]
        FG = NewSolution[1][0]
        VG = NewSolution[2][0]

    
    haloinfo = HaloData(NewSolution,C0,K0)
    fdm = haloinfo[3]
    index = haloinfo[-1]
    if index >= len(fdm):
        HalfFraction = fdm[-1]
    else:
        HalfFraction = fdm[index]
    
    if HalfFraction < DesiredHalfFraction:
        while HalfFraction < DesiredHalfFraction:
            #print('Decreasing C ' + str(HalfFraction))
            C0 = C0*0.9
            K0 = BaryonFraction*DMOMass*C0**3/(8*pi)
            for i in range(len(REXT)):
                VEXT = VEXT+[Vext(C0,K0,REXT[i])]
            newparams = SG.PhysicalState(FG,VG,N,0.8,1.0,VEXT,REXT,99.9)[0]
            FG = newparams[1][0]
            VG = newparams[2][0]
            #newdata = KSolver(Fguess,Vguess,1.25,1,C0,K0)
            haloinfo = HaloData(newparams,C0,K0)
            fdm = haloinfo[3]
            index = haloinfo[-1]
            if index > len(fdm):
                HalfFraction = fdm[-1]
            else:
                HalfFraction = fdm[index]
            
    elif HalfFraction > DesiredHalfFraction:
        while HalfFraction > DesiredHalfFraction:
            #print('Increasing C ' + str(HalfFraction))
            C0 = C0*1.1
            K0 = BaryonFraction*DMOMass*C0**3/(8*pi)
            for i in range(len(REXT)):
                VEXT = VEXT+[Vext(C0,K0,REXT[i])]
            newparams = SG.PhysicalState(FG,VG,N,0.8,1.0,VEXT,REXT,99.9)[0]
            FG = newparams[1][0]
            VG = newparams[2][0]
            #newdata = KSolver(Fguess,Vguess,1.25,1,C0,K0)
            haloinfo = HaloData(newparams,C0,K0)
            fdm = haloinfo[3]
            index = haloinfo[-1]
            if index >= len(fdm):
                HalfFraction = fdm[-1]
            else:
                HalfFraction = fdm[index]
    return(newparams,C0,K0)
        
    



#the following will test the routine for the first 25 states.  computation time per state
#scales linearly with state, so overall computation time scales quadratically with the total number computed
'''
Fguess =0
Vguess =0
Cguess =0
Kguess =0
Solution = []
fraction = 0.111

for i in range(25):
    print(i)
    if i ==0:
        Data = ExponentialGalaxy(False,False,0,0.65,fraction,False,False)
    else: 
        Data = ExponentialGalaxy(Fguess,Vguess,i,0.65,fraction,Cguess,Kguess)
    Fguess = Data[0][1][0]
    Vguess = Data[0][2][0]
    Cguess = Data[1]
    Kguess = Data[2]
    Solution_i = [Fguess,Vguess,Cguess,Kguess]
    Solution = Solution+[Solution_i]


#savetxt('SolutionFile.txt',Solution)

plot(Data[0][0],Data[0][1])
show()
'''