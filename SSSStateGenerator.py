from matplotlib.pylab import *
from matplotlib.pyplot import plot





Upsilon = 100 #This is the computational scale of the boson mass
dr =.02 #This is the length resolution of the ODE solver

#The following will give you the values of various scales in comparison to computational units
#to unitize any quantities, multiply the computed value by these scale factors.  
#Mass =1E-22#quote mass in eV/c^2
Mass = 1E-22
Wavenumber = Mass / (1.973E-7)  #gives wavenumber in inverse meter (M/hbarc where hbarc =197.3evnm)
Wavenumber = Wavenumber * (9.461E15) #converts to inverse light years
Conversion = Wavenumber / Upsilon #Conversion factor for length scale
Lengthscale = 1.0/Conversion  #length scale in Ly
Lengthscale2 = Lengthscale / 3262.  #length scale in kpc
Lengthscale3 = Lengthscale*9.461E15 #length in meter
Timescale = Lengthscale #timescale in years
Timescale2= 3.154E7*Timescale #seconds
Massscale = Lengthscale*(6.41E12) #In solar masses
Massscale2 = Massscale*1.988E30 #in kg
DensityScale = Massscale2/(Lengthscale3)**3 
DensityScale2 = Massscale/(Lengthscale2)**3
bigG = 6.674E-11 #Newton constant in SI units


#Phi defined for mathematical convenience
def Phi(m,r):
    if r==0:
        P = 1.
    else:
        P = 1.-2*m/r
    return P
def Phir(m,r,mr):
    if r == 0:
        Pr = 0
    else:
        Pr = -2.*((r*mr-m)/(r*r))
    return Pr
    
    
##The following 3 definitions are for the Einstein Klein Gordon system in spherical symmetry    
def x(r,f,fr,P,m,v,ALPHA,Vext,w): #Derivative of M
    VV = v+ALPHA*Vext
    X = 4.*pi*(r*r/(Upsilon*Upsilon))*((Upsilon*Upsilon + w*w*exp(-2.*VV))*f*f + P*fr*fr)
    return X
def y(r,f,fr,P,m,v,ALPHA,Vext,w):#Derivative of V
    
    VV = v+ALPHA*Vext
    if r == 0:
        Y = (0 - 4.*pi*(r/(Upsilon*Upsilon))*((Upsilon*Upsilon-w*w*exp(-2.*VV))*f*f - P*fr*fr))/P
    else:    
        Y = (m/(r*r) - 4.*pi*(r/(Upsilon*Upsilon))*((Upsilon*Upsilon-w*w*exp(-2.*VV))*f*f - P*fr*fr))/P
    return Y
    
def z(r,f,fr,P,Pr,v,vr,ALPHA,Vext,w): #Second Derivative of F
    VV = v+ALPHA*Vext
    if r == 0:
        Z = (Upsilon*Upsilon-w*w*exp(-2.*VV))*f/P 
    else:
        Z = (Upsilon*Upsilon-w*w*exp(-2.*VV))*f/P -2.*fr/r - vr*fr - Pr*fr/(2.*P)
    return Z
    
    
#The following are for the Runge Kutta Solving Routine, specialized to the ODEs 
def VectorPrime(Vector,r,ALPHA,Vext,w): #Takes an Input format Vector = [M,V,F,Fr]
    m = Vector[0]
    v = Vector[1]
    f = Vector[2]
    fr = Vector[3]
    P = Phi(m,r)
    MR = x(r,f,fr,P,m,v,ALPHA,Vext,w)
    VR = y(r,f,fr,P,m,v,ALPHA,Vext,w)
    Pr = Phir(m,r,MR)
    FRR = z(r,f,fr,P,Pr,v,VR,ALPHA,Vext,w)
    return array((MR,VR,fr,FRR))  #Returns the Derivative vector
    

#Helper function that counts how many zeros a function has
def crossings(F):
    n = 0 
    for j in range(len(F)-1):
        if sign(F[j+1]) != sign(F[j]):
            n = n+1
        else:
            n = n
    return n


#This is the RK 4th order routine - It is built to handle the inclusion of an external potential
#The external potential is included via continuation parameter, V_tot = V+ alpha*V_ext . 
#Alpha = 0 corresponds to no inclusion, Alpha = 1 corresponds to full inclusion
#External Potential VEXT must be specified as a vector corresponding to a vector of radii, REXT 
# -- If alpha=0, one still needs to insert placeholder values for VEXT and REXT
#If w>Upsilon then this routine will not converge, only use for w<Upsilon.
def RK4Solver(F0,V0,Scale,ALPHA,VEXT,REXT,w):#Takes F0, V0 as the initial conditions, and solves the ODEs for a specified frequency 
   
    dR = REXT[1]-REXT[0] #for the external data
    #These specify the initial conditions and some helpful quantities
    M = [0]
    V = [V0]
    F = [F0]
    Fr = [0]
    R = [0]
    k = [ Upsilon*Upsilon-w*w*exp(-2*V[0])]
    Density = [2*F0*F0]
    A = [ sqrt(F[0]**2 + Fr[0]**2/(abs(k[0])))  ]
    Lambda = [ 2*pi/ sqrt(abs(-k[0]/Phi(M[0],R[0])))]
    Vrot = [0]
        
    SlopeCondition = [-F[0]]
    
    #root = [4*k[0]]
    n = 0
    while k[-1]<0: #this solves the ODE until exponential behavior occurs.  Importantly, if the solution oscillates infinitely, this will keep solving forever - make sure to check that w>
                   # in that case, one should consider a different routine.  
        
        #following if statements make sure the correct value of the 
        #external potential is used for the solver... 
        if ALPHA ==0:
            index = 0
        else:
            temp = int(R[-1]/dR) #this finds the position in the array
            if temp >(len(VEXT)-1):
                index = len(VEXT)-1
            else:
                index = temp
                
                
        #The followinUpsiloUpsilon is the actual RK4 forward step        
        Vector = array((M[-1],V[-1],F[-1],Fr[-1]))
        K1 = VectorPrime(Vector,R[-1],ALPHA,VEXT[index],w)
        Vector2 = array(K1/2*dr)+Vector
        K2 = VectorPrime(Vector2,R[-1]+dr/2,ALPHA,VEXT[index],w)
        Vector3 = array(K2/2*dr) + Vector
        K3 = VectorPrime(Vector3,R[-1]+dr/2,ALPHA,VEXT[index],w)
        Vector4 = array(K3*dr) + Vector
        K4 = VectorPrime(Vector4,R[-1]+dr,ALPHA,VEXT[index],w)
        NewVector = Vector + array(1./6*K1 + 1./3*K2 + 1./3*K3 + 1./6*K4)*dr
        M = M + [NewVector[0]]
        V = V + [NewVector[1]]
        F = F + [NewVector[2]]
        Fr = Fr + [NewVector[3]]
        R = R + [R[-1]+dr]
        #useful quantities computed
        k = k + [Upsilon*Upsilon-w*w*exp(-2.*V[-1])]
        A = A + [sqrt(F[-1]**2 + Fr[-1]**2/(abs(k[-1])))]   
        Lambda = Lambda + [2*pi/sqrt(abs(-k[-1]/Phi(M[-1],R[-1])))]
        SlopeCondition = SlopeCondition +[-Fr[-1]/(1.0/R[-1]+sqrt(Upsilon*Upsilon-w*w))-F[-1]]
        Vrot = Vrot + [sqrt(M[-1]/R[-1])]
        
        #print 'radius ' + str(R[-1]) + ' k value ' + str(k[-1])
    RDM = R[-1] #Records the k=0 radius
    
    #continues the RK4 method to a given radius past RDM 
    while R[-1]<(Scale+RDM):#Continues solving past the point where exponential behavior starts, up to a tolerance Scale
        if ALPHA ==0:
            index = 0
        else:
            temp = int(R[-1]/dR)
            if temp >(len(VEXT)-1):
                index = len(VEXT)-1
            else:
                index = temp
        Vector = array((M[-1],V[-1],F[-1],Fr[-1]))
        K1 = VectorPrime(Vector,R[-1],ALPHA,VEXT[index],w)
        Vector2 = array(K1/2*dr)+Vector
        K2 = VectorPrime(Vector2,R[-1]+dr/2,ALPHA,VEXT[index],w)
        Vector3 = array(K2/2*dr) + Vector
        K3 = VectorPrime(Vector3,R[-1]+dr/2,ALPHA,VEXT[index],w)
        Vector4 = array(K3*dr) + Vector
        K4 = VectorPrime(Vector4,R[-1]+dr,ALPHA,VEXT[index],w)
        NewVector = Vector + array(1./6*K1 + 1./3*K2 + 1./3*K3 + 1./6*K4)*dr
        M = M + [NewVector[0]]
        V = V + [NewVector[1]]
        F = F + [NewVector[2]]
        Fr = Fr + [NewVector[3]]
        R = R + [R[-1]+dr]
        k = k + [Upsilon*Upsilon-w*w*exp(-2*V[-1])]
        A = A + [sqrt(F[-1]**2 + Fr[-1]**2/(abs(k[-1])))]   
        Lambda = Lambda + [2*pi/sqrt(abs(-k[-1]/Phi(M[-1],R[-1])))]
        SlopeCondition = SlopeCondition +[-Fr[-1]/(1.0/R[-1]+sqrt(Upsilon*Upsilon-w*w))-F[-1]]
        Vrot = Vrot + [sqrt(M[-1]/R[-1])]
    N = crossings(F)
    Ns = crossings(SlopeCondition)
    Mtotal = M[-1]
    if R[-1]>0:
        Vinf = V[-1]+Mtotal/R[-1]
    else:
        Vinf = 0
    
    return (R,F,V,M,N,Ns,Fr,SlopeCondition,RDM,Mtotal,Vrot,Lambda,A,Vinf)



def ResultScaling(R,F,V,M,C): #ReScales the entire solution given a constant C.  
    R2 = 1/C*R
    F2 = C*C*F
    V2 = C*C*V
    M2 = C*M
    return(R2,F2,V2,M2)


#This is another version of the RK4 solver that is built to stop after a certain number
#of oscillations of the function.  After NSTOP zeros occur, the solver will stop
def RK4Solver2(F0,V0,Scale,ALPHA,VEXT,REXT,w,NSTOP):
    dR = REXT[1]-REXT[0]
    M = [0]
    V = [V0]
    F = [F0]
    Fr = [0]
    R = [0]
    k = [ Upsilon*Upsilon-w*w*exp(-2*V[0])]
    Density = [2*F0*F0]
    A = [ sqrt(F[0]**2 + Fr[0]**2/(abs(k[0])))  ]
    Lambda = [ 2*pi/ sqrt(abs(-k[0]/Phi(M[0],R[0])))]
    Vrot = [0]
        
    SlopeCondition = [-F[0]]
    zeros = 0 
    root = [4*k[0]]
    n = 0
    while k[-1]<0 and zeros<NSTOP:
        if ALPHA ==0:
            index = 0
        else:
            temp = int(R[-1]/dR)
            if temp >(len(VEXT)-1):
                index = len(VEXT)-1
            else:
                index = temp
        Vector = array((M[-1],V[-1],F[-1],Fr[-1]))
        K1 = VectorPrime(Vector,R[-1],ALPHA,VEXT[index],w)
        Vector2 = array(K1/2*dr)+Vector
        K2 = VectorPrime(Vector2,R[-1]+dr/2,ALPHA,VEXT[index],w)
        Vector3 = array(K2/2*dr) + Vector
        K3 = VectorPrime(Vector3,R[-1]+dr/2,ALPHA,VEXT[index],w)
        Vector4 = array(K3*dr) + Vector
        K4 = VectorPrime(Vector4,R[-1]+dr,ALPHA,VEXT[index],w)
        NewVector = Vector + array(1./6*K1 + 1./3*K2 + 1./3*K3 + 1./6*K4)*dr
        M = M + [NewVector[0]]
        V = V + [NewVector[1]]
        F = F + [NewVector[2]]
        if sign(F[-1])!=sign(F[-2]):
            zeros = zeros+1
        Fr = Fr + [NewVector[3]]
        R = R + [R[-1]+dr]
        k = k + [Upsilon*Upsilon-w*w*exp(-2.*V[-1])]
        A = A + [sqrt(F[-1]**2 + Fr[-1]**2/(abs(k[-1])))]   
        Lambda = Lambda + [2*pi/sqrt(abs(-k[-1]/Phi(M[-1],R[-1])))]
        SlopeCondition = SlopeCondition +[-Fr[-1]/(1.0/R[-1]+sqrt(Upsilon*Upsilon-w*w))-F[-1]]
        Vrot = Vrot + [sqrt(M[-1]/R[-1])]
        
        #print 'radius ' + str(R[-1]) + ' k value ' + str(k[-1])
    RDM = R[-1]
    
  
     
    
    while R[-1]<(Scale+RDM) and zeros <NSTOP:
        if ALPHA ==0:
            index = 0
        else:
            temp = int(R[-1]/dR)
            if temp >(len(VEXT)-1):
                index = len(VEXT)-1
            else:
                index = temp
        Vector = array((M[-1],V[-1],F[-1],Fr[-1]))
        K1 = VectorPrime(Vector,R[-1],ALPHA,VEXT[index],w)
        Vector2 = array(K1/2*dr)+Vector
        K2 = VectorPrime(Vector2,R[-1]+dr/2,ALPHA,VEXT[index],w)
        Vector3 = array(K2/2*dr) + Vector
        K3 = VectorPrime(Vector3,R[-1]+dr/2,ALPHA,VEXT[index],w)
        Vector4 = array(K3*dr) + Vector
        K4 = VectorPrime(Vector4,R[-1]+dr,ALPHA,VEXT[index],w)
        NewVector = Vector + array(1./6*K1 + 1./3*K2 + 1./3*K3 + 1./6*K4)*dr
        M = M + [NewVector[0]]
        V = V + [NewVector[1]]
        F = F + [NewVector[2]]
        if sign(F[-1])!=sign(F[-2]):
            zeros = zeros+1
        Fr = Fr + [NewVector[3]]
        R = R + [R[-1]+dr]
        k = k + [Upsilon*Upsilon-w*w*exp(-2*V[-1])]
        A = A + [sqrt(F[-1]**2 + Fr[-1]**2/(abs(k[-1])))]   
        Lambda = Lambda + [2*pi/sqrt(abs(-k[-1]/Phi(M[-1],R[-1])))]
        SlopeCondition = SlopeCondition +[-Fr[-1]/(1.0/R[-1]+sqrt(Upsilon*Upsilon-w*w))-F[-1]]
        Vrot = Vrot + [sqrt(M[-1]/R[-1])]
    N = crossings(F)
    Ns = crossings(SlopeCondition)
    Mtotal = M[-1]
    if R[-1]>0:

        Vinf = V[-1]+Mtotal/R[-1]
    else:
        Vinf = 0
    return (R,F,V,M,N,Ns,Fr,SlopeCondition,RDM,Mtotal,Vrot,Lambda,A,Vinf)



#Following returns the actual values of F0 and V0
def Nfind(F0,V0,N,Scale,ALPHA,VEXT,REXT,w): #finds a solution of order N to a Scale past RDM.  Requires a reasonable guess of F0, V0.  
    data0 = RK4Solver2(F0,V0,0,ALPHA,VEXT,REXT,w,N+2)
    Ns = data0[5] 
    
    sc = 0   
    #The following if statements contain a section search for solutions of order N
    factor = .01
    if Ns<N:
        while Ns!=N:
            #print(str(Ns-N))
            Ftest = F0*(1.0-factor) #Decreasing F0 results in more oscillations
            datatest = RK4Solver2(Ftest,V0,0,ALPHA,VEXT,REXT,w,N+2)
            Ns = datatest[5]
            if Ns>N: #the factor was too large 
                factor = factor*3.0/4.0
            elif Ns<N:#the factor was too small
                factor = factor*2.0 
        F0 = Ftest
        
    elif Ns>N:
        while Ns!=N:
            #print(str(Ns-N))
            Ftest = F0*(1.0+factor) #Increasing F0 results in less oscillations
            datatest = RK4Solver2(Ftest,V0,0,ALPHA,VEXT,REXT,w,N+2)
            Ns = datatest[5]
            if Ns<N:
                factor = factor*3.0/4.0
            elif Ns>N:
                factor = factor*2.0
        F0 = Ftest
    
    #Now the solution is extended out to the desired scale via continuation, 
    while sc < Scale:
        
        sc = sc+5*dr #right now extending 5 steps at a time
        data = RK4Solver2(F0,V0,sc,ALPHA,VEXT,REXT,w,N+2)
        Ns = data[5]
        #following routing could probably be improved but works
        while Ns!= N:
            
            if Ns<N:
                F0 = F0*(1-.01/(sc+1))   
            elif Ns>N:
                F0 = F0*(1+.01/(sc+1))
            data = RK4Solver2(F0,V0,sc,ALPHA,VEXT,REXT,w,N+2)
            Ns=data[5]
    return (F0,V0) 


#This takes a good guess for the bound state, and uses bisection to find the 
#exact values for the initial conditions   
def BoundFind(F0,V0,N,Scale,ALPHA,VEXT,REXT,w):
    #print('finding left')
    dataleft = Nfind(F0,V0,N,Scale,ALPHA,VEXT,REXT,w) #get a state of order N
    #print('finding right')
    dataright = Nfind(dataleft[0],dataleft[1],N+1,Scale,ALPHA,VEXT,REXT,w)#get a state of order N+1
    #print('found')
    
    #perform a bisection search in F0
    Fleft = dataleft[0]
    Fright = dataright[0]
    #print('Bisecting')
    t = 0 
    #bisects until the answers differ by no more than 0.01%
    Fguess = (Fleft + Fright)/2.0
    PercentDiff = abs((Fguess-Fleft)/Fleft)*100
    while PercentDiff > 0.01:
       # t = t+1
        Fguess = (Fleft + Fright)/2.0
        newdata = RK4Solver(Fguess,V0,Scale,ALPHA,VEXT,REXT,w)
        Ns = newdata[5]
        if Ns==N:
            Fleft = Fguess
        elif Ns==(N+1):
            Fright = Fguess
        PercentDiff = abs((Fguess-Fleft)/Fleft)*100
            
    finaldata = RK4Solver(Fleft,V0,Scale,ALPHA,VEXT,REXT,w)
    #return (Fleft,V0,finaldata[-1])
    return (finaldata,Fleft,V0,finaldata[-1])

def PhysicalState(F0,V0,N,Scale,ALPHA,VEXT,REXT,w):

    data = BoundFind(F0,V0,N,Scale,ALPHA,VEXT,REXT,w)
    F0 = data[1]
    #print(F0)
    
    Vinf = data[-1]
   
    
    
    if Vinf>0:
        while Vinf>0:
            #print(Vinf)
            V0 = 0.99*V0
            data = BoundFind(F0,V0,N,Scale,ALPHA,VEXT,REXT,w)
            F0 = data[1]
            Vinf = data[-1]


                    
    elif Vinf<0:
        while Vinf<.01*V0:
            V0 = 1.01*V0
            data = BoundFind(F0,V0,N,Scale,ALPHA,VEXT,REXT,w)
            F0 = data[1]
            PercentDiff = 100*abs((Vinf-data[-1])/Vinf)
            Vinf = data[-1]
            if PercentDiff<.01:
                break
    if Vinf<0:
        while Vinf<.01*V0:
            V0 = 1.01*V0
            data = BoundFind(F0,V0,N,Scale,ALPHA,VEXT,REXT,w)
            F0 = data[1]
            PercentDiff = 100*abs((Vinf-data[-1])/Vinf)
            Vinf = data[-1]
            if PercentDiff<.01:
                break
    #print('Vinf is' + str(Vinf))     
    return data
                   




#this routine approximates the PhysicalState functio in the low-field non-relativistic regime by 
#taking advantage of scaling relationships.  This is adapted to easily produce an excited state with a defined total mass, and V_infinity =0
def DMOExcitedState(N,Mtot,Steps,Scaling,FG,VG):
    #if Scaling = False, then this will just return the DMO state for w=99.9 and the Mtot value will not be scaled
    #if Scaling = True then this will scale the solution to have total mass of Mtot

     #returns the profile of an Nth excited state with central density Rho.  (R,Psi,V,M)
     #solution will be computed for Steps number of distance steps passed the decay radius.
      
    #These guess are chosen specifically to be near the w=99.9 solution for N=0.  They SHOULD NOT BE NAIVELY CHANGED at risk of breaking the routine
    if FG == False:
        Fguess = .028839362371243148
        Vguess = -.001940598
    else:
        Fguess = FG
        Vguess = VG
    Vext = [1,2]  # this is just a place holder for the routine since it requires one to specify external potentials
    Rext = [1,2]  # I am using alpha = 0 so no actual inclusion of Vext or Rext happens

    #computes a bound state order N
    Rmax = Steps*dr
    data0 = BoundFind(Fguess,Vguess,N,Rmax,0,Vext,Rext,99.9)
    #shifts the potential so that vinf = 0 convention is satisfied
    Vinf = data0[-1]
    w2 = 99.9*exp(-Vinf)
    Finit = data0[1]
    Vinit = data0[2]-Vinf
    #re solves to ensure the solution works.  
    data1 = BoundFind(Finit,Vinit,N,Rmax,0,Vext,Rext,w2)
    #Now rescales the solution size to the desired mass
    factor = Mtot/data1[0][-5]
    Finaldata = ResultScaling(array(data1[0][0]),array(data1[0][1]),array(data1[0][2]),array(data1[0][3]),factor)
    w2 = -factor*factor*(Upsilon-w2)+Upsilon
    #returns the wavefunction and frequency
    if Scaling == True:
        output = (Finaldata,w2)
    else:
        output = (data1[0],99.9)
    return output


#The below are testing routines 
'''
#the following blocks should produce a state with 4 crossings
Ft = .028839362371243148
Vt = -.001940598
Vext = [1,2]  #placeholder value, alpha=0
Rext = [1,2]  #placeholder value alpha=0
Test = PhysicalState(Ft,Vt,4,1.2,0,Rext,Vext,99.9)
plot(Test[0][0],Test[0][1])
show()


freq = []
for i in range(5):
    Test = DMOExcitedState(i,.01,0,True,False,False)
    plot(Test[0][0],Test[0][1])
    freq = freq+[Test[1]]
show()'''

