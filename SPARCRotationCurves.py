
#enable the below import statement if you wish to render plots in the script
#from matplotlib.pylab import plot, show

import pandas as pd


#This is the original data file from SPARC Survey at http://astroweb.cwru.edu/SPARC/MassModels_Lelli2016c.mrt
#Please cite SPARC appropriately if you use this script
Data = pd.read_csv("RotationCurveDataSPARC.txt")
#This is how the data file is formated
#Position 0-10 = Name of Galaxy
#Position 12-17 = Distance in Mpc
#Position 19-24 = Galaxy Radius in kpc
#Position 26-31 = Circular velocity observed
#Position 33-37 = Uncertainty in circular velocity
#Position 39-44 = Gas Velocity Contribution 
#Position 46-51 = Disk Velocity Contribution
#Position 53-58 = Bulge Velocity Contribution
#Position 60-66 = Disk Surface Brightness
#Position 68-75 = Bulge Surface Brightness

#Creating a DataFrame to store the data, columns defined as above
FramedData = pd.DataFrame(columns = ['Name', 'Distance', 'Radius','Vcirc','SigmaV','Vgas','Vdisk','Vbulge','DiskBrightness','BulgeBrightness'],
index=range(len(Data)))

#Choosing the correct data for each column
for i in range(len(Data)):
    FramedData.iloc[i][0] = str(Data.iloc[i][0][0:11])
    FramedData.iloc[i][1] = float(Data.iloc[i][0][12:18])
    FramedData.iloc[i][2] = float(Data.iloc[i][0][19:25])
    FramedData.iloc[i][3] = float(Data.iloc[i][0][26:32])
    FramedData.iloc[i][4] = float(Data.iloc[i][0][33:38])
    FramedData.iloc[i][5] = float(Data.iloc[i][0][39:45])
    FramedData.iloc[i][6] = float(Data.iloc[i][0][46:52])
    FramedData.iloc[i][7] = float(Data.iloc[i][0][53:59])
    FramedData.iloc[i][8] = float(Data.iloc[i][0][60:67])
    FramedData.iloc[i][9] = float(Data.iloc[i][0][68:76])

#Unique galaxies in the dataframe
Names = FramedData['Name'].unique()

#Selects the data for entries with Name of galaxy
def GalaxySelect(galaxy):
    galaxy = FramedData.query('Name == @galaxy')
    return galaxy
#Finds the Mass Components for galaxies named galaxy
def MassModel(galaxy):
    Galaxy = GalaxySelect(galaxy)
    Name = Galaxy['Name']
    Radii = Galaxy['Radius']
    Velocity = Galaxy['Vcirc']
    Sigma = Galaxy['SigmaV']
    Gas  = Galaxy['Vgas']
    Bulge = Galaxy['Vbulge']
    Disk = Galaxy['Vdisk']
    return(Name, Radii, Velocity, Sigma, Gas, Bulge, Disk)


#Uncomment below to test the script by plotting a galaxy rotation curve
Filtered = FramedData.query('Vbulge>Vdisk')
Names2 = Filtered['Name'].unique()

for i in range(len(Names2)):
    print(i)
    data = MassModel(Names2[i])
    plot(data[1].iloc[:],data[2].iloc[:])
show()

'''
TestGal = MassModel(Names[0])
plot(TestGal[1].iloc[:],TestGal[2].iloc[:])
show()
'''


