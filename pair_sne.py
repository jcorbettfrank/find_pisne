#%%
from lightCurve import lightCurve
import pandas as pd 
import numpy as np
import os
class pairSNE(lightCurve):

    def __init__(self, name):
        lightCurve.__init__(self, name)
    def loadData(self, dir):
        fn = [f for f in os.listdir(dir) 
                if f[-4:]=='spec']
        numT = len(fn)
        if numT == 0:
            self.rawDat = None
            self.tKey = None
            self.lambdaKey = None
            return None 
        self.tKey = [float(f.split("_")[1][1:]) for f in fn]
        fn = [os.path.join(dir,f) for f in fn]
        #lets read first file to see how many bins we have in spectra
        arr = np.loadtxt(fn[0],skiprows=1,usecols=[0,1,4])
        numBins = arr.shape[0]
        #3D array will hold all the data; Each 2d plane holds two cols
        #specific luminosity (erg/s/angstrom) and number of photons of that wavelength
        #rows are labeled by wavelength.
        #Each plane corresponds to a  day. We will store the days and wavelengths once
        #since they are the same for each day
        holder = np.zeros((numBins, 2, numT))
        self.lambdaKey = arr[:,0] #angstroms, rows,
        #have tKey  which is in days and correponds to planes 
         
        #now load in the actual data
        for idx, f in enumerate(fn):
            arr = np.loadtxt(f, skiprows=1, usecols=[0,1,4])
            holder[:,:,idx] = arr[:,1:]
        self.rawDat = holder
         
#%%
if __name__=='__main__':
    data_dir = os.path.join(os.path.join(os.getcwd(), "data"), 'B200') 
    pair = pairSNE("red_p")
    pair.loadData(data_dir)
#%%
if __name__=='__main__':
    print(pair.rawDat)
    import matplotlib.pyplot as plt
    plt.plot(pair.rawDat[:,0,100])
    print("name is {}".format(pair.name))

# %%
