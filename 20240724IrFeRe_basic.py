import os
import pylab as plt
import numpy as np
# import pandas as pd
import mass
import sys
from mass.off import ChannelGroup, Channel, getOffFileListFromOneFile
from mass.calibration import _highly_charged_ion_lines
import mass.calibration.hci_models
#import mass.calibration._highly_charged_ion_lines
import lmfit

filename = r"/Users/gmondeel/Documents/Mass/NeLikeData/20240724_off/0000/20240724_run0000_chan1.off"

plt.ion()
data = ChannelGroup(getOffFileListFromOneFile(filename, maxChans=999), verbose=True, channelClass = Channel, excludeStates=['IGNORE', 'STOP'])
data.setDefaultBinsize(1.)
mass.line_models.VALIDATE_BIN_SIZE = False 
ds = data.firstGoodChannel()

#Make a dictionary with lists of aliases for each element.
#Pass this anywhere with the 'states=None' argument, e.g., data.plotHist(..., states=statesDict["W"])
statesDict = {
    "Re_OFF" : ["E_OFF", "H_OFF"],
    "Ir_ON" : ["I_ON"],
    "Ir_OFF" : ["I_OFF"],
    "Cal": ["G_OFF"],
    "CalOn": ["G_OFF", "I_ON"], #states where the calibration source was turned on
    "SciOrCal":["E_OFF", "H_OFF", "I_OFF", "G_OFF"]
}


###view filtValue plot
#data[1].plotHist(np.arange(0,55000,10), "filtValue", coAddStates=False, states=statesDict["CalOn"])

def getEnergy(lineName):
    return mass.getmodel(lineName).spect.nominal_peak_energy
def getPH(lineName): #See lines with mass.spectra and mass.STANDARD_FEATURES
    energy=getEnergy(lineName)
    return ds.recipes['energy'].f.energy2ph(energy)

ds.calibrationPlanInit("filtValue")
calStates = statesDict["Cal"]

ds.calibrationPlanAddPoint(31300, "MnKAlpha", states=calStates)
ds.calibrationPlanAddPoint(43120, "ZnKAlpha", states=calStates)

data.alignToReferenceChannel(referenceChannel=ds,
                            binEdges=np.arange(3000, 60000, 3), attr="filtValue")#, _rethrow=True)
data.cutAdd("cutForLearnDC", lambda energyRough: np.logical_and(energyRough > 8000, energyRough < 11000), setDefault=False)
data.learnPhaseCorrection(indicatorName="filtPhase", uncorrectedName="filtValue", correctedName = "filtValuePC", states=statesDict["SciOrCal"])
data.learnDriftCorrection(uncorrectedName="filtValuePC", indicatorName="pretriggerMean", correctedName="filtValueDC",
                            states=statesDict["Cal"], cutRecipeName="cutForLearnDC")#, _rethrow=True)

data.calibrateFollowingPlan("filtValueDC", calibratedName="energy", dlo=70, dhi=70, approximate=False, overwriteRecipe=True)
data[6].markBad("bad")

data.linefit('ZnKAlpha', states=statesDict['Cal'])