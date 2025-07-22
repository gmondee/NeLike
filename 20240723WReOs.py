import os
import pylab as plt
import numpy as np
#import pandas as pd
import mass
import sys
from mass.off import ChannelGroup, Channel, getOffFileListFromOneFile
from mass.calibration import _highly_charged_ion_lines
import mass.calibration.hci_models
#import mass.calibration._highly_charged_ion_lines
import lmfit

import CustomFluorescenceLines

#Laptop
#filename = r"C:\Users\Grant Mondeel\Box\CfA\TES\Ne-Like\20240723_off\0000\20240723_run0000_chan1.off"

#PC
#filename = r"C:\Users\lamat\Box\CfA\TES\Ne-Like\20240723_off\0000\20240723_run0000_chan1.off"

#Mac
filename = r"/Users/gmondeel/Documents/Mass/NeLikeData/20240723_off/0000/20240723_run0000_chan1.off"

plt.ion()
data = ChannelGroup(getOffFileListFromOneFile(filename, maxChans=999), verbose=True, channelClass = Channel, excludeStates=['A', 'IGNORE', 'STOP'])
data.setDefaultBinsize(1.)
mass.line_models.VALIDATE_BIN_SIZE = False ###with Galen's cal, had to go down to 0.025eV bins...
ds = data.firstGoodChannel()
#data.setOutputDir(baseDir=baseOutputDirectory, deleteAndRecreate=False)

###Some residual based cuts
#data.learnResidualStdDevCut() #see if this is works

###view some pulses
# off = mass.off.OffFile(filename)
# x,y = off.recordXY(1)

# plt.figure()
# plt.plot(off.basis[:,::-1]) # put the higher signal to noise components in the front
# plt.xlabel("sample number")
# plt.ylabel("weight")
# plt.legend(["pulseMean","derivativeLike","pulseLike","extra0","extra1"])
# plt.title("basis components")

# plt.figure()
# plt.plot(x,y)
# plt.xlabel("time (s)")
# plt.ylabel("reconstructed signal (arbs)")
# plt.title("example reconstructed pulse")

###Optional: assignes names to the states based on what was being observed
# data.experimentStateFile.aliasState("START", "Cal")
# data.experimentStateFile.aliasState("B", "W 1")
# data.experimentStateFile.aliasState("D", "W 2")
# data.experimentStateFile.aliasState("E", "Re 1")
# data.experimentStateFile.aliasState("F", "Re 2")
# data.experimentStateFile.aliasState("G", "Os 1")


#Make a dictionary with lists of aliases for each element.
#Pass this anywhere with the 'states=None' argument, e.g., data.plotHist(..., states=statesDict["W"])
statesDict = {
    "W_ON"  : ["D_ON"],
    "W_OFF"  : ["B_OFF", "D_OFF"],
    "Re_ON" : ["E_ON", "F_ON"],
    "Re_OFF" : ["E_OFF", "F_OFF"],
    "Os_ON" : ["G_ON"],
    "Os_OFF" : ["G_OFF"],
    "Cal": ["E_ON"],#"START_OFF"],
    "CalOn": ["D_ON", "E_ON", "F_ON", "G_ON"], #states where the calibration source was turned on except START because of huge drift
    "SciOrCal":["B_OFF", "D_OFF", "E_OFF", "F_OFF", "G_OFF", "D_ON", "E_ON", "F_ON", "G_ON"]
}


###view filtValue plot
data[1].plotHist(np.arange(0,55000,10), "filtValue", coAddStates=False, states=statesDict["CalOn"])

###Check how many pulses are in the data file
# totalPulses = 0
# for chan in data.keys():
#     print("channel", chan,"has",len(data[chan].getAttr('filtValue','I')))
#     totalPulses += len(data[chan].getAttr('filtValue', 'I'))
# print("Combining all channels gives",totalPulses,"pulses.")

def getEnergy(lineName):
    return mass.getmodel(lineName).spect.nominal_peak_energy
def getPH(lineName): #See lines with mass.spectra and mass.STANDARD_FEATURES
    energy=getEnergy(lineName)
    return ds.recipes['energy'].f.energy2ph(energy)
#Sc_keys = [k for k in mass.STANDARD_FEATURES.keys() if 'Sc' in k]

ds.calibrationPlanInit("filtValue")
calStates = statesDict["Cal"]

### Identify calibration lines
### W lines from beiersdorfer 2012

#too low e
# # ds.calibrationPlanAddPoint(7034, "MgKAlpha", states=calStates)
# ds.calibrationPlanAddPoint(8724, "AlKAlpha", states=calStates)
# ds.calibrationPlanAddPoint(14863, "ClKAlpha", states=calStates)
# # ds.calibrationPlanAddPoint(13267, "SKAlpha", states=calStates)
# ds.calibrationPlanAddPoint(22563, "ScKAlpha", states=calStates)
# ds.calibrationPlanAddPoint(24414, "ScKBeta", states=calStates)
# ds.calibrationPlanAddPoint(26826, "VKAlpha", states=calStates)
# # ds.calibrationPlanAddPoint(29045, "VKBeta", states=calStates) #Blended with CrKAlpha


ds.calibrationPlanAddPoint(31263, "MnKAlpha", states=calStates)
ds.calibrationPlanAddPoint(33600, "FeKAlpha", states=calStates)
ds.calibrationPlanAddPoint(35955, "CoKAlpha", states=calStates)
# ds.calibrationPlanAddPoint(36515, "FeKBeta", states=calStates)#low counts
# ds.calibrationPlanAddPoint(39106, "CoKBeta", states=calStates)#low counts
ds.calibrationPlanAddPoint(40767, "CuKAlpha", states=calStates)

ds.calibrationPlanAddPoint(41840, "WM2+3G_23", states=statesDict["W_OFF"])
# ds.calibrationPlanAddPoint(41820, "WM2", states=statesDict["W_OFF"]) # 8299.22 eV, might not be not enough counts
# ds.calibrationPlanAddPoint(41858, "W3G", states=statesDict["W_OFF"]) # 8307.51 eV, has more counts but might get confused by M2 during fit
ds.calibrationPlanAddPoint(43193, "ZnKAlpha", states=calStates)
ds.calibrationPlanAddPoint(45173, "W3D", states=statesDict["W_OFF"]) # 9126.25 eV, only a few counts per channel but should be isolated
# # ds.calibrationPlanAddPoint(44335, "CuKBeta", states=calStates)
ds.calibrationPlanAddPoint(46908, "ZnKBeta", states=calStates)
ds.calibrationPlanAddPoint(48105, "GeKAlphaCustom", states=calStates)
ds.calibrationPlanAddPoint(52207, "GeKBetaCustom", states=calStates)






### Check calibration on just one channel
# ds.calibrateFollowingPlan("filtValue", calibratedName="energy", dlo=50, dhi=50, approximate=True, overwriteRecipe=True)
# ds.diagnoseCalibration()
# ds.plotHist(np.arange(800, 12000, 1), "energy", states=["Cal"], coAddStates=False)

data.alignToReferenceChannel(referenceChannel=ds,
                            binEdges=np.arange(6000, 60000, 5), attr="filtValue")#, _rethrow=True)
data.cutAdd("cutForLearnDC", lambda energyRough: np.logical_and(
energyRough > 8000, energyRough < 10000), setDefault=False)
data.cutAdd("cutForLearnDC_W", lambda energyRough: np.logical_and(
energyRough > 9000, energyRough < 9300), setDefault=False)

### Mass corrections.
# data.learnPhaseCorrection(indicatorName="filtPhase", uncorrectedName="filtValue", correctedName = "filtValuePC", states="Cal")#, cutRecipeName="cutForPC")
data.learnPhaseCorrection(indicatorName="filtPhase", uncorrectedName="filtValue", correctedName = "filtValuePC", states=statesDict["SciOrCal"])
data.learnDriftCorrection(uncorrectedName="filtValuePC", indicatorName="pretriggerMean", correctedName="filtValueDC",
                            states=statesDict["CalOn"], cutRecipeName="cutForLearnDC")#, _rethrow=True)
# data.learnDriftCorrection(uncorrectedName="filtValueDC", indicatorName="pretriggerMean", correctedName="filtValueDC_W",
#                             states=statesDict["W_OFF"], cutRecipeName="cutForLearnDC_W")#, _rethrow=True)
data.calibrateFollowingPlan("filtValueDC", calibratedName="energy", dlo=30, dhi=40, approximate=True, overwriteRecipe=True)
#data.qualityCheckLinefit("ZnKAlpha", positionToleranceFitSigma=3, worstAllowedFWHM=10, states=statesDict["Cal"], dlo=70, dhi=40)
ds.diagnoseCalibration()
data[6].markBad("bad")
# data.qualityCheckLinefit("GeKBeta", positionToleranceFitSigma=2, worstAllowedFWHM=12, states="Cal", dlo=50, dhi=40)

# calLines = ["AlKAlpha","ScKAlpha","VKAlpha","MnKAlpha",
#             "FeKAlpha","CoKAlpha","CuKAlpha"]#,"ZnKAlpha","GeKAlphaCustom"]

# calLines = ["ClKAlpha", "FeKAlpha","ZnKAlpha","GeKAlphaCustom"]
# for line in calLines:
#     data.qualityCheckLinefit(line, positionToleranceFitSigma=5, worstAllowedFWHM=15, states=statesDict["CalOn"], dlo=30, dhi=30)
# plt.close()

data.qualityCheckLinefit("ZnKAlpha", positionToleranceFitSigma=7, worstAllowedFWHM=15, states=statesDict["CalOn"], dlo=30, dhi=30)
# data.qualityCheckLinefit("W3D", positionToleranceFitSigma=7, worstAllowedFWHM=15, states=statesDict["W_OFF"], dlo=30, dhi=30) #shouldnt filter on this because too low counts
data.plotHist(np.arange(800, 13000, 1.), "energy", states=statesDict["Cal"], coAddStates=False)

fig = plt.figure()
ax = fig.gca()
chan_range = data.keys()[0:5]
for Ds in chan_range:
    data[Ds].plotHist(np.arange(800, 13000, 1.), "energy", states=statesDict["CalOn"], coAddStates=True, axis=ax)
plt.legend(chan_range)


# iral1_line = mass.fluorescence_lines.SpectralLine.quick_monochromatic_line("Ir", 626.8, 0, 0)
# iral1_model = iral1_line.model(has_linear_background=False)
# iral1_prefix = 'ir_Al1_like_'
# iral1_model.prefix = iral1_prefix
# iral1_model.set_param_hint(name=f'{iral1_prefix}peak_ph', min=624.5, max=630) #set some reasonable limits
# iral1_model.set_param_hint(name=f'{iral1_prefix}fwhm', min=0.001, max=6)
# iral1_params = iral1_model.guess(ir_energies, bin_centers=ir_bin_centers, dph_de=1) #make initial guess based on data
# iral1_params[f'{iral1_prefix}dph_de'].set(value=1.0, vary=False)
# iral1_result = iral1_model.fit(ir_energies, params=iral1_params, bin_centers=ir_bin_centers) #we pass these params on to the composite model
# iral1_result.plotm()


# #zinc keys
# zn = [k for k, v in mass.STANDARD_FEATURES.items() if "Zn" in k]
# zn.extend([k for k, v in mass.spectra.items() if "Zn" in k])
# #germanium keys 
# ge = [k for k, v in mass.STANDARD_FEATURES.items() if "Ge" in k]
# ge.extend([k for k, v in mass.spectra.items() if "Ge" in k])

W = [k for k, v in mass.STANDARD_FEATURES.items() if "W" in k]
W.extend([k for k, v in mass.spectra.items() if "W" in k])

# def keys_energies(keys):
#     for key in keys:
#         try:
#             e = mass.STANDARD_FEATURES[key]
#         except:
#             e = mass.spectra[key]
#         print(f'{key}\t\t{e} \n')

if True:
    binsize=1
    WBinCenters, WData = data.hist(np.arange(800, 13000, binsize), "energy", states=statesDict["W_OFF"])
    import csv
    rows = zip(WBinCenters, WData)
    with open("W_20240723.txt", "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Bin center (eV)", "Counts per 2 eV bin"])
        for row in rows:
            writer.writerow(row)

    OsBinCenters, OsData = data.hist(np.arange(800, 13000, binsize), "energy", states=statesDict["Os_OFF"])
    rows = zip(OsBinCenters, OsData)
    with open("Os_20240723.txt", "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Bin center (eV)", "Counts per 2 eV bin"])
        for row in rows:
            writer.writerow(row)

    ReBinCenters, ReData = data.hist(np.arange(800, 13000, binsize), "energy", states=statesDict["Re_OFF"])
    rows = zip(ReBinCenters, ReData)
    with open("Re_20240723.txt", "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Bin center (eV)", "Counts per 2 eV bin"])
        for row in rows:
            writer.writerow(row)

data.plotHist(np.arange(800, 13000, 2.), "energy", states=statesDict["W_OFF"], coAddStates=True)
plt.title("20240723 W")