import os
import pylab as plt
import numpy as np
import pandas as pd
import mass
import sys
from mass.off import ChannelGroup, Channel, getOffFileListFromOneFile
from mass.calibration import _highly_charged_ion_lines
import mass.calibration.hci_models
#import mass.calibration._highly_charged_ion_lines
import lmfit

import CustomFluorescenceLines

#Laptop
#filename = r"C:\Users\Grant Mondeel\Box\CfA\TES\Ne-Like\20240724_off\0000\20240724_run0000_chan1.off"

#PC
filename = r"C:\Users\lamat\Box\CfA\TES\Ne-Like\20240724_off\0000\20240724_run0000_chan1.off"

plt.ion()
data = ChannelGroup(getOffFileListFromOneFile(filename, maxChans=999), verbose=True, channelClass = Channel, excludeStates=['IGNORE', 'STOP'])
data.setDefaultBinsize(1.)
mass.line_models.VALIDATE_BIN_SIZE = False ###with Galen's cal, had to go down to 0.025eV bins...
ds = data.firstGoodChannel()
#data.setOutputDir(baseDir=baseOutputDirectory, deleteAndRecreate=False)

###Some residual based cuts
#data.learnResidualStdDevCut() #see if this is works



#Make a dictionary with lists of aliases for each element.
#Pass this anywhere with the 'states=None' argument, e.g., data.plotHist(..., states=statesDict["W"])


statesDict = {
    "Re_OFF" : ["E_OFF", "H_OFF"],
    "Ir_ON" : ["I_ON"],
    "Ir_OFF" : ["I_OFF"],
    "Cal": ["G_OFF"],
    "CalOn": ["START_OFF", "G_OFF", "I_ON"] #states where the calibration source was turned on
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
ds.calibrationPlanAddPoint(8724, "AlKAlpha", states=calStates)
ds.calibrationPlanAddPoint(22563, "ScKAlpha", states=calStates)
# ds.calibrationPlanAddPoint(24414, "ScKBeta", states=calStates)
ds.calibrationPlanAddPoint(26826, "VKAlpha", states=calStates)
ds.calibrationPlanAddPoint(29016, "VKBeta", states=calStates)
ds.calibrationPlanAddPoint(31300, "MnKAlpha", states=calStates)
ds.calibrationPlanAddPoint(33600, "FeKAlpha", states=calStates)
ds.calibrationPlanAddPoint(36015, "CoKAlpha", states=calStates)
# ds.calibrationPlanAddPoint(36575, "FeKBeta", states=calStates)
# ds.calibrationPlanAddPoint(39140, "CoKBeta", states=calStates)
# ds.calibrationPlanAddPoint(40827, "CuKAlpha", states=calStates)
ds.calibrationPlanAddPoint(43223, "ZnKAlpha", states=calStates)
# # ds.calibrationPlanAddPoint(44335, "CuKBeta", states=calStates)
# ds.calibrationPlanAddPoint(47040, "ZnKBeta", states=calStates)
ds.calibrationPlanAddPoint(48125, "GeKAlphaCustom", states=calStates)
#ds.calibrationPlanAddPoint(52189, "GeKBeta", states=calStates)
ds.calibrationPlanAddPoint(52314, "GeKBeta", states="START_OFF")
### Check calibration on just one channel
# ds.calibrateFollowingPlan("filtValue", calibratedName="energy", dlo=50, dhi=50, approximate=True, overwriteRecipe=True)
# ds.diagnoseCalibration()
# ds.plotHist(np.arange(800, 12000, 1), "energy", states=["Cal"], coAddStates=False)

data.alignToReferenceChannel(referenceChannel=ds,
                            binEdges=np.arange(3000, 60000, 3), attr="filtValue")#, _rethrow=True)
data.cutAdd("cutForLearnDC", lambda energyRough: np.logical_and(
energyRough > 8000, energyRough < 11000), setDefault=False)

### Mass corrections.
# data.learnPhaseCorrection(indicatorName="filtPhase", uncorrectedName="filtValue", correctedName = "filtValuePC", states="Cal")#, cutRecipeName="cutForPC")
data.learnPhaseCorrection(indicatorName="filtPhase", uncorrectedName="filtValue", correctedName = "filtValuePC", states=statesDict["Cal"])
data.learnDriftCorrection(uncorrectedName="filtValuePC", indicatorName="pretriggerMean", correctedName="filtValueDC",
                            states=statesDict["CalOn"], cutRecipeName="cutForLearnDC")#, _rethrow=True)
data.calibrateFollowingPlan("filtValueDC", calibratedName="energy", dlo=70, dhi=40, approximate=True, overwriteRecipe=True)
data.qualityCheckLinefit("ZnKAlpha", positionToleranceFitSigma=3, worstAllowedFWHM=10, states=statesDict["Cal"], dlo=70, dhi=40)
ds.diagnoseCalibration()
data[6].markBad("bad")
# data.qualityCheckLinefit("GeKBeta", positionToleranceFitSigma=2, worstAllowedFWHM=12, states="Cal", dlo=50, dhi=40)

# calLines = ["AlKAlpha","ScKAlpha","VKAlpha","MnKAlpha",
#             "FeKAlpha","CoKAlpha","CuKAlpha"]#,"ZnKAlpha","GeKAlphaCustom"]
# for line in calLines:
#     data.qualityCheckLinefit(line, positionToleranceFitSigma=2, worstAllowedFWHM=14, states="Cal", dlo=50, dhi=50)
# plt.close()

data.plotHist(np.arange(800, 13000, 1.), "energy", states=statesDict["CalOn"], coAddStates=False)

fig = plt.figure()
ax = fig.gca()
chan_range = data.keys()[:]
for Ds in chan_range:
    data[Ds].plotHist(np.arange(8000, 13000, 5.), "energy", states=statesDict["CalOn"], coAddStates=False, axis=ax)
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

# def keys_energies(keys):
#     for key in keys:
#         try:
#             e = mass.STANDARD_FEATURES[key]
#         except:
#             e = mass.spectra[key]
#         print(f'{key}\t\t{e} \n')

# IrBinCenters, IrData = data.hist(np.arange(800, 13000, 1.), "energy", states="Ir 1")
# np.savetxt("IrSpect.txt", [IrBinCenters, IrData])

if False:
    import csv
    IrBinCenters, IrData = data.hist(np.arange(800, 13000, 2.), "energy", states=statesDict["Ir_OFF"])
    rows = zip(IrBinCenters, IrData)
    with open("Ir_20240724.txt", "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Bin center (eV)", "Counts per 2 eV bin"])
        for row in rows:
            writer.writerow(row)

# CalBinCenters, CalData = data.hist(np.arange(800, 13000, 1.), "energy", states=["Cal", "START"])
# np.savetxt("IrCalSpect.txt", [CalBinCenters, CalData])
fig = plt.figure()
ax = fig.gca()
data.plotHist(np.arange(800, 13000, 1.), "energy", states=statesDict["CalOn"], coAddStates=True, axis=ax)
data.plotHist(np.arange(800, 13000, 1.), "energy", states=statesDict["Ir_OFF"], coAddStates=False, axis=ax)
plt.title("20240724 Ir")