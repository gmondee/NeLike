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

#Laptop
filename = r"C:\Users\Grant Mondeel\Box\my EUV\tes\realtime\realtime\Spring2024\OffFiles\20181205_chan1.off"

#PC
#filename = r"C:\Users\lamat\Box\my EUV\tes\realtime\realtime\Spring2024\OffFiles\20181205_chan1.off"

plt.ion()
data = ChannelGroup(getOffFileListFromOneFile(filename, maxChans=999), verbose=True, channelClass = Channel, excludeStates=['START','A'])
data.setDefaultBinsize(0.3)
mass.line_models.VALIDATE_BIN_SIZE = False ###with Galen's cal, had to go down to 0.025eV bins...
ds = data.firstGoodChannel()
#data.setOutputDir(baseDir=baseOutputDirectory, deleteAndRecreate=False)

###Some residual based cuts
#data.learnResidualStdDevCut() #for some reason, this erases the majority of the pulses since the off files aren't good

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


###view filtValue plot
#data[93].plotHist(np.arange(0,13000,10), "filtValue", coAddStates=False)

###Check how many pulses are in the data file
# totalPulses = 0
# for chan in data.keys():
#     print("channel", chan,"has",len(data[chan].getAttr('filtValue','I')))
#     totalPulses += len(data[chan].getAttr('filtValue', 'I'))
# print("Combining all channels gives",totalPulses,"pulses.")


#Optional: assignes names to the states based on what was being observed
data.experimentStateFile.aliasState("B", "Ne")
data.experimentStateFile.aliasState("C", "W 1")
data.experimentStateFile.aliasState("D", "Os")
data.experimentStateFile.aliasState("E", "Ar")
data.experimentStateFile.aliasState("F", "Re")
data.experimentStateFile.aliasState("G", "W 2")
data.experimentStateFile.aliasState("H", "CO2")
data.experimentStateFile.aliasState("I", "Ir")

## Galen's full-spectrum energy calibration
if False:
    ds.calibrationPlanInit("filtValue")
    ds.calibrationPlanAddPoint(2128, "O He-Like 1s2s+1s2p", states="CO2")
    ds.calibrationPlanAddPoint(2421, "O H-Like 2p", states="CO2")
    ds.calibrationPlanAddPoint(2864, "O H-Like 3p", states="CO2")
    ds.calibrationPlanAddPoint(3404, "Ne He-Like 1s2s+1s2p", states="Ne")
    ds.calibrationPlanAddPoint(3768, "Ne H-Like 2p", states="Ne")
    ds.calibrationPlanAddPoint(5716, "W Ni-2", states=["W 1", "W 2"])
    ds.calibrationPlanAddPoint(6413, "W Ni-4", states=["W 1", "W 2"])
    ds.calibrationPlanAddPoint(7641, "W Ni-7", states=["W 1", "W 2"])
    ds.calibrationPlanAddPoint(10256, "W Ni-17", states=["W 1", "W 2"])
    # ds.calibrationPlanAddPoint(10700, "W Ni-20", states=["W 1", "W 2"]) #we exclude this line here to fit it after the calibration to evaluate our energy calibration
    ds.calibrationPlanAddPoint(11125, "Ar He-Like 1s2s+1s2p", states="Ar")
    ds.calibrationPlanAddPoint(11728, "Ar H-Like 2p", states="Ar")
    #ds.plotHist(np.arange(0, 4000, 1), "energyRough", coAddStates=False)

    #align the other channels to ds. You can find which channel this is by ds.channum
    data.alignToReferenceChannel(referenceChannel=ds,
                                binEdges=np.arange(500, 20000, 4), attr="filtValue")#, _rethrow=True)
    data.cutAdd("cutForLearnDC", lambda energyRough: np.logical_and(
    energyRough > 1000, energyRough < 3500), setDefault=False, _rethrow=True) #this focuses on the W states from 1000-3500eV. For Ir, we will focus on ~620eV on state I.
    data.learnDriftCorrection(uncorrectedName="filtValue", indicatorName="pretriggerMean", correctedName="filtValueDC",
                                states=["W 1", "W 2"], cutRecipeName="cutForLearnDC")#, _rethrow=True)
    data.calibrateFollowingPlan("filtValueDC", calibratedName="energy", dlo=10, dhi=10, approximate=False, overwriteRecipe=True)
    ds.diagnoseCalibration()
    data.linefit("W Ni-20", states=["W 1", "W 2"])#, dlo=5, dhi=5) #fit the W line in the W states and compare it to theory
    data.plotHist(np.arange(400, 4000, 0.5), "energy", coAddStates=False)

    ## energy2ph to build new calibration, ~500eV-800eV. Science line at 621eV.
    # mass.getmodel('O H-Like 2p').spect.nominal_peak_energy
    # energy_dict = dict()
    # energy_dict["O He-Like 1s2"]


def getEnergy(lineName):
    return mass.getmodel(lineName).spect.nominal_peak_energy
def getPH(lineName):
    energy=getEnergy(lineName)
    return ds.recipes['energy'].f.energy2ph(energy)
#Custom line ratios of O He-Like 1s2s+1s2p WXYZ region. Only using W Z separately and XY blended into one.
mass.calibration.fluorescence_lines.addline(
    element="O",
    material="Highly Charged Ion",
    linetype=" He-Like 1s2s+1s2p Ir",
    reference_short='NIST ASD',
    fitter_type=mass.calibration.fluorescence_lines.line_models.GenericLineModel,
    reference_plot_instrument_gaussian_fwhm=0.5,
    nominal_peak_energy=573.94777,
    energies=np.array([560.983, 568.551, 573.94777]), lorentzian_fwhm=np.array([0.1, 0.1, 0.1]),
    reference_amplitude=np.array([120, 77, 1000]),
    reference_amplitude_type=mass.calibration.fluorescence_lines.LORENTZIAN_PEAK_HEIGHT,
    ka12_energy_diff=None,
    reference_measurement_type="Experiment"
)
mass.calibration.fluorescence_lines.addline(
    element="O",
    material="Highly Charged Ion",
    linetype=" He-Like 1s2s+1s2p CO2",
    reference_short='NIST ASD',
    fitter_type=mass.calibration.fluorescence_lines.line_models.GenericLineModel,
    reference_plot_instrument_gaussian_fwhm=0.5,
    nominal_peak_energy=573.94777,
    energies=np.array([560.983, 568.551, 573.94777]), lorentzian_fwhm=np.array([0.1, 0.1, 0.1]),
    reference_amplitude=np.array([77, 120, 1000]),
    reference_amplitude_type=mass.calibration.fluorescence_lines.LORENTZIAN_PEAK_HEIGHT,
    ka12_energy_diff=None,
    reference_measurement_type="Experiment"
)
mass.calibration.fluorescence_lines.addline(
    element="O",
    material="Highly Charged Ion",
    linetype=" He-Like 1s2s+1s2p Total",
    reference_short='NIST ASD',
    fitter_type=mass.calibration.fluorescence_lines.line_models.GenericLineModel,
    reference_plot_instrument_gaussian_fwhm=0.5,
    nominal_peak_energy=573.94777,
    energies=np.array([560.983, 568.551, 573.94777]), lorentzian_fwhm=np.array([0.1, 0.1, 0.1]),
    reference_amplitude=np.array([88, 134, 1000]),
    reference_amplitude_type=mass.calibration.fluorescence_lines.LORENTZIAN_PEAK_HEIGHT,
    ka12_energy_diff=None,
    reference_measurement_type="Experiment"
)
#ds.recipes['energy'].f.energy2ph ## convert energy to pulse height from the calibration energy result named "energy"
## *.rowcount -> *.subframecount in massGui

#See lines with mass.spectra and mass.STANDARD_MODELS
## Build a narrow-range calibration from Dr. Takacs' line IDs
if True:
    # ###for ch55, looked good but only has 1 count in CO2 for 1s4p
    ## we can still try using all of these lines for the fit.
    ds=data[55]
    ds.calibrationPlanInit("filtValue")
    #ds.calibrationPlanAddPoint(2064, "O He-Like 1s2s+1s2p", states=["CO2","Ir"]) #this model has the wrong intensity ratios for this dataset
    
    #ds.calibrationPlanAddPoint(2016, "O7 1s.2s 3S J=1", states=["CO2", "Ir"]) #W line
    #ds.calibrationPlanAddPoint(2063, "O He-Like 1s2p 1P1", states=["CO2", "Ir"]) #Fitting the Z line (strongest, but blended)
    #ds.calibrationPlanAddPoint(2064, "O He-Like 1s2s+1s2p CO2", states="CO2") #Grant's custom model
    #ds.calibrationPlanAddPoint(2064, "O He-Like 1s2s+1s2p Ir", states="Ir") #Grant's custom model
    ds.calibrationPlanAddPoint(2064, "O He-Like 1s2s+1s2p Total", states=["CO2", "Ir"])
    ds.calibrationPlanAddPoint(2350, "O H-Like 2p", states=["CO2", "Ir"])
    #ds.calibrationPlanAddPoint(2385, "O7 1s.3p 3P* J=1", states="CO2") ##These lines are present but have low counts
    #ds.calibrationPlanAddPoint(2501, "O He-Like 1s4p 1P1", states="CO2") ##These lines are present but have low counts
    #ds.calibrationPlanAddPoint(2558, "O7 1s.5p 3P* J=2", states="CO2") ##These lines are present but have low counts
    ds.calibrationPlanAddPoint(2778, "O H-Like 3p", states=["CO2", "Ir"])
    #ds.calibrationPlanAddPoint(2926, "O H-Like 4p", states=["CO2","Ir"]) ##These lines are present but have low counts
    ds.calibrationPlanAddPoint(2990, "O8 5p 2P* J=1/2", states=["CO2","Ir"]) ##These lines are present but have low counts
    ds.calibrationPlanAddPoint(3293, "Ne He-Like 1s2s+1s2p", states="Ne")
    ds.calibrationPlanAddPoint(3642, "Ne H-Like 2p", states="Ne")

    data.alignToReferenceChannel(referenceChannel=ds,
                                binEdges=np.arange(500, 20000, 3), attr="filtValue")#, _rethrow=True)
    data.cutAdd("cutForLearnDC", lambda energyRough: np.logical_and(
    energyRough > 500, energyRough < 670), setDefault=False) #For Ir, we will focus on ~620eV on state I.

    ### Mass corrections. Phase correct appears to lower resolution and causes Ne 1s2s+1s2p to be off by ~2 eV, so I don't use it
    #data.learnPhaseCorrection(indicatorName="filtPhase", uncorrectedName="filtValue", correctedName = "filtValuePC")#, states=["CO2","Ir"])#, cutRecipeName="cutForPC")
    data.learnDriftCorrection(uncorrectedName="filtValue", indicatorName="pretriggerMean", correctedName="filtValueDC",
                                states=["CO2","Ir"], cutRecipeName="cutForLearnDC")#, _rethrow=True)
    data.calibrateFollowingPlan("filtValueDC", calibratedName="energy", dlo=4, dhi=4, approximate=True, overwriteRecipe=True)
    ds.diagnoseCalibration()
    #data.qualityCheckLinefit("O H-Like 2p", positionToleranceFitSigma=3, worstAllowedFWHM=5, states="Ir", dlo=5, dhi=5)
    ###If I've used the O He-like complex (~570 eV) then we need to cut channels that misidentify the line
    cutBinEdges=np.arange(583, 587.5, 0.3)
    for chan in data.keys():
        counts = sum(data[chan].hist(cutBinEdges, "energy", states=["CO2","Ir"])[1])
        if counts>5:
            data[chan].markBad(reason="failed cal of O He-like complex")
    data.plotHist(np.arange(400, 1000, 0.3), "energy", states=["CO2","Ir"], coAddStates=False)

    ##Finds channel with most counts (301)
    # maxCounts = 0
    # maxDs = 0
    # for Ds in data:
    #     if len(data[Ds].energy)>maxCounts:
    #         maxCounts=len(data[Ds].energy)
    #         maxDs=Ds



ir_bin_edges = np.arange(550, 715, 0.3) #original: 616, 626, 0.3
ir_bin_centers, ir_energies = data.hist(ir_bin_edges, "energy", states=["Ir"])
#ir_energies = plt.hist(ir_hist, )

### Calibration 1: only fit the D2 line with a gaussian at ~622eV
irD2_line = mass.fluorescence_lines.SpectralLine.quick_monochromatic_line("Ir", 622, 0, 0)
irD2_model = irD2_line.model()
irD2_prefix = 'ir_D2_'
irD2_model.prefix = irD2_prefix
irD2_model.nominal_peak_energy = 621.16952 #eV
irD2_model.set_param_hint(name=f'{irD2_prefix}peak_ph', min=617, max=624) #set some reasonable limits
irD2_model.set_param_hint(name=f'{irD2_prefix}fwhm', min=0.001, max=20)
irD2_model.set_param_hint(name=f'{irD2_prefix}background', min=0.000, max=20)

ir_params = irD2_model.guess(ir_energies, bin_centers=ir_bin_centers, dph_de=1) #make initial guess based on data


ir_params[f'{irD2_prefix}dph_de'].set(value=1.0, vary=False)

#ir_result = ir_model.fit(ir_energies, params=ir_params, bin_centers=ir_bin_centers) #perform single-line fits based on data
#ir_result.plotm()

### Calibration 2: fit the D2 line and the blended Ir line just above the D2 line
# iral1_bin_edges = np.arange(616, 640, 0.3)
# iral1_bin_centers, iral1_energies = data.hist(ir_bin_edges, "energy", states=["Ir"])
iral1_line = mass.fluorescence_lines.SpectralLine.quick_monochromatic_line("Ir", 626.8, 0, 0)
iral1_model = iral1_line.model(has_linear_background=False)
iral1_prefix = 'ir_Al1_like_'
iral1_model.prefix = iral1_prefix
iral1_model.set_param_hint(name=f'{iral1_prefix}peak_ph', min=624.5, max=630) #set some reasonable limits
iral1_model.set_param_hint(name=f'{iral1_prefix}fwhm', min=0.001, max=6)
iral1_params = iral1_model.guess(ir_energies, bin_centers=ir_bin_centers, dph_de=1) #make initial guess based on data
iral1_params[f'{iral1_prefix}dph_de'].set(value=1.0, vary=False)
#iral1_result = iral1_model.fit(ir_energies, params=iral1_params, bin_centers=ir_bin_centers) #we pass these params on to the composite model
#iral1_result.plotm()

# irCompositeModel = ir_model + iral1_model
# irCompositeParams = ir_params + iral1_params
# nominal_separation = abs(ir_model.nominal_peak_energy - 626)
# irCompositeParams[f'{iral1_prefix}fwhm'].expr = f'{ir_prefix}fwhm'  #we can set a fwhm ratio. here, al-like ir has many fewer counts.
# #irCompositeParams[f'{iral1_prefix}peak_ph'].expr = f'{ir_prefix}peak_ph + {nominal_separation}' #if we know one of the lines well, we can set the separation

# irCompositeResult = irCompositeModel.fit(ir_energies, params=irCompositeParams, bin_centers=ir_bin_centers)
# irCompositeResult.plotm()
# for model_component in irCompositeResult.eval_components().values():
#     plt.plot(ir_bin_centers, model_component)
# print(irCompositeResult.best_values)

### 3: Add Mg-like Ir and other Al-like Ir line
iral2_line = mass.fluorescence_lines.SpectralLine.quick_monochromatic_line("Ir", 638, 0, 0)
iral2_model = iral2_line.model(has_linear_background=False)
iral2_prefix = 'ir_Al2_like_'
iral2_model.prefix = iral2_prefix
iral2_model.set_param_hint(name=f'{iral2_prefix}peak_ph', min=636, max=640) #set some reasonable limits
iral2_model.set_param_hint(name=f'{iral2_prefix}fwhm', min=0.001, max=6)
iral2_params = iral2_model.guess(ir_energies, bin_centers=ir_bin_centers, dph_de=1) #make initial guess based on data
iral2_params[f'{iral2_prefix}dph_de'].set(value=1.0, vary=False)
#iral2_result = iral2_model.fit(ir_energies, params=iral2_params, bin_centers=ir_bin_centers) #we pass these params on to the composite model

irMg_line = mass.fluorescence_lines.SpectralLine.quick_monochromatic_line("Ir", 633.5, 0, 0)
irMg_model = irMg_line.model(has_linear_background=False)
irMg_prefix = 'ir_Mg_like_'
irMg_model.prefix = irMg_prefix
irMg_model.set_param_hint(name=f'{irMg_prefix}peak_ph', min=630.5, max=636) #set some reasonable limits
irMg_model.set_param_hint(name=f'{irMg_prefix}fwhm', min=0.001, max=6)
irMg_params = irMg_model.guess(ir_energies, bin_centers=ir_bin_centers, dph_de=1) #make initial guess based on data
irMg_params[f'{irMg_prefix}dph_de'].set(value=1.0, vary=False)
#irMg_result = irMg_model.fit(ir_energies, params=irMg_params, bin_centers=ir_bin_centers) #we pass these params on to the composite model

irNa_line = mass.fluorescence_lines.SpectralLine.quick_monochromatic_line("Ir", 647, 0, 0) ###Blended with O He-like 1s4p
irNa_model = irNa_line.model(has_linear_background=False)
irNa_prefix = "irNa_"
irNa_model.prefix = irNa_prefix
irNa_model.set_param_hint(name=f'{irNa_prefix}peak_ph', min=640, max=652) #set some reasonable limits
#irNa_model.set_param_hint(name=f'{irNa_prefix}fwhm', min=0.001, max=6)
irNa_params = irNa_model.guess(ir_energies, bin_centers=ir_bin_centers, dph_de=1) #make initial guess based on data
irNa_params[f'{irNa_prefix}dph_de'].set(value=1.0, vary=False)

#O2p = mass.fluorescence_lines.SpectralLine. #get the mass model
O2p_model = mass.spectra["O H-Like 2p"].model(has_linear_background=False)
O2p_prefix = "O2p_"
O2p_model.prefix = O2p_prefix
O2p_model.set_param_hint(name=f'{O2p_prefix}peak_ph', min=650, max=657) #set some reasonable limits
#O2p_model.set_param_hint(name=f'{O2p_prefix}fwhm', min=0.001, max=6)
O2p_params = O2p_model.guess(ir_energies, bin_centers=ir_bin_centers, dph_de=1) #make initial guess based on data
O2p_params[f'{O2p_prefix}dph_de'].set(value=1.0, vary=False)
#O2p_result = O2p_model.fit(ir_energies, params=O2p_params, bin_centers=ir_bin_centers) #we pass these params on to the composite model

O1s3p_model = mass.spectra["O He-Like 1s3p 1P1"].model(has_linear_background=False)
O1s3p_prefix = "O3p_"
O1s3p_model.prefix = O1s3p_prefix
O1s3p_model.set_param_hint(name=f'{O1s3p_prefix}peak_ph', min=660, max=670) #set some reasonable limits
O1s3p_model.set_param_hint(name=f'{O1s3p_prefix}fwhm', min=2, max=4.5)
#O3p_model.set_param_hint(name=f'{O3p_prefix}fwhm', min=0.001, max=6)
O1s3p_params = O1s3p_model.guess(ir_energies, bin_centers=ir_bin_centers, dph_de=1) #make initial guess based on data
O1s3p_params[f'{O1s3p_prefix}dph_de'].set(value=1.0, vary=False)
#O3p_result = O3p_model.fit(ir_energies, params=O3p_params, bin_centers=ir_bin_centers) #we pass these params on to the composite model

iral3_line = mass.fluorescence_lines.SpectralLine.quick_monochromatic_line("Ir", 686.30, 0, 0) #3p1 to 3d1 transition of Al-like Ir
iral3_model = iral3_line.model(has_linear_background=False)
iral3_prefix = 'ir_Al3_like_'
iral3_model.prefix = iral3_prefix
iral3_model.set_param_hint(name=f'{iral3_prefix}peak_ph', min=680, max=691) #set some reasonable limits
iral3_model.set_param_hint(name=f'{iral3_prefix}fwhm', min=0.001, max=6)
iral3_params = iral3_model.guess(ir_energies, bin_centers=ir_bin_centers, dph_de=1) #make initial guess based on data
iral3_params[f'{iral3_prefix}dph_de'].set(value=1.0, vary=False)
#iral3_result = iral3_model.fit(ir_energies, params=iral3_params, bin_centers=ir_bin_centers) #we pass these params on to the composite model

O1s4p_model = mass.spectra["O He-Like 1s4p 1P1"].model(has_linear_background=False)
O1s4p_prefix = "O1s4p_"
O1s4p_model.prefix = O1s4p_prefix
O1s4p_model.set_param_hint(name=f'{O1s4p_prefix}peak_ph', min=695, max=699) #set some reasonable limits
O1s4p_model.set_param_hint(name=f'{O1s4p_prefix}fwhm', min=2, max=4.5)
#O1s4p_model.set_param_hint(name=f'{O1s4p_prefix}fwhm', min=0.001, max=6)
O1s4p_params = O1s4p_model.guess(ir_energies, bin_centers=ir_bin_centers, dph_de=1) #make initial guess based on data
O1s4p_params[f'{O1s4p_prefix}dph_de'].set(value=1.0, vary=False)
#O1s4p_result = O1s4p_model.fit(ir_energies, params=O1s4p_params, bin_centers=ir_bin_centers) #we pass these params on to the composite model

irSi_line = mass.fluorescence_lines.SpectralLine.quick_monochromatic_line("Ir", 701, 0, 0) ###Blended with O He-like 1s4p
irSi_model = irSi_line.model(has_linear_background=False)
irSi_prefix = "irSi_"
irSi_model.prefix = irSi_prefix
irSi_model.set_param_hint(name=f'{irSi_prefix}peak_ph', min=699, max=703) #set some reasonable limits
#irSi_model.set_param_hint(name=f'{irSi_prefix}fwhm', min=0.001, max=6)
irSi_params = irSi_model.guess(ir_energies, bin_centers=ir_bin_centers, dph_de=1) #make initial guess based on data
irSi_params[f'{irSi_prefix}dph_de'].set(value=1.0, vary=False)
#irSi_result = irSi_model.fit(ir_energies, params=irSi_params, bin_centers=ir_bin_centers) #we pass these params on to the composite model

O1s5p_model = mass.spectra["O7 1s.5p 3P* J=1"].model(has_linear_background=False)
O1s5p_prefix = "O1s5p_"
O1s5p_model.prefix = O1s5p_prefix
O1s5p_model.set_param_hint(name=f'{O1s5p_prefix}peak_ph', min=710, max=714) #set some reasonable limits
O1s5p_model.set_param_hint(name=f'{O1s5p_prefix}fwhm', min=2, max=5)
#O1s5p_model.set_param_hint(name=f'{O1s5p_prefix}fwhm', min=0.001, max=6)
O1s5p_params = O1s5p_model.guess(ir_energies, bin_centers=ir_bin_centers, dph_de=1) #make initial guess based on data
O1s5p_params[f'{O1s5p_prefix}dph_de'].set(value=1.0, vary=False)
#O1s5p_result = O1s5p_model.fit(ir_energies, params=O1s5p_params, bin_centers=ir_bin_centers) #we pass these params on to the composite model

O_he_model = mass.spectra["O He-Like 1s2s+1s2p Ir"].model(has_linear_background=False)
O_he_prefix = "O_he_"
O_he_model.prefix = O_he_prefix
O_he_model.set_param_hint(name=f'{O_he_prefix}peak_ph', min=550, max=600) #set some reasonable limits
#O_he_model.set_param_hint(name=f'{O_he_prefix}fwhm', min=0.001, max=6)
O_he_params = O_he_model.guess(ir_energies, bin_centers=ir_bin_centers, dph_de=1) #make initial guess based on data
O_he_params[f'{O_he_prefix}dph_de'].set(value=1.0, vary=False)
O_he_result = O_he_model.fit(ir_energies, params=O_he_params, bin_centers=ir_bin_centers) #we pass these params on to the composite model

# Composite model
irCompositeModel = irD2_model + iral1_model +iral2_model + irMg_model + irNa_model + O2p_model +O1s3p_model +  O1s4p_model + irSi_model + iral3_model + O1s5p_model + O_he_model
irCompositeParams = ir_params + iral1_params +iral2_params + irMg_params + irNa_params + O2p_params + O1s3p_params + O1s4p_params + irSi_params +iral3_params + O1s5p_params+ O_he_params
#nominal_separation = abs(irD2_model.nominal_peak_energy - 626)
irCompositeParams[f'{iral1_prefix}fwhm'].expr = f'{irD2_prefix}fwhm'  #we can set a fwhm
irCompositeParams[f'{iral2_prefix}fwhm'].expr = f'{irD2_prefix}fwhm'
irCompositeParams[f'{irMg_prefix}fwhm'].expr = f'{irD2_prefix}fwhm'
irCompositeParams[f'{irNa_prefix}fwhm'].expr = f'{irD2_prefix}fwhm'
irCompositeParams[f'{iral3_prefix}fwhm'].expr = f'{irD2_prefix}fwhm'
irCompositeParams[f'{O_he_prefix}fwhm'].expr = f'{irD2_prefix}fwhm'
irCompositeParams[f'{O1s4p_prefix}fwhm'].expr = f'{irD2_prefix}fwhm'
irCompositeParams[f'{O1s5p_prefix}fwhm'].expr = f'{irD2_prefix}fwhm'
irCompositeParams[f'{irSi_prefix}fwhm'].expr = f'{irD2_prefix}fwhm'
#irCompositeParams[f'{O2p_prefix}fwhm'].expr = f'{ir_prefix}fwhm' #the O H-like 2p model is a voigt, so there is a lorentzian fwhm and sigma. harder to set?
#irCompositeParams[f'{O1s4p_prefix}peak_ph'].expr = f'{O2p_prefix}peak_ph + {O1s4p_model.nominal_peak_energy - O2p_model.nominal_peak_energy}' #if we know one of the lines well, we can set the separation

irCompositeResult = irCompositeModel.fit(ir_energies, params=irCompositeParams, bin_centers=ir_bin_centers)
irCompositeResult.plotm()
legend=["Data", "Best fit"]
for model_prefix in irCompositeResult.eval_components():
    plt.plot(ir_bin_centers, irCompositeResult.eval_components()[model_prefix])
    legend.append(model_prefix[:-1])
plt.legend(legend)
print(irCompositeResult.best_values)
ir_background = irCompositeResult.best_values['ir_D2_background']


if True: #Use a model to get the He-like oxygen line energies
    OHe_bin_edges = np.arange(556, 580, 0.3) 
    OHe_bin_centers, OHe_energies = data.hist(OHe_bin_edges, "energy", states=["Ir"])

    O_He_Z_line = mass.fluorescence_lines.SpectralLine.quick_monochromatic_line("Ir", 561, 0, 0) 
    O_He_Z_model = O_He_Z_line.model(has_linear_background=True)
    O_He_Z_prefix = "O_He_Z_"
    O_He_Z_model.prefix = O_He_Z_prefix
    O_He_Z_model.set_param_hint(name=f'{O_He_Z_prefix}peak_ph', min=557, max=564) #set some reasonable limits
    O_He_Z_model.set_param_hint(name=f'{O_He_Z_prefix}fwhm', min=0.001, max=4.5)
    O_He_Z_params = O_He_Z_model.guess(ir_energies, bin_centers=ir_bin_centers, dph_de=1) #make initial guess based on data
    O_He_Z_params[f'{O_He_Z_prefix}dph_de'].set(value=1.0, vary=False)
    #O_He_Z_result = O_He_Z_model.fit(ir_energies, params=O_He_Z_params, bin_centers=ir_bin_centers) #we pass these params on to the composite model

    O_He_W_line = mass.fluorescence_lines.SpectralLine.quick_monochromatic_line("Ir", 573.95, 0, 0) 
    O_He_W_model = O_He_W_line.model(has_linear_background=False)
    O_He_W_prefix = "O_He_W_"
    O_He_W_model.prefix = O_He_W_prefix
    O_He_W_model.set_param_hint(name=f'{O_He_W_prefix}peak_ph', min=570, max=578) #set some reasonable limits
    O_He_W_model.set_param_hint(name=f'{O_He_W_prefix}fwhm', min=0.001, max=4.5)
    O_He_W_params = O_He_W_model.guess(ir_energies, bin_centers=ir_bin_centers, dph_de=1) #make initial guess based on data
    O_He_W_params[f'{O_He_W_prefix}dph_de'].set(value=1.0, vary=False)
    #O_He_W_result = O_He_W_model.fit(ir_energies, params=O_He_W_params, bin_centers=ir_bin_centers) #we pass these params on to the composite model

    O_He_YX_line = mass.fluorescence_lines.SpectralLine.quick_monochromatic_line("Ir", 568, 0, 0) 
    O_He_YX_model = O_He_YX_line.model(has_linear_background=False)
    O_He_YX_prefix = "O_He_YX_"
    O_He_YX_model.prefix = O_He_YX_prefix
    O_He_YX_model.set_param_hint(name=f'{O_He_YX_prefix}peak_ph', min=564, max=570) #set some reasonable limits
    O_He_YX_model.set_param_hint(name=f'{O_He_YX_prefix}fwhm', min=0.001, max=8)
    O_He_YX_params = O_He_YX_model.guess(ir_energies, bin_centers=ir_bin_centers, dph_de=1) #make initial guess based on data
    O_He_YX_params[f'{O_He_YX_prefix}dph_de'].set(value=1.0, vary=False)
    #O_He_YX_result = O_He_YX_model.fit(ir_energies, params=O_He_YX_params, bin_centers=ir_bin_centers) #we pass these params on to the composite model

    OHeComposite = O_He_Z_model + O_He_W_model + O_He_YX_model
    OHeParams = O_He_Z_params + O_He_W_params + O_He_YX_params
    OHeParams[f'{O_He_Z_prefix}background'].set(value=ir_background, vary=False) 
    OHeParams[f'{O_He_Z_prefix}peak_ph'].expr = f'{O_He_W_prefix}peak_ph - {573.94777 - 561.072349}'
    #OHeCompositeResult = OHeComposite.fit(OHe_energies, params=OHeParams, bin_centers=cal.f3(OHe_bin_centers, *cal.popt)) #poly
    #OHeCompositeResult = OHeComposite.fit(OHe_energies, params=OHeParams, bin_centers=calLine.f1(OHe_bin_centers, *calLine.popt))
    OHeCompositeResult = OHeComposite.fit(OHe_energies, params=OHeParams, bin_centers=OHe_bin_centers)
    OHeCompositeResult.plotm()
    legend=["Data", "Best fit"]
    for model_prefix in OHeCompositeResult.eval_components():
        #plt.plot(calLine.f1(OHe_bin_centers, *calLine.popt), OHeCompositeResult.eval_components()[model_prefix])
        plt.plot(OHe_bin_centers, OHeCompositeResult.eval_components()[model_prefix])
        legend.append(model_prefix[:-1])
    plt.legend(legend)
    print(OHeCompositeResult.best_values)

    Zenergy = OHeCompositeResult.params['O_He_Z_peak_ph'].value
    Zunc = OHeCompositeResult.params['O_He_Z_peak_ph'].stderr
    Wenergy = OHeCompositeResult.params['O_He_W_peak_ph'].value
    Wunc = OHeCompositeResult.params['O_He_W_peak_ph'].stderr

#  For Hunter's polynomial calibraition [Mass energy, mass unc,   lit energy, lit unc]. Needs >= 5 lines.
# calArr = [[561.73,      0.18,       561.072349,	0.03], #outlier, makes residuals worse if it is forced to be included
#           [573.837,	    0.04,	    573.96106,	0.03], 
#           [653.658297,	0.02471115,	653.6801,	0.0000017],
#           [665.594,	    0.039,	    665.5744,	0.04],
#           [697.88,  	0.075,	    697.7859,	0.0013],
#           [712.7,   	0.13,	    712.7221,	0.0017],
#           [774.615,	    0.049,	    774.634218,	0.0000013],
#           [816.823,	    0.061,	    816.9745429,0.0000013],
#           [836.6,   	0.088,  	836.5722377,0.0000013]]   
# for i, row in enumerate(calArr): #swap the order to use w/ hunter's code 
#     newLine = [row[2], row[3], row[0], row[1]]
#     calArr[i] = newLine

# ##Quickly plot all of the lines for a polynomial fit
# linesForPoly = ['O7 1s.2s 3S J=1', 'O7 1s.2p 1P* J=1', "O H-Like 2p", "O He-Like 1s3p 1P1", "O H-Like 3p",'O H-Like 4p','O8 5p 2P* J=1/2']
linesForPoly = ["O H-Like 2p", "O He-Like 1s3p 1P1", "O H-Like 3p",'O H-Like 4p','O8 5p 2P* J=1/2']
#data.linefit("O He-Like 1s2s+1s2p Ir", states="Ir", dlo=4, dhi=4)
lineEnergies = np.empty(0)
lineErrs = np.empty(0)
for line in linesForPoly:
    lineFit = data.linefit(line, states="Ir", dlo=4, dhi=4, plot=False)
    lineEnergy = lineFit.params['peak_ph'].value
    lineErr = lineFit.params['peak_ph'].stderr
    lineEnergies=np.append(lineEnergies,lineEnergy)
    lineErrs=np.append(lineErrs,lineErr)
    # print(line, lineEnergy)
    # print(line, mass.spectra[line].peak_energy)#, mass.getmodel(line).spect.nominal_peak_energy)
#cal_csv = pd.read_csv('Mass2Poly.csv')
cal_csv = pd.read_csv('Mass2PolyPeakEnergies.csv')

#massE = np.array(cal_csv['Mass energy'])[list(cal_csv['For Calibration'])]
massE = np.append(lineEnergies, 0) #add the 0,0 point
massE = np.append([Zenergy, Wenergy],massE)

#massUnc = np.array(cal_csv['Mass uncertainty'])[list(cal_csv['For Calibration'])]
massUnc = np.append(lineErrs, 0.01)
massUnc = np.append([Zunc, Wunc], massUnc)

refE = np.array(cal_csv['Ref energy'])[list(cal_csv['For Calibration'])]
refUnc = np.array(cal_csv['Ref uncertainty'])[list(cal_csv['For Calibration'])]

calArr=[]
for mE, mU, rE, rU in zip(massE, massUnc, refE, refUnc):
    calArr.append([rE,rU,mE,mU])

import calibrators as cals
###Polynomial calibration
# cal = cals.Calibrator(calArr)
# cal.fit(verbose=True, outlier_sigma=None)
# cal.residual_plot()
###Linear calibration
# calLine = cals.Calibrator(calArr)
# calLine.fit(verbose=True, outlier_sigma=None, funcs='f1')
# calLine.residual_plot(funcs='f1')

if False: #3rd order polynomial
    cal = cals.Calibrator(calArr)
    cal.fit(verbose=True, outlier_sigma=None)
    cal.residual_plot()
    irPoly = cal.f3(ir_bin_centers, *cal.popt)
    irCompositeResult = irCompositeModel.fit(ir_energies, params=irCompositeParams, bin_centers=irPoly)
    irCompositeResult.plotm()
    legend=["Data", "Best fit"]
    for model_prefix in irCompositeResult.eval_components():
        plt.plot(irPoly, irCompositeResult.eval_components()[model_prefix])
        legend.append(model_prefix[:-1])
    plt.legend(legend)
else: #linear
    calLine = cals.Calibrator(calArr)
    calLine.fit(verbose=True, outlier_sigma=2, funcs='f1')
    calLine.residual_plot(funcs='f1')
    irLine = calLine.f1(ir_bin_centers, *calLine.popt)
    irCompositeResult = irCompositeModel.fit(ir_energies, params=irCompositeParams, bin_centers=irLine)
    irCompositeResult.plot() #irCompositeResult.plotm() for results
    legend=["Data", "Best fit"]
    for model_prefix in irCompositeResult.eval_components():
        plt.plot(irLine, irCompositeResult.eval_components()[model_prefix])
        #legend.append(model_prefix[:-1])
    plt.legend(legend)

# compositeModel = model1 + model2 #sum the two models together
# compositeParams = result1.params + result2.params #store the parameters from each individual model's fit
# #compositeParams[f'{prefix1}fwhm'].expr = f'{prefix2}fwhm'
# compositeParams[f'{prefix1}peak_ph'].expr = f'{prefix2}peak_ph - {nominal_separation}' #make relations between lines (expect a known separation)
# compositeParams.add(name='ampRatio', value=1, vary=True)
# compositeParams[f'{prefix1}integral'].expr = f'{prefix2}integral * ampRatio'
# compositeResult = compositeModel.fit(
#     sim_energies, params=compositeParams, bin_centers=bin_centers) #fit the composite model
# compositeResult._validate_bins_per_fwhm(minimum_bins_per_fwhm=3)
# compositeResult.plot()
# print(compositeResult.best_values)


#Narrow fit of the D2 and first Mg-like line for best position uncertainty. Uses linear calibration.
if False:
    D2_bin_edges = np.arange(616, 628, 0.3)
    D2_bin_centers, D2_energies = data.hist(D2_bin_edges, "energy", states=["Ir"])

    D2MgModel = irD2_model + iral1_model
    D2CompositeParams = ir_params + iral1_params
    D2CompositeParams[f'{iral1_prefix}fwhm'].expr = f'{irD2_prefix}fwhm'  #we can set a fwhm
    #D2CompositeResult = D2MgModel.fit(D2_energies, params=D2CompositeParams, bin_centers=cal.f3(D2_bin_centers, *cal.popt)) #poly cal
    D2CompositeResult = D2MgModel.fit(D2_energies, params=D2CompositeParams, bin_centers=calLine.f1(D2_bin_centers, *calLine.popt)) #linear
    D2CompositeResult.plotm()
    legend=["Data", "Best fit"]
    for model_prefix in D2CompositeResult.eval_components():
        #plt.plot(cal.f3(D2_bin_centers, *cal.popt), D2CompositeResult.eval_components()[model_prefix]) #poly
        plt.plot(calLine.f1(D2_bin_centers, *calLine.popt), D2CompositeResult.eval_components()[model_prefix]) #linear
        legend.append(model_prefix[:-1])
    plt.legend(legend)
    print(D2CompositeResult.best_values)
else: #D2 line, and first Mg- and Al- like lines. Uses linear calibration
    # D2_bin_edges = np.arange(616, 642, 0.3) 
    D2_bin_edges = np.arange(616, 642, 0.3)  
    D2_bin_centers, D2_energies = data.hist(D2_bin_edges, "energy", states=["Ir"])

    D2MgModel = irD2_model + iral1_model +irMg_model + iral2_model
    D2CompositeParams = ir_params + iral1_params + irMg_params + iral2_params
    D2CompositeParams[f'{iral1_prefix}fwhm'].expr = f'{irD2_prefix}fwhm'  #we can set a fwhm
    D2CompositeParams[f'{iral2_prefix}fwhm'].expr = f'{irD2_prefix}fwhm'
    D2CompositeParams[f'{irMg_prefix}fwhm'].expr = f'{irD2_prefix}fwhm'
    # fwhm=4.4
    # D2CompositeParams[f'{irD2_prefix}fwhm'].set(min=4.3, max=fwhm, vary=True)
    # D2CompositeParams[f'{iral1_prefix}fwhm'].set(min=4.3, max=fwhm, vary=True)
    # D2CompositeParams[f'{iral2_prefix}fwhm'].set(min=4.3, max=fwhm, vary=True)
    # D2CompositeParams[f'{irMg_prefix}fwhm'].set(min=4.3, max=fwhm, vary=True)
    # D2CompositeParams[f'{irD2_prefix}fwhm'].set(value=fwhm, vary=True)
    # D2CompositeParams[f'{iral1_prefix}fwhm'].set(value=fwhm, vary=True)
    # D2CompositeParams[f'{iral2_prefix}fwhm'].set(value=fwhm, vary=True)
    # D2CompositeParams[f'{irMg_prefix}fwhm'].set(value=fwhm, vary=True)
    
    D2CompositeParams[f'{irD2_prefix}background'].set(value=ir_background, vary=False) #Use the background seen between 500-600 eV instead of letting it be a fitted parameter

    #D2CompositeResult = D2MgModel.fit(D2_energies, params=D2CompositeParams, bin_centers=cal.f3(D2_bin_centers, *cal.popt)) #poly cal
    #D2CompositeResult = D2MgModel.fit(D2_energies, params=D2CompositeParams, bin_centers=calLine.f1(D2_bin_centers, *calLine.popt)) #linear
    D2CompositeResult = D2MgModel.fit(D2_energies, params=D2CompositeParams, bin_centers=D2_bin_centers) #none
    
    D2CompositeResult.plotm()
    legend=["Data", "Best fit"]
    for model_prefix in D2CompositeResult.eval_components():
        #plt.plot(cal.f3(D2_bin_centers, *cal.popt), D2CompositeResult.eval_components()[model_prefix]) #poly
        #plt.plot(calLine.f1(D2_bin_centers, *calLine.popt), D2CompositeResult.eval_components()[model_prefix]) #linear
        plt.plot(D2_bin_centers, D2CompositeResult.eval_components()[model_prefix])#none
        legend.append(model_prefix[:-1])
    plt.legend(legend)
    print(D2CompositeResult.best_values)
    D2energy = D2CompositeResult.params['irD2_peak_ph'].value
    D2unc = D2CompositeResult.params['irD2_peak_ph'].stderr
    print(f'Bin size: 0.3 eV \n   {D2energy=} \n   {D2unc=}')

    if False: ###Option to include H-like O 2p->1s. I found that the O2p tail is vanishingly small even for the closest Ir line
        D2_bin_edges = np.arange(616, 660, 0.3)
        D2_bin_centers, D2_energies = data.hist(D2_bin_edges, "energy", states=["Ir"])

        D2MgModel = irD2_model + iral1_model +irMg_model + iral2_model + O2p_model
        D2CompositeParams = ir_params + iral1_params + irMg_params + iral2_params + O2p_params
        D2CompositeParams[f'{iral1_prefix}fwhm'].expr = f'{irD2_prefix}fwhm'  #we can set a fwhm
        D2CompositeParams[f'{iral2_prefix}fwhm'].expr = f'{irD2_prefix}fwhm'
        D2CompositeParams[f'{irMg_prefix}fwhm'].expr = f'{irD2_prefix}fwhm'
        # D2CompositeParams[f'{irD2_prefix}fwhm'].set(max=4)
        # D2CompositeParams[f'{iral1_prefix}fwhm'].set(max=4)
        # D2CompositeParams[f'{iral2_prefix}fwhm'].set(max=4)
        # D2CompositeParams[f'{irMg_prefix}fwhm'].set(max=4)
        D2CompositeParams[f'{irD2_prefix}background'].set(value=ir_background, vary=False)

        #D2CompositeResult = D2MgModel.fit(D2_energies, params=D2CompositeParams, bin_centers=cal.f3(D2_bin_centers, *cal.popt)) #poly cal
        #D2CompositeResult = D2MgModel.fit(D2_energies, params=D2CompositeParams, bin_centers=calLine.f1(D2_bin_centers, *calLine.popt)) #linear
        D2CompositeResult = D2MgModel.fit(D2_energies, params=D2CompositeParams, bin_centers=D2_bin_centers) #none
        
        D2CompositeResult.plotm()
        legend=["Data", "Best fit"]
        for model_prefix in D2CompositeResult.eval_components():
            #plt.plot(cal.f3(D2_bin_centers, *cal.popt), D2CompositeResult.eval_components()[model_prefix]) #poly
            #plt.plot(calLine.f1(D2_bin_centers, *calLine.popt), D2CompositeResult.eval_components()[model_prefix]) #linear
            plt.plot(D2_bin_centers, D2CompositeResult.eval_components()[model_prefix]) 
            legend.append(model_prefix[:-1])
        plt.legend(legend)
        print(D2CompositeResult.best_values)



# D2_bin_edges = np.arange(616, 624, 0.3) #original: 616, 626, 0.3
# D2_bin_centers, D2_energies = data.hist(D2_bin_edges, "energy", states=["Ir"])

# irD2_model
# D2Result = irD2_model.fit(D2_energies, params=ir_params, bin_centers=cal.f3(D2_bin_centers, *cal.popt))
# D2Result.plotm()
# legend=["Data", "Best fit"]
# for model_prefix in D2Result.eval_components():
#     plt.plot(D2_bin_centers, D2Result.eval_components()[model_prefix])
#     legend.append(model_prefix[:-1])
# plt.legend(legend)
# print(D2Result.best_values)





#readp = pd.read_pickle(r'C:\Users\lamat\Box\my EUV\tes\realtime\realtime\src\mass\mass\calibration\nist_asd.pickle')



# spect_bin_edges = np.arange(400, 2000, 0.3) 
# spect_bin_centers, spect_energies = data.hist(spect_bin_edges, "energy", states=["Ir"])
