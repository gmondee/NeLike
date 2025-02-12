import mass
import h5py
import numpy as np
import time

def data_load5LagInfo(self, model_hdf5_path):
    with h5py.File(model_hdf5_path,"r") as h5:
        models = {int(ch) : mass.pulse_model.PulseModel.fromHDF5(h5[ch]) for ch in h5.keys()}
    for channum, ds in self.items():
        # define recipes for "filtValue5Lag", "peakX5Lag" and "cba5Lag"
        # where cba refers to the coefficiencts of a polynomial fit to the 5 lags of the filter
        ds.model = models[ds.channum]
        filter_5lag = ds.model.f_5lag
        ds.add5LagRecipes(filter_5lag)

def ds_predictFV(self, energy, attr="energy"):
    return self.recipes[attr].f.energy2ph(energy)
mass.off.ChannelGroup.load5LagInfo = data_load5LagInfo
mass.off.Channel.predictFV = ds_predictFV

def load_timing(timing_path):
    cal_times = np.genfromtxt(timing_path,skip_header=1,usecols=0)
    onoff = np.genfromtxt(timing_path,skip_header=1,usecols=2,dtype='str')

    if onoff[-1]=='OFF':
        cal_times = np.append(cal_times, time.time())

    cal_status = [[0, cal_times[0]]] # calibration source was not during initial states
    cal_status += [[cal_times[i],cal_times[i+1]] for i,s in enumerate(onoff) if s=='OFF']

    # return a list of pairs of [off_time, on time]
    return cal_status


def make_cal_aware_esf(orig_esf_path, new_esf_path, cal_status_orig, off_expand_s):
    # I don't like this code, it's hard to reasonabout the correctness.
    # make alternate experiment state file
    # i have two sources of timestamps
    # cal_status: list of [off_time, on_time] pairs
    # experiment_state_file: allLabels and unixnanos form pairs of [state_name, on_time]
    # experiment_state_file contains state names like A, B, C, D, IGNORE
    # here we merge them into a joint list of new states with names like
    # A_ON, A_OFF, B_ON, B_OFF, etc
    esf = mass.off.ExperimentStateFile(orig_esf_path)
    i=0 #index into cal_status
    j=0 #index into esf.allLabels and esf.unixnanos
    # the loop will walk through time starting from 0, and create states as it goes
    # keeping track of if the cal is on or off at any given time
    cal_on = False
    current_orig_state_name = None
    cal_status = cal_status_orig + [[np.inf, np.inf]]
    new_states = [] # pairs of [state_name, start_nano]
    while j<len(esf.allLabels):
        cal_off_time_nano, cal_on_time_nano = np.array(cal_status[i])*1e9
        cal_on_time_nano-=off_expand_s*1e9
        cal_off_time_nano+=off_expand_s*1e9
        next_state_name, next_state_start_nano = esf.allLabels[j], esf.unixnanos[j]
        if cal_on == True:
            # look for either cal_off or next state
            if cal_off_time_nano < next_state_start_nano:
                cal_on = False
                new_states.append((f"{current_orig_state_name}_OFF", cal_off_time_nano))
            elif next_state_start_nano < cal_off_time_nano:
                new_states.append((f"{next_state_name}_ON", next_state_start_nano))
                current_orig_state_name = next_state_name
                j+=1
            else: # they're equal, probalby wont happen
                raise Exception("not yet implemented")
        else: # cal_on == False
            # look for either cal_off or next state
            if cal_on_time_nano < next_state_start_nano:
                cal_on = True
                new_states.append((f"{current_orig_state_name}_ON", cal_on_time_nano))
                i+=1
            elif next_state_start_nano < cal_on_time_nano:
                new_states.append((f"{next_state_name}_OFF", next_state_start_nano))
                current_orig_state_name = next_state_name
                j+=1
            else: # they're equal, do both
                cal_on = True
                new_states.append((f"{next_state_name}_ON", cal_on_time_nano))
                current_orig_state_name = next_state_name
                i+=1
                j+=1


    with open(new_esf_path,"w") as f:
        f.write("# header\n")
        for (state_name, start_nano) in new_states:
            f.write(f"{int(start_nano)}, {state_name}\n")

import mass
from mass import algorithms

def ds_learnCalibrationPlanFromEnergiesAndPeaks(self, attr, states, ph_fwhm, line_names, maxacc):
    peak_ph_vals, _peak_heights = algorithms.find_local_maxima(self.getAttr(attr, indsOrStates=states), ph_fwhm)
    _name_e, energies_out, opt_assignments = algorithms.find_opt_assignment(peak_ph_vals, line_names, maxacc=maxacc)

    self.calibrationPlanInit(attr)
    for ph, name in zip(opt_assignments, _name_e):
        self.calibrationPlanAddPoint(ph, name, states=states)
mass.off.Channel.learnCalibrationPlanFromEnergiesAndPeaks = ds_learnCalibrationPlanFromEnergiesAndPeaks