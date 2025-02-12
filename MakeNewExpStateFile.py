from local_util import make_cal_aware_esf, load_timing
import mass
from mass.off import ChannelGroup, Channel, getOffFileListFromOneFile

#Laptop
#filename = r"C:\Users\Grant Mondeel\Box\CfA\TES\Ne-Like\20240723_off\0000\20240723_run0000_chan1.off"

#PC
filename = r"C:\Users\lamat\Box\CfA\TES\Ne-Like\20240723_off\0000\20240723_run0000_chan1.off"


data = ChannelGroup(getOffFileListFromOneFile(filename, maxChans=999), verbose=True, channelClass = Channel)

###Separate states into calibration source on/off
timing_file = r"C:\Users\lamat\Box\CfA\TES\Ne-Like\20240723_off\0000\time_20240723.txt"
cal_status = load_timing(timing_file)
make_cal_aware_esf(orig_esf_path=data.experimentStateFile.filename, 
                   new_esf_path=r"C:\Users\lamat\Box\CfA\TES\Ne-Like\20240723_off\0000\20240723_run0000_experiment_state_new.txt", 
                   cal_status_orig=cal_status, off_expand_s=3)

