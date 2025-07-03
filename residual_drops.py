import matplotlib.pyplot as plt 
import numpy as np 
from lmfit import Parameters, minimize, fit_report
from lmfit.models import VoigtModel, ConstantModel

def midpoints(x):
    return (x[:-1]+x[1:])/2

start , end = 2950, 2980
x = np.arange(start, end, 1)

filename = '20240419_cal'
date = filename[:-4]
extensions = ['nothing', 'Mg', 'Al', 'Cl', 'Sc' ,'V']

e = []
de = []
for ext in extensions:
    data = [0,0]
    for state in ["B_OFF", "C_OFF"]:
        data = np.vstack((data,np.genfromtxt(f"{filename}/{date}_{state}_{ext}.csv", delimiter=',')))
    dat, bins = np.histogram(data,x)
    dat = dat/np.max(dat)  # Normalize the histogram
    plt.plot(midpoints(bins), dat, label=f'{ext} dropped')
    mod = VoigtModel(prefix='v_')
    params = mod.make_params()
    params[f'v_center'].set(min=start, max=end, vary=True)
    params[f'v_amplitude'].set(value=np.max(dat), min=0, vary=True)
    params[f'v_sigma'].set(value=5, min=0, vary=True)
    params[f'v_gamma'].set(min=0, vary=True)
    result = mod.fit(dat, params, x=midpoints(bins))
    plt.plot(midpoints(bins), result.best_fit, label=f'{ext} fit')
    e.append(result.params['v_center'].value)
    de.append(result.params['v_center'].stderr)

ref_e = e[0]
e = np.array(e[1:]) - ref_e
print(np.std(e))
ref_de = de[0]
de = np.sqrt(np.array(de[1:])** 2 + ref_de ** 2)

plt.legend()
plt.show()

plt.errorbar(extensions[1:], e, yerr=de, fmt='o', label='Residuals')
plt.show()