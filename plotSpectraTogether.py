import matplotlib as mpl
mpl.use('QtAgg')
mpl.rcParams["image.interpolation"] = "none"
import matplotlib.pyplot as plt
plt.ion()
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import lmfit
import pathlib

spectraPath = pathlib.Path(r'/Users/gmondeel/Documents/Mass/NeLike/Spectra')
spectraFiles = list(spectraPath.glob("*.txt"))

numSpectra = len(spectraFiles)
fig = plt.figure()
ax = plt.gca()
#Bin center (eV),Counts per 2 eV bin
verticalBuffer = 1000
startInd = 2500
#xs = np.arange(800, 13000, 2.)
for i, spectrumFile in enumerate(spectraFiles):
  spectrum= pd.read_csv(spectrumFile, delimiter=",")
  xs = np.array(spectrum["Bin center (eV)"].to_list())[startInd:]
  counts = np.array(spectrum["Counts per 2 eV bin"].to_list())[startInd:]
  countsNorm = counts*(verticalBuffer/max(counts))
  ax.plot(xs, countsNorm-i*verticalBuffer, label=str(spectraFiles[i]).split('/')[-1][:-4])
plt.legend()

fig = plt.figure()
ax = plt.gca()
#Bin center (eV),Counts per 2 eV bin
verticalBuffer = 1000
startInd = 2500
#xs = np.arange(800, 13000, 2.)
cmap=plt.get_cmap("brg", len(spectraFiles))
for i, spectrumFile in enumerate(spectraFiles):
  if "Os" in str(spectrumFile):
    elemScale=1
  elif "Re" in str(spectrumFile):
    elemScale=2
  elif "W" in str(spectrumFile):
    elemScale=3
  elif "Ir" in str(spectrumFile):
    elemScale=0
  
  color = cmap(int(str(spectrumFile)[-5])-3) #0,1,2


  spectrum= pd.read_csv(spectrumFile, delimiter=",")
  xs = np.array(spectrum["Bin center (eV)"].to_list())[startInd:]
  counts = np.array(spectrum["Counts per 2 eV bin"].to_list())[startInd:]
  countsNorm = counts*(verticalBuffer/max(counts))
  ax.plot(xs, countsNorm-elemScale*verticalBuffer, label=str(spectraFiles[i]).split('/')[-1][:-4], color=color)
plt.legend()
