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
from lmfit.models import GaussianModel
import pathlib

spectraPath = pathlib.Path(r'/Users/gmondeel/Documents/Mass/NeLike')
spectraFiles = list(spectraPath.glob("*2024*.txt"))

def fit_multiple_gaussians(x, y, centroids, sigma_guess=1.0):
    """
    Fit multiple Gaussian models to the data using lmfit.

    Parameters:
        x (np.array): The x-data.
        y (np.array): The y-data.
        centroids (list or np.array): List of initial guesses for the Gaussian centers.
        sigma_guess (float): Initial guess for the Gaussian width (standard deviation).

    Returns:
        result: lmfit ModelResult object from the fit.
        composite_model: the lmfit composite model object.
    """
    model = None
    params = lmfit.Parameters()
    shared_sigma_name='shared_sigma'
    shared_sigma=1.5
    params.add('shared_sigma', value=shared_sigma, min=0.1,max=13/2.35)
    for i, center in enumerate(centroids):
        prefix = f'g{int(centroids[i])}_'
        gauss = GaussianModel(prefix=prefix)
        
        if model is None:
            model = gauss
        else:
            model += gauss

        # Create initial parameter guesses
        amp_guess = y[np.abs(x - center).argmin()]
        param_guesses = gauss.make_params()
        param_guesses[f'{prefix}amplitude'].set(value=amp_guess, min=1)
        param_guesses[f'{prefix}center'].set(min=center-15, max=center+15)
        param_guesses[f'{prefix}sigma'].set(expr=shared_sigma_name)

        params.update(param_guesses)
    lin = lmfit.models.LinearModel(prefix='lin_')
    lin_param_guess = lin.make_params()
    lin_param_guess['lin_slope'].set(value=0, vary=False)
    # lin_param_guess['lin_intercept'].set(min=0, max=3)

    model+=lin
    params.update(lin_param_guess)
    # Fit the composite model
    result = model.fit(y, params, x=x)

    return result, model


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
startInd = 3500
#xs = np.arange(800, 13000, 2.)
cmap=plt.get_cmap("brg", len(spectraFiles))

fits={}
for i, spectrumFile in enumerate(spectraFiles):
  if "Os" in str(spectrumFile):
    elemScale=1
    centroidLit=[]
  elif "Re" in str(spectrumFile):
    elemScale=2
    centroidLit=[]
  elif "W" in str(spectrumFile):
    elemScale=3
    centroidLit=[8299.22, 8307.51, 8450.14, 8996.31, 9126.25, 9689.29, 10317.23, 10408.69, 10706.85, 10967.75]
  elif "Ir" in str(spectrumFile):
    elemScale=0
    centroidLit=[]
  
  color = cmap(int(str(spectrumFile)[-5])-3) #0,1,2


  spectrum= pd.read_csv(spectrumFile, delimiter=",")
  xs = np.array(spectrum["Bin center (eV)"].to_list())[startInd:]
  counts = np.array(spectrum["Counts per 2 eV bin"].to_list())[startInd:]
  countsNorm = counts*(verticalBuffer/max(counts))
  ax.plot(xs, countsNorm-elemScale*verticalBuffer, label=str(spectraFiles[i]).split('/')[-1][:-4], color=color)
  if centroidLit:
    multiGaussFit, multiGaussModel = fit_multiple_gaussians(x=xs, y=counts+3, centroids=centroidLit) #i add 1 to counts because least squares doesnt do well with low counts(galen?)
    fits[str(spectraFiles[i]).split('/')[-1][:-4]]=[multiGaussFit, multiGaussModel]
    # Create a finer x-grid for a smoother fit curve
    x_fine = np.linspace(xs.min(), xs.max(), len(xs)*3)  # increase number of points

    # Evaluate the best-fit model on the fine grid
    y_fit_smooth = multiGaussModel.eval(multiGaussFit.params, x=x_fine)

    # Plot original data and smooth fit
    plt.figure(figsize=(10, 6))
    plt.plot(xs, counts+1, label='Data', alpha=0.6)
    plt.plot(x_fine, y_fit_smooth, label='Smooth Fit', color='red', lw=2)
    for j, c in enumerate(centroidLit):
      plt.axvline(x=c, color='gray', linestyle='--', alpha=0.7)
      plt.text(c, plt.ylim()[1]*0.9, f'{c:.1f}', rotation=45,
              verticalalignment='top', horizontalalignment='center',
              fontsize=9, color='gray')
      center = multiGaussFit.params[f'g{int(centroidLit[j])}_center'].value
      center_err = multiGaussFit.params[f'g{int(centroidLit[j])}_center'].stderr
      label = f"E={center:.2f}±{center_err:.2f}" if center_err else f"E= {center:.2f}"

      plt.text(center, plt.ylim()[1]*0.85,
             label,
             rotation=45,
             verticalalignment='top',
             horizontalalignment='center',
             fontsize=9,
             color='blue')
    plt.legend()
    plt.title(str(spectraFiles[i]).split('/')[-1][:-4])
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid(True)

ax.legend()
