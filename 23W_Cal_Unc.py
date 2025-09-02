from calibrators import cals
import matplotlib.pyplot as plt
import numpy as np

### Determine which lines to include in the calibration uncertainty estimate
#These are the W64+ lines from Beiersdorfer 2012; there are also lines from other charges states
W64LitCentroids = [8299.22, 8307.51, 8450.14, 8996.31, 9126.25, 9689.29, 10317.23, 10408.69, 10706.85, 10967.75]
W64LitUncs = [0.40, 0.40, 0.60, 0.50, 0.50, 0.50, 0.50, 0.40, 0.90, 1.20]
#W63+ lines
W63LitCentroids = [9088.8, 9096.0, 10350.1, 10367.3, 10381.2]
W63LitUncs = [0.5, 0.7, 1.0, 1.0, 1.0]
#W62+ lines
W62LitCentroids = [9055.3, 10350.1]
W62LitUncs = [0.6,1.2]
#W65+ lines
W65LitCentroids = [8455.4, 8466.0, 8519.7, 8592.7, 9256.3, 9256.3, 9291.5, 10562.2, 10562.2, 10562.2]
W65LitUncs = [0.6, 0.6, 0.6, 0.7, 0.6, 0.6, 1.0, 0.5, 0.5, 0.5]
#W66+ lines
W66LitCentroids = [8592.7, 8599.4, 8607.1, 8644.5, 8654.4, 9413.6] 
W66LitUncs = [0.7, 1.0, 1.0, 0.6, 0.6, 1.2]

### calArr should look like [literature E, lit Unc, experiment E, exp Unc]


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
