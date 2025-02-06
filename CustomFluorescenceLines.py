import mass
import numpy as np

LORENTZIAN_PEAK_HEIGHT = 999
LORENTZIAN_INTEGRAL_INTENSITY = 9999
VOIGT_PEAK_HEIGHT = 99999

### Ge K Alpha feature, based on the Zn K Alpha feature shape.
mass.calibration.fluorescence_lines.addline(
    element="Ge",
    material="metal",
    linetype="KAlphaCustom",
    reference_short="Zn Hack",
    reference_plot_instrument_gaussian_fwhm=None,
    nominal_peak_energy=9886.52,
    energies=[9886.52,  9855.42],
    lorentzian_fwhm=np.array((3, 3)) * 1.1,
    reference_amplitude=np.array((10, 5.3)),
    reference_amplitude_type=LORENTZIAN_PEAK_HEIGHT,
    ka12_energy_diff=31.1,
)