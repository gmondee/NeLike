import mass
import numpy as np

LORENTZIAN_PEAK_HEIGHT = 999
LORENTZIAN_INTEGRAL_INTENSITY = 9999
VOIGT_PEAK_HEIGHT = 99999

# ### Ge K Alpha feature, based on the Zn K Alpha feature shape.
mass.calibration.fluorescence_lines.addline(
    element="Ge",
    material="metal",
    linetype="KAlphaCustomOld",
    reference_short="Zn Hack",
    reference_plot_instrument_gaussian_fwhm=None,
    nominal_peak_energy=9886.52,
    energies=[9886.52,  9855.42],
    lorentzian_fwhm=np.array((3, 3)) * 1.1,
    reference_amplitude=np.array((10, 5.3)),
    reference_amplitude_type=LORENTZIAN_PEAK_HEIGHT,
    ka12_energy_diff=31.1,
)

mass.calibration.fluorescence_lines.addline(
    element="Ge",
    material="metal",
    linetype="KAlphaCustom",
    reference_short="Zn Hack",
    reference_plot_instrument_gaussian_fwhm=2.5,
    nominal_peak_energy=9886.5,#9886.5,#9855.3,
    energies=[9886.47,  9882.68, 9855.32, 9852.73],
    lorentzian_fwhm=np.array((3.024, 3.68, 3.008, 3.26)),
    reference_amplitude=np.array((100, 3.39, 50.54, 3.18)),
    reference_amplitude_type=LORENTZIAN_PEAK_HEIGHT,
    ka12_energy_diff=9886.47-9855.32,
    position_uncertainty=.25
)

mass.calibration.fluorescence_lines.addline(
    element="Ge",
    material="metal",
    linetype="KBetaCustom",
    reference_short="Dean 2020",
    reference_plot_instrument_gaussian_fwhm=2.5,
    nominal_peak_energy=2/3*(10982.169)+1/3*(10977.961),
    energies=np.array((10982.169, 10977.961, 10975.24, 10988.193, 11099.732)),  # Table 3 C_i
    lorentzian_fwhm=np.array((4.2, 4.22, 4.2, 4.2, 4.2, 4.22)),  # Table 3 W_i
    reference_amplitude=np.array((100, 44.9, 11.9, 8.08, 3.03)),  # Table 3 Fraction
    reference_amplitude_type=LORENTZIAN_INTEGRAL_INTENSITY,
    position_uncertainty=.1
)

# mass.calibration.fluorescence_lines.addline(
#     element="Ge",
#     material="metal",
#     linetype="KAlphaCustom",
#     reference_short="Zn Hack",
#     reference_plot_instrument_gaussian_fwhm=None,
#     nominal_peak_energy=9886.52,
#     energies=[9886.52, 9886.52, 9855.42, 9855.42],
#     lorentzian_fwhm=np.array((2.285, 3.358, 2.667, 3.571)) * 1.1,
#     reference_amplitude=np.array((1000, 1, 500, 1)),
#     reference_amplitude_type=LORENTZIAN_PEAK_HEIGHT,
#     ka12_energy_diff=31.1,
# )

# mass.calibration.fluorescence_lines.addline(
#     element="Test",
#     material="metal",
#     linetype="KAlphaCustom",
#     reference_short="Zn Hack",
#     reference_plot_instrument_gaussian_fwhm=None,
#     nominal_peak_energy=9886.52,
#     energies=[9886.52, 9855.42],
#     lorentzian_fwhm=np.array((3.358, 3.571)) * 1.1,
#     reference_amplitude=np.array((1000,  500)),
#     reference_amplitude_type=LORENTZIAN_PEAK_HEIGHT,
#     ka12_energy_diff=31.1,
# )


# W lines from beiersdorfer 2012
mass.calibration.fluorescence_lines.addline(
    element="W",
    material="metal",
    linetype="M2",
    reference_short="Zn Hack", #Beiersdorfer 2012",
    # fitter_type=mass.calibration.line_models.GenericLineModel,
    reference_plot_instrument_gaussian_fwhm=0.5,
    nominal_peak_energy=8299.22,
    energies=np.array((8299.22,)),
    lorentzian_fwhm=np.array((1,)),  # Table 3 W_i
    reference_amplitude=np.array((1,)),  # Table 3 Fraction
    reference_amplitude_type=LORENTZIAN_INTEGRAL_INTENSITY,
    position_uncertainty=.4
)

mass.calibration.fluorescence_lines.addline(
    element="W",
    material="metal",
    linetype="3G",
    reference_short="Zn Hack", #Beiersdorfer 2012",
    reference_plot_instrument_gaussian_fwhm=0.5,
    nominal_peak_energy=8307.51,
    energies=np.array((8307.51,)),
    lorentzian_fwhm=np.array((1,)),  # Table 3 W_i
    reference_amplitude=np.array((1,)),  # Table 3 Fraction
    reference_amplitude_type=LORENTZIAN_INTEGRAL_INTENSITY,
    position_uncertainty=.4
)

mass.calibration.fluorescence_lines.addline(
    element="W",
    material="metal",
    linetype="3D",
    reference_short="Zn Hack", #Beiersdorfer 2012",
    reference_plot_instrument_gaussian_fwhm=0.5,
    nominal_peak_energy=9126.25,
    energies=np.array((9126.25,)),
    lorentzian_fwhm=np.array((1,)),  # Table 3 W_i
    reference_amplitude=np.array((1,)),  # Table 3 Fraction
    reference_amplitude_type=LORENTZIAN_INTEGRAL_INTENSITY,
    position_uncertainty=.5
)

mass.calibration.fluorescence_lines.addline( #1:1.5 intensity
    element="W",
    material="metal",
    linetype="M2+3G_23",
    reference_short="Dean 2020",
    reference_plot_instrument_gaussian_fwhm=0.5,
    nominal_peak_energy=1/2*(8299.22)+1/2*(8307.51),
    energies=np.array((8299.22, 8307.51)), 
    lorentzian_fwhm=np.array((4,4)), 
    reference_amplitude=np.array((100,150)),  # Table 3 Fraction
    reference_amplitude_type=LORENTZIAN_INTEGRAL_INTENSITY,
    position_uncertainty=.5
)

mass.calibration.fluorescence_lines.addline( #1:1 intensity
    element="W",
    material="metal",
    linetype="M2+3G_25",
    reference_short="Dean 2020",
    reference_plot_instrument_gaussian_fwhm=0.5,
    nominal_peak_energy=1/2*(8299.22)+1/2*(8307.51),
    energies=np.array((8299.22, 8307.51)), 
    lorentzian_fwhm=np.array((4,4)), 
    reference_amplitude=np.array((100,100)),  # Table 3 Fraction
    reference_amplitude_type=LORENTZIAN_INTEGRAL_INTENSITY,
    position_uncertainty=.5
)
# element=element,
# material="Highly Charged Ion",
# linetype=linetype,
# reference_short='NIST ASD',
# fitter_type=line_models.GenericLineModel,
# reference_plot_instrument_gaussian_fwhm=0.5,
# nominal_peak_energy=nominal_peak_energy,
# energies=energies,
# lorentzian_fwhm=widths,
# reference_amplitude=ratios,
# reference_amplitude_type=LORENTZIAN_PEAK_HEIGHT,
# ka12_energy_diff=None