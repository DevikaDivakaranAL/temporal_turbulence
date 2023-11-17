from ComputeWindSpeed import *
from ComputeRefractiveIndexStructureParameter import *
from ComputeTurbulenceStrength import *
from ComputeBeamEffects import *
from ComputeScintillationIndex import *
from ComputeProbabilityDensityFunction import *
from ComputeIntensityTimeSeries import *
from ComputeIntensityParameter import *
from ComputeTimeParameter import *
import matplotlib.pyplot as plt


def main_fade_simulation(link_geometry, daynightmodel, wind_height_min, wind_height_max, sunset, sunrise, time, divergence, M2, wavelength, fadeProb, sltAltitude, altAltitude, elevationAngle, ground_wind=8, slew_rate=0.1, hv_ground_cst= 1.7e-14, lower_limit = 0.004, upper_limit = 1.6, transmission_losses = -2, altApertureDiameter = 0.3, printResults = False, C_r = 0, integration_step_multiplier = 1, nr_transmitters = 1, compute_only_fades = True):
    """Main function of the simulation. Creates object of each class involved in the simulation and calls their main method in order to compute all their attributes. It then return the values of interest.
    If requested by the user can print a table of the results.
    :param ground_wind: float,   wind speed at ground in m/s
    :param slew_rate: float, Slew rate in deg/s
    :param daynightmodel: float, day or night
    :param sunset: sunset time: hr.min format
    :param sunrise: sunrise time: hr.min format
    :param time:  time at which hap model is used in night time
    :param link_geometry: string,   can be 'uplink' or 'downlink'
    :param divergence: float, the beam divergence in rad
    :param M2: float, dimensionless parameter that quantifies the laser beam quality
    :param wavelength: float, beam wavelength in m
    :param wind_height_min: minimum height at which wind speed is calculated , in m
    :param wind_height_max: maximum height at which wind speed is calculated , in m
    :param fadeProb: float, fade statistics probability of interest.
    :param sltAltitude: float, SLT altitude in m
    :param altAltitude: float, ALT altitude in m
    :param elevationAngle: elevation angle in degrees
    :param rms_wind_speed: float, rms wind speed computed from wind profile (in m/s)
    :param hv_ground_cst: float,    Hufnagel Valley constant A_0 in m**(-2/3)
    :param altApertureDiameter: float, ALT aperture diameter in m
    :param lower_limit: float, lower limit of the turbulence eddies (in m)
    :param upper_limit: float,    upper limit of the turbulent eddies (in m)
    :param transmission_losses: the transmission losses at the receiver (in dB)
    :param printResults: boolean, True if user wants a table of results to be printed, False otherwise
    :param C_r: float, scaling constant related to the size of biggest turbulence eddies
    :param integration_step_multiplier: float, factor that multiplies the integration step. Can be increased to make computation faster (or decreased for more precision)
    :param nr_transmitters: int, number of independent trasmitters N, by default = 1.
    :param compute_only_fades:  boolean, True if user wants also surge statitistics to be computed, False otherwise
    :return: the losses in dB due to fades with respect to mean irradiance considering turbulence, the losses in dB due to change of mean irradiance from diffractive only to turbulent case
    """

    wavelength = wavelength*1e-9                                     #Simulation uses SI units for all parameters
    W_0 =  M2 * wavelength / (np.pi * divergence)                   #Computes the beam waist at the receiver from divergence and Beam Propagation Ratio
    F_0 = np.inf                                                    #Lasers are collimated beams, the phase front radius of curvature at the transmitter is always infinite

    # Initialize models
    windModel = ComputeWindSpeed(Vg=ground_wind, slew=slew_rate, wind_height_min=wind_height_min, wind_height_max=wind_height_max, geometry=link_geometry)
    #windModel.compute_wind_speed()
    c2nModel = ComputeRefractiveIndexStructureParameter(windModel, daynightmodel, wind_height_max=wind_height_max, sunset=sunset, sunrise=sunrise, time=time, hv_ground_cst=hv_ground_cst)
    c2nModel.compute_c2n_fct()
    turbulenceStrengthModel = ComputeTurbulenceStrength(c2nModel, elevation=elevationAngle, geometry=link_geometry,
                                                        ALT_altitude=altAltitude, SLT_altitude=sltAltitude,
                                                        wavelength=wavelength, W_0=W_0, F_0=F_0)
    turbulenceStrengthModel.compute_turbulence_strength()
    beamEffectsModel = ComputeBeamEffects(c2nModel, turbulenceStrengthModel)
    beamEffectsModel.compute_beam_effects()
    scintillationIndexModel = ComputeScintillationIndex(c2nModel, turbulenceStrengthModel, beamEffectsModel, C_r=C_r,
                                                        D=altApertureDiameter)
    scintillationIndexModel.compute_scintillation_index()
    intensityTimeSeriesModel = ComputeIntensityTimeSeries(c2nModel, scintillationIndexModel, windModel,
                                                          turbulenceStrengthModel, lower_limit, upper_limit,
                                                          transmission_losses, altApertureDiameter)
    intensityTimeSeriesModel.compute_intensity_time_series()
    pdfModel = ComputeProbabilityDensityFunction(turbulenceStrengthModel, beamEffectsModel,
                                                 scintillationIndexModel, intensityTimeSeriesModel)  # Instantiate a ComputeProbabilityDensityFunction object
    pdfModel.compute_pdf()  # Compute all pdfModel object attributes

    intensityParametersModel = ComputeIntensityParameters(turbulenceStrengthModel, scintillationIndexModel, pdfModel,
                                                          probability=fadeProb,
                                                          integration_step_multiplier=integration_step_multiplier,
                                                          nr_transmitters=nr_transmitters,
                                                          compute_only_fades=compute_only_fades)
    # Instantiate a ComputeIntensityParameters object
    intensityParametersModel.compute_intensity_parameters()
    timeParametersModel = ComputeTimeParameters(c2nModel, turbulenceStrengthModel, scintillationIndexModel, pdfModel,
                                                intensityParametersModel)
    # Instantiate a ComputeTimeParameters object
    timeParametersModel.compute_time_parameters()

    # Perform computations

    if printResults == True:
        print('------------------------------------------ Results ------------------------------------------')
        print('Link geometry:                                         {}'.format(turbulenceStrengthModel.geometry))
        print("Gateway altitude:                                      {:.2e} m".format(
            turbulenceStrengthModel.ALT_altitude))
        print("Satellite altitude:                                    {:.2e} m".format(
            turbulenceStrengthModel.SLT_altitude))
        print('Link length:                                           {:.2e} m'.format(turbulenceStrengthModel.R))
        print('Wavelength:                                            {:.2e} m'.format(
            turbulenceStrengthModel.wavelength))
        print('Turbulence regime:                                     {}'.format(
            turbulenceStrengthModel.turbulence_strength))
        print('Rytov variance   :                                     {:.2e}'.format(
            turbulenceStrengthModel.rytov_variance))
        print('Beam size at transmitter (W_0):                        {:.2e} m'.format(turbulenceStrengthModel.W_0))
        print('Beam size at receiver (diffractive):                   {:.2e} m'.format(turbulenceStrengthModel.W))
        print('Beam size at receiver (with turbulence):               {:.2e} m'.format(beamEffectsModel.W_eff))
        print('Coherence length (Fried parameter):                    {:.2e} m'.format(turbulenceStrengthModel.r_0))
        print('Beam wander r_c**2:                                           {} m'.format(
            "None" if beamEffectsModel.r2_c is None else str('{:.2e}'.format(beamEffectsModel.r2_c))))
        print("Scintillation index (untracked if uplink):             {:.2e}".format(
            scintillationIndexModel.scintillation_index))
        print('Average intensity at receiver:                         {:.2e} W/m**2'.format(pdfModel.mean_intensity))
        print('Scintillation loss (fade)                               {} dB'.format(
            "None" if intensityParametersModel.F_T is None else str('-{:.2f}'.format(intensityParametersModel.F_T))))
        print('Probability of such fade loss                          {}'.format(
            "None" if intensityParametersModel.fade_probability is None else str(
                '{:.4f}'.format(intensityParametersModel.fade_probability))))
        print('Avg fades per second at {}dB:                        {}'.format(
            "None" if intensityParametersModel.F_T is None else str('{:.2f}'.format(intensityParametersModel.F_T)),
            "None" if timeParametersModel.avg_number_fades is None else str(
                '{:.2e}'.format(timeParametersModel.avg_number_fades))))
        print('Avg fade time at {} dB:                              {} s'.format(
            "None" if intensityParametersModel.F_T is None else str('{:.2f}'.format(intensityParametersModel.F_T)),
            "None" if timeParametersModel.avg_fade_duration is None else str(
                '{:.2e}'.format(timeParametersModel.avg_fade_duration))))
        print('Surge level                                            {} dB'.format(
            "None" if intensityParametersModel.S_T is None else str('{:.2f}'.format(intensityParametersModel.S_T))))
        print('Probability of such surge gain                          {}'.format(
            "None" if intensityParametersModel.surge_probability is None else str(
                '{:.4f}'.format(intensityParametersModel.surge_probability))))
        print('Avg surges per second at {}dB:                       {}'.format(
            "None" if intensityParametersModel.S_T is None else str('{:.2f}'.format(intensityParametersModel.S_T)),
            "None" if timeParametersModel.avg_number_surges is None else str(
                '{:.2e}'.format(timeParametersModel.avg_number_surges))))
        print('Avg surge time at {} dB (aperture averaged):         {} s'.format(
            "None" if intensityParametersModel.S_T is None else str('{:.2f}'.format(intensityParametersModel.S_T)),
            "None" if timeParametersModel.avg_surges_duration is None else str(
                '{:.2e}'.format(timeParametersModel.avg_surges_duration))))

        if link_geometry == 'downlink':
            print(
                '------------------------------------------ Downlink only ---------------------------------------------')
            print("Receiver aperture                                      {:.2e} m".format(scintillationIndexModel.D))
            print('Scintillation index aperture averaged:                 {}'.format(
                "None" if scintillationIndexModel.scintillation_index_averaged is None else str(
                    '{:.2e}'.format(scintillationIndexModel.scintillation_index_averaged))))
            print('Scintillation loss aperture averaged(fade)            {} dB'.format(
                "None" if intensityParametersModel.F_T_averaged is None else str(
                    '-{:.2f}'.format(intensityParametersModel.F_T_averaged))))
            print('Probability of such fade loss                          {}'.format(
                "None" if intensityParametersModel.fade_probability_averaged is None else str(
                    '{:.4f}'.format(intensityParametersModel.fade_probability_averaged))))
            print('Avg fades per second at {}dB (aperture averaged):    {}'.format(
                "None" if intensityParametersModel.F_T_averaged is None else str(
                    '{:.2f}'.format(intensityParametersModel.F_T_averaged)),
                "None" if timeParametersModel.avg_number_fades_aperture is None else str(
                    '{:.2e}'.format(timeParametersModel.avg_number_fades_aperture))))
            print('Avg fade time at {} dB (aperture averaged):          {} s'.format(
                "None" if intensityParametersModel.F_T_averaged is None else str(
                    '{:.2f}'.format(intensityParametersModel.F_T_averaged)),
                "None" if timeParametersModel.avg_fade_duration_aperture is None else str(
                    '{:.2e}'.format(timeParametersModel.avg_fade_duration_aperture))))
            print('Surge gains, aperture averaged                         {} dB'.format(
                "None" if intensityParametersModel.S_T_averaged is None else str(
                    '{:.2f}'.format(intensityParametersModel.S_T_averaged))))
            print('Probability of such surge gain                          {}'.format(
                "None" if intensityParametersModel.surge_probability_averaged is None else str(
                    '{:.4f}'.format(intensityParametersModel.surge_probability_averaged))))
            print('Avg surges per second at {}dB (aperture averaged):   {}'.format(
                "None" if intensityParametersModel.S_T_averaged is None else str(
                    '{:.2f}'.format(intensityParametersModel.S_T_averaged)),
                "None" if timeParametersModel.avg_number_surges_aperture is None else str(
                    '{:.2e}'.format(timeParametersModel.avg_number_surges_aperture))))
            print('Avg surge time at {} dB (aperture averaged):         {} s'.format(
                "None" if intensityParametersModel.S_T_averaged is None else str(
                    '{:.2f}'.format(intensityParametersModel.S_T_averaged)),
                "None" if timeParametersModel.avg_surges_duration_aperture is None else str(
                    '{:.2e}'.format(timeParametersModel.avg_surges_duration_aperture))))

        elif link_geometry == 'uplink':
            print(
                '------------------------------------------ Uplink only ---------------------------------------------')
            print('Scintillation index tracked                            {}'.format(
                "None" if scintillationIndexModel.scintillation_index_tracked is None else str(
                    '{:.2e}'.format(scintillationIndexModel.scintillation_index_tracked))))
            print('Scintillation loss tracked                            {} dB'.format(
                "None" if intensityParametersModel.F_T_tracked is None else str(
                    '-{:.2f}'.format(intensityParametersModel.F_T_tracked))))
            print('Probability of such fade loss                          {}'.format(
                "None" if intensityParametersModel.fade_probability_tracked is None else str(
                    '{:.4f}'.format(intensityParametersModel.fade_probability_tracked))))
            print('Avg fades per second at {} dB (tracked):             {}'.format(
                "None" if intensityParametersModel.F_T_tracked is None else '{:.2f}'.format(
                    intensityParametersModel.F_T_tracked),
                "None" if timeParametersModel.avg_number_fades_tracked is None else str(
                    '{:.2e}'.format(timeParametersModel.avg_number_fades_tracked))))
            print('Avg fade time at {} dB  (tracked) :                  {} s'.format(
                "None" if intensityParametersModel.F_T_tracked is None else str(
                    '{:.2f}'.format(intensityParametersModel.F_T_tracked)),
                "None" if timeParametersModel.avg_fade_duration_tracked is None else str(
                    '{:.2e}'.format(timeParametersModel.avg_fade_duration_tracked))))
            print('Scintillation gains tracked (surges)                   {} dB'.format(
                "None" if intensityParametersModel.S_T_tracked is None else str(
                    '{:.2f}'.format(intensityParametersModel.S_T_tracked))))
            print('Probability of such surge gains                       {}'.format(
                "None" if intensityParametersModel.surge_probability_tracked is None else str(
                    '{:.4f}'.format(intensityParametersModel.surge_probability_tracked))))
            print('Avg surges per second at {} dB (tracked):            {}'.format(
                "None" if intensityParametersModel.S_T_tracked is None else '{:.2f}'.format(
                    intensityParametersModel.S_T_tracked),
                "None" if timeParametersModel.avg_number_surges_tracked is None else str(
                    '{:.2e}'.format(timeParametersModel.avg_number_surges_tracked))))
            print('Avg surge time at {} dB  (tracked) :                {} s'.format(
                "None" if intensityParametersModel.S_T_tracked is None else str(
                    '{:.2f}'.format(intensityParametersModel.S_T_tracked)),
                "None" if timeParametersModel.avg_surges_duration_tracked is None else str(
                    '{:.2e}'.format(timeParametersModel.avg_surges_duration_tracked))))

        if nr_transmitters > 1:
            print(
                '------------------------------------------ Multiple transmitters only ---------------------------------------------')
            print('Scintillation loss (fade) with {} transmitters         -{:.2f} dB'.format(
                intensityParametersModel.nr_transmitters, intensityParametersModel.F_T_multiple))
            print('Probability of such fade loss                          {:.4f}'.format(
                intensityParametersModel.fade_probability_multiple))
            print('Scintillation gains (surges) with {} transmitters       {:.2f} dB'.format(
                intensityParametersModel.nr_transmitters, intensityParametersModel.S_T_multiple))
            print('Probability of such surge gains                        {:.4f}'.format(
                intensityParametersModel.surge_probability_multiple))

        print('---------------------------------------------------------------------------------------')
        print('Change of mean loss                                     {:.4f}'.format(pdfModel.mean_change_loss))

    scintillation_loss = None
    if link_geometry == 'downlink':
        if turbulenceStrengthModel.turbulence_strength == 'weak':
            scintillation_loss = intensityParametersModel.F_T_averaged
        else:
            scintillation_loss = intensityParametersModel.F_T
    else:
        scintillation_loss = intensityParametersModel.F_T

    if link_geometry == 'uplink':
        scintillation_surge = intensityParametersModel.S_T_tracked
    else:
        scintillation_surge = intensityParametersModel.S_T_averaged

    if nr_transmitters > 1:
        scintillation_loss = intensityParametersModel.F_T_multiple
        scintillation_surge = intensityParametersModel.S_T_multiple

    return pdfModel.mean_change_loss, scintillation_loss, scintillation_surge, turbulenceStrengthModel.elevation, windModel.rms_wind_speed, c2nModel.c2n_value, turbulenceStrengthModel.rytov_variance, scintillationIndexModel.scintillation_index




