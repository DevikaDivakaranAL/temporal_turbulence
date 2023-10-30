from ComputeWindSpeed import *
from ComputeRefractiveIndexStructureParameter import *
from ComputeTurbulenceStrength import *
from ComputeBeamEffects import *
from ComputeScintillationIndex import *
# from ComputeProbabilityDensityFunction import *
# from ComputeIntensityParameters import *
# from ComputeTimeParameters import *
import matplotlib.pyplot as plt


def main_fade_simulation(link_geometry, HVmodel, heights, divergence, M2, wavelength, fadeProb, sltAltitude, altAltitude,elevationAngle,ground_wind=8, slew_rate=0.1, hv_ground_cst= 1.7e-14, altApertureDiameter = 0.3, printResults = False, C_r = 0, integration_step_multiplier = 1, nr_transmitters = 1, compute_only_fades = True):
    """Main function of the simulation. Creates object of each class involved in the simulation and calls their main method in order to compute all their attributes. It then return the values of interest.
    If requested by the user can print a table of the results.
    :param ground_wind: float,   wind speed at ground in m/s
    :param slew_rate: float, Slew rate in deg/s
    :param C2n_model: float, day or night
    :param link_geometry: string,   can be 'uplink' or 'downlink'
    :param divergence: float, the beam divergence in rad
    :param M2: float, dmensionless parameter that quantifies the laser beam quality
    :param wavelength: float, beam wavelength in m
    :param fadeProb: float, fade statistics probability of interest.
    :param sltAltitude: float, SLT altitude in m
    :param altAltitude: float, ALT altitude in m
    :param elevationAngle: elevation angle in degrees
    :param rms_wind_speed: float, rms wind speed computed from wind profile (in m/s)
    :param hv_ground_cst: float,    Hufnagel Valley constant A_0 in m**(-2/3)
    :param altApertureDiameter: float, ALT aperture diameter in m
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
    windModel = ComputeWindSpeed(Vg=ground_wind, Ws=slew_rate, height=heights, geometry=link_geometry)

    c2nModel = ComputeRefractiveIndexStructureParameter(windModel, HVmodel, hv_ground_cst=hv_ground_cst)

    turbulenceStrengthModel = ComputeTurbulenceStrength(c2nModel, elevation=elevationAngle, geometry=link_geometry,
                                                        ALT_altitude=altAltitude, SLT_altitude=sltAltitude,
                                                        wavelength=wavelength, W_0=W_0, F_0=F_0)
    beamEffectsModel = ComputeBeamEffects(c2nModel, turbulenceStrengthModel)
    scintillationIndexModel = ComputeScintillationIndex(c2nModel, turbulenceStrengthModel, beamEffectsModel, C_r=C_r,
                                                        D=altApertureDiameter)
    # Perform computations
    windModel.compute_wind_speed()
    c2nModel.compute_c2n_fct()
    turbulenceStrengthModel.compute_turbulence_strength()
    beamEffectsModel.compute_beam_effects()
    scintillationIndexModel.compute_scintillation_index()

    # print('------------------------------------------ Results ------------------------------------------')
    # print('Link geometry:                                         {}'.format(turbulenceStrengthModel.geometry))
    # print('Height:                                                {:.2e} m'.format(windModel.height))
    # print('wind speed:                                            {:.2e} m/s'.format(windModel.rms_wind_speed))
    # print('Slew rate:                                             {:.2e} '.format(windModel.Ws))
    # print("Gateway altitude:                                      {:.2e} m".format(turbulenceStrengthModel.ALT_altitude))
    # print("Satellite altitude:                                    {:.2e} m".format(turbulenceStrengthModel.SLT_altitude))
    # print('Link length:                                           {:.2e} m'.format(turbulenceStrengthModel.R))
    # print('Wavelength:                                            {:.2e} m'.format(turbulenceStrengthModel.wavelength))
    print('Turbulence regime:                                     {}'.format(turbulenceStrengthModel.turbulence_strength))
    # print('Rytov variance   :                                     {:.2e}'.format(turbulenceStrengthModel.rytov_variance))
    # print('Beam size at transmitter (W_0):                        {:.2e} m'.format(turbulenceStrengthModel.W_0))
    # print('Beam size at receiver (diffractive):                   {:.2e} m'.format(turbulenceStrengthModel.W))
    # print('Beam size at receiver (with turbulence):               {:.2e} m'.format(beamEffectsModel.W_eff))
    # print('Coherence length (Fried parameter):                    {:.2e} m'.format(turbulenceStrengthModel.r_0))
    # print('Beam wander r_c**2:                                           {} m'.format("None" if beamEffectsModel.r2_c is None else str('{:.2e}'.format(beamEffectsModel.r2_c))))
    # print("Scintillation index (untracked if uplink):             {:.2e}".format( scintillationIndexModel.scintillation_index))
    # if link_geometry == 'downlink':
    #     print('------------------------------------------ Downlink only ---------------------------------------------')
    #     print("Receiver aperture                                      {:.2e} m".format(scintillationIndexModel.D))
    #     print('Scintillation index aperture averaged:                 {}'.format("None" if scintillationIndexModel.scintillation_index_averaged is None else str('{:.2e}'.format(scintillationIndexModel.scintillation_index_averaged))))
    #
    # else:
    #     print('------------------------------------------ Uplink only ---------------------------------------------')
    #     print('Scintillation index tracked                            {}'.format("None" if scintillationIndexModel.scintillation_index_tracked is None else str('{:.2e}'.format(scintillationIndexModel.scintillation_index_tracked))))


    return windModel.height, windModel.rms_wind_speed, c2nModel.c2n, turbulenceStrengthModel.rytov_variance, scintillationIndexModel.scintillation_index

