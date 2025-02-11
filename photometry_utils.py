import numpy as np
from rubin_sim.phot_utils import signaltonoise
from pyLIMA import event, toolbox

def flux_to_mag(ZP, flux):
    """
    Function converts input flux value into magnitudes, using the zeropoint provided
    """
    return ZP - 2.5*np.log10(flux)

def mag_to_flux(ZP, mag):
    """
    Function converts input flux value into magnitudes, using the zeropoint provided
    """
    return 10**((ZP - mag)/2.5)

def generate_event_photometry_arrays(time_stamps, target_mag, mag_offset, ZP, dataSlice, LSST_BandPass, photParams):
    """
    Set up data arrays for lightcurves in all passbands.  At this point we just have the timestamps,
    and the photometry will be added later, so we fill the magnitude and error columns with the limiting magnitude values
    """

    photometry = {}
    baseline_flux = {}
    flux_parameters = {}

    for f, ts in time_stamps.items():
        # Timestamps from the LSST opsim
        arr = np.zeros((len(ts), 3))
        arr[:, 0] = ts

        # Baseline (unmagnified) target photometry, given the simulation parameters
        mag_baseline = target_mag[f] + mag_offset
        arr[:, 1] = mag_baseline

        # Calculate the corresponding source and blend fluxes, assuming no blending
        flux_baseline = mag_to_flux(ZP, mag_baseline)
        flux_parameters[f] = {'f_source': flux_baseline, 'f_blend': 0.0}
        baseline_flux[f] = flux_baseline

        # Photometric uncertainties
        # These can be calculated using the Rubin_sim package's built in functions, using the m5 value for the first image
        m5 = dataSlice['fiveSigmaDepth'][np.where(dataSlice['filter'] == f)]
        magerr = signaltonoise.calc_mag_error_m5(mag_baseline, LSST_BandPass[f], m5[0], photParams)[0]
        arr[:, 2] = toolbox.brightness_transformation.error_magnitude_to_error_flux(magerr, flux_baseline)

        # Store photometry
        photometry[f] = arr

    return photometry, baseline_flux, flux_parameters


def update_phot_uncertainties(event_model, dataSlice, filters, LSST_BandPass, photParams, ZP):
    for j, tel in enumerate(event_model.event.telescopes):
        flux = tel.lightcurve_flux['flux']
        mag = flux_to_mag(ZP, flux)

        m5 = dataSlice['fiveSigmaDepth'][np.where(dataSlice['filter'] == filters[j])]
        magerr = [signaltonoise.calc_mag_error_m5(m, LSST_BandPass[filters[j]], m5[k], photParams)[0] for k, m in enumerate(mag)]

        fluxerr = toolbox.brightness_transformation.error_magnitude_to_error_flux(magerr, flux)
        tel.lightcurve_flux['err_flux'] = fluxerr
        tel.lightcurve_magnitude['err_mag'] = magerr

    return event_model