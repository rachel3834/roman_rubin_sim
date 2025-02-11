import numpy as np
from rubin_sim.phot_utils import signaltonoise, PhotometricParameters
def get_slicer_timestamps(filters, dataSlice):
    """
    Function to extract the timestamps of observations for a specific target location from an opsim DB.
    """

    print('Total number of visits to this sky location over 10 yrs:')
    time_stamps = {}
    for f in filters:
        data = []
        # Observations in the current filter:
        fdx = np.where(dataSlice['filter'] == f)[0]
        print('  ' + f + ' n_visits= ' + str(len(fdx)))
        for i in fdx:
            jd = dataSlice['observationStartMJD'][i] + 2400000.5
            data.append(jd)
        time_stamps[f] = np.array(data)

    return time_stamps

def set_photometric_parameters(exptime, nexp, readnoise=None):
    """
    Load photometric noise properties from the rubin_sim utilities
    """
    # readnoise = None will use the default (8.8 e/pixel). Readnoise should be in electrons/pixel.
    photParams = PhotometricParameters(exptime=exptime, nexp=nexp, readnoise=readnoise)
    return photParams

