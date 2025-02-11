import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
import rubin_sim.maf as maf
import rubin_sim.phot_utils as phot_utils
from rubin_sim.data import get_data_dir
from rubin_sim.data import get_baseline
from rubin_sim.phot_utils import signaltonoise, PhotometricParameters
from pyLIMA.simulations import simulator
from pyLIMA import event, toolbox
from pyLIMA import telescopes
from pyLIMA.models import PSPL_model
from pyLIMA.outputs import pyLIMA_plots
import pyLIMA.fits.objective_functions
from pyLIMA.toolbox import fake_telescopes, plots
import opsim_utils
import photometry_utils
import pylima_utils

##############
## Source star parameters
# The following parameters describe the source star for the simulated event

# Target coordinates
target = SkyCoord(281.1125, -13.88194, frame='icrs', unit=(u.deg, u.deg))

ra = target.ra.deg
dec = target.dec.deg

print('Target sky location: RA=' + str(ra) + 'deg, Dec=' + str(dec))

# Target distance in parsecs, modulus in mag
distance = 22.4e6
distance_modulus = 31.75
mag_offset = 0.0

# Extinction in that direction [mags]
Av = 0.29

# Source star unmagnified brightness in all filters
target_mag = {'u': 24.31, 'g': 21.81, 'r': 22.81, 'i': 23.31, 'z': 24.31, 'y': 25.11}
print('Source star unlensed magnitudes:' + repr(target_mag))

# Photometric zeropoint
ZP = pyLIMA.toolbox.brightness_transformation.ZERO_POINT
print('Photometric zeropoint in use = ' + str(ZP) + ' mag')

##############
## Event parameters
# Here we define the parameters of the microlensing event

# Microlensing event parameters.  Event t0 is set in July 16, 2026
#event_parameters = {'t0': 2461238.325405, 'u0': 0.005, 'tE': 70.0, 'piEN': 0.3, 'piEE': 0.1}

# Alternative parameter set where the event peaks in year 7, when the number of
# visits to the Roman field is augmented
event_parameters = {'t0': 2463065.8608590458, 'u0': 0.005, 'tE': 70.0, 'piEN': 0.3, 'piEE': 0.1} # Year 7
print('Simulated event parameters: ' + repr(event_parameters))

##############
## Rubin OpSim
# Load the LSST filterset and expected throughput in each filter
LSST_BandPass = {}
lsst_filterlist = 'ugrizy'
for f in lsst_filterlist:
    LSST_BandPass[f] = phot_utils.Bandpass()
    # print(os.path.join(get_data_dir(), 'throughputs', 'baseline', f'total_{f}.dat'))
    LSST_BandPass[f].read_throughput(os.path.join(get_data_dir(), 'throughputs', 'baseline', f'total_{f}.dat'))
photParams = opsim_utils.set_photometric_parameters(30,1)

# Rubin limiting magnitudes, typical for crowded regions
mag_limit = [23.7, 24.97, 24.52, 24.13, 23.56, 22.55]

# In order to simulate microlensing event lightcurves with time sampling realistic
# to that of the Rubin LSST, use load information from an LSST Operational Simuation
# database
print('Loading LSST OpSim database...')

# Load the opsim database
baseline_file = '/Users/rstreet/LSST/SCOC/baseline_v4.0_10yrs.db'
name = os.path.basename(baseline_file).replace('.db','')

outDir = 'temp'
resultsDb = maf.db.ResultsDb()

# This section selects which columns from the database will be loaded - in this case, the information on the
# bandpass and MJD of each exposure, together with the 5sigma limiting magnitude.
metric = maf.metrics.PassMetric(cols=['filter', 'observationStartMJD', 'fiveSigmaDepth'])
sql = ''

# Rubin filterset, and the single-visit limiting magnitude in each of those filters
filters = ['u', 'g', 'r', 'i', 'z', 'y']
mag_limit = {'u': 23.7, 'g': 24.97, 'r': 24.52, 'i': 24.13, 'z': 23.56, 'y': 22.55}

# Create a slicer for the survey HEALpixel corresponding to this target
slicer = maf.slicers.UserPointsSlicer(ra=[ra], dec=[dec])
bundleList = [maf.MetricBundle(metric, slicer, sql)]

# Run the chosen metrics for this slicer to get the opsim results for this target
bg = maf.MetricBundleGroup(
    bundleList, baseline_file, out_dir=outDir)
bg.run_all()
dataSlice = bundleList[0].metric_values[0]

# For the purposes of simulation we extract the timestamps of the sequence of exposures into an array
# for ease of later handling.
time_stamps = opsim_utils.get_slicer_timestamps(filters, dataSlice)


############
## Event lightcurve simulation

# Generate an unlensed lightcurve for the source star with realistic time-sampling
# and photometric noise
photometry, baseline_flux, flux_parameters = photometry_utils.generate_event_photometry_arrays(time_stamps, target_mag, mag_offset, ZP, dataSlice, LSST_BandPass, photParams)

# Inject the microlensing event into the simulated source photometry
sim_event, event_model = pylima_utils.generate_sim_event(ra, dec, photometry, event_parameters)

# Create a list of the parameters to be simulated for this type of event, then append the source flux parameters
fs = []
flux_params = []

for f in filters:
    fs.append(flux_parameters[f]['f_source'])
    # The flux parameters PyLIMA seems to expect by default (per telescope) are [
    flux_params.append(flux_parameters[f]['f_source'])
    flux_params.append(flux_parameters[f]['f_source'] + flux_parameters[f]['f_blend'])

param_order = ['t0', 'u0', 'tE', 'piEN', 'piEE']
param_set = []
for key in param_order:
    param_set.append(event_parameters[key])
param_set += flux_params

# Simulate the event lightcurve
pyL_params = event_model.compute_pyLIMA_parameters(param_set)
simulator.simulate_lightcurve_flux(event_model, pyL_params)

# Update the photometric uncertainties to reflect the changed fluxes
event_model = photometry_utils.update_phot_uncertainties(event_model, dataSlice, filters, LSST_BandPass, photParams, ZP)


################
## Output simulated lightcurves
CWD = os.getcwd()

for k, tel in enumerate(event_model.event.telescopes):

    with open(os.path.join(CWD, 'output', 'lsst_event_' + tel.name + '.dat'), 'w') as f:
        f.write('# JD    flux    err_flux    mag   err_mag\n')
        for i in range(0, len(tel.lightcurve_flux),1):
            f.write(
                str(tel.lightcurve_flux['time'][i].value) + ' '
                + str(tel.lightcurve_flux['flux'][i]) + ' '
                + str(tel.lightcurve_flux['err_flux'][i]) + ' '
                + str(tel.lightcurve_magnitude['mag'][i].value) + ' '
                + str(tel.lightcurve_magnitude['err_mag'][i]) + '\n'
            )
