from pyLIMA.simulations import simulator
from pyLIMA import event, toolbox
from pyLIMA import telescopes
from pyLIMA.models import PSPL_model
from pyLIMA.outputs import pyLIMA_plots
import pyLIMA.fits.objective_functions
from pyLIMA.toolbox import fake_telescopes, plots

def generate_sim_event(ra, dec, photometry, event_parameters):
    sim_event = event.Event(ra=ra, dec=dec)
    sim_event.name = 'Rubin microlensing event'

    # Create lightcurves in each Rubin passband, using the timestamps of observations in those passbands
    LSST_u = telescopes.Telescope(name='u', camera_filter='u', location='Earth',
                                  light_curve=photometry['u'].astype(float),
                                  light_curve_names=['time', 'mag', 'err_mag'], light_curve_units=['JD', 'mag', 'mag'])

    LSST_g = telescopes.Telescope(name='g', camera_filter='g', location='Earth',
                                  light_curve=photometry['g'].astype(float),
                                  light_curve_names=['time', 'mag', 'err_mag'], light_curve_units=['JD', 'mag', 'mag'])

    LSST_r = telescopes.Telescope(name='r', camera_filter='r', location='Earth',
                                  light_curve=photometry['r'].astype(float),
                                  light_curve_names=['time', 'mag', 'err_mag'], light_curve_units=['JD', 'mag', 'mag'])

    LSST_i = telescopes.Telescope(name='i', camera_filter='i', location='Earth',
                                  light_curve=photometry['i'].astype(float),
                                  light_curve_names=['time', 'mag', 'err_mag'], light_curve_units=['JD', 'mag', 'mag'])

    LSST_z = telescopes.Telescope(name='z', camera_filter='z', location='Earth',
                                  light_curve=photometry['z'].astype(float),
                                  light_curve_names=['time', 'mag', 'err_mag'], light_curve_units=['JD', 'mag', 'mag'])

    LSST_y = telescopes.Telescope(name='y', camera_filter='y', location='Earth',
                                  light_curve=photometry['y'].astype(float),
                                  light_curve_names=['time', 'mag', 'err_mag'], light_curve_units=['JD', 'mag', 'mag'])

    sim_event.telescopes.append(LSST_u)
    sim_event.telescopes.append(LSST_g)
    sim_event.telescopes.append(LSST_r)
    sim_event.telescopes.append(LSST_i)
    sim_event.telescopes.append(LSST_z)
    sim_event.telescopes.append(LSST_y)

    # Create a PyLIMA event model with the given event event_parameters
    event_model = PSPL_model.PSPLmodel(sim_event, parallax=['Full', event_parameters['t0']])

    return sim_event, event_model