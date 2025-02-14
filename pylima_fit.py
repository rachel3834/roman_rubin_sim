import numpy as np
import matplotlib.pyplot as plt
import os

from pyLIMA import event
from pyLIMA import telescopes

from pyLIMA.fits import TRF_fit
from pyLIMA.models import PSPL_model

#We need to define 
lsst_event = event.Event(ra=281.1125, dec=-13.88194)  
lsst_event.name = 'LSST_1'


datasets = [i for i in os.listdir('./output/') if '.dat' in i]
max_flux = []
for data in datasets:

    lightcurve = np.loadtxt('./output/'+data)[:,[0,3,4]]

    
    telescope_lsst = telescopes.Telescope(name=data[:-4], 
                                     camera_filter=data[-4], 
                                     lightcurve = lightcurve,
                                     lightcurve_names = ['time','mag','err_mag'],
                                     lightcurve_units = ['JD','mag','err_mag'],
                                     location = 'Earth')  
    lsst_event.telescopes.append(telescope_lsst)
    max_flux.append(lightcurve[lightcurve[:,1].argmin(),0])

pspl = PSPL_model.PSPLmodel(lsst_event, parallax=['Full', np.mean(max_flux)])
trf  = TRF_fit.TRFfit(pspl)

trf.fit_parameters['t0'][1] = [2460000,2470000]
trf.fit_parameters['tE'][1] = [0,2000]
trf.fit()


### Truth is
#event_parameters = {'t0': 2463065.8608590458, 'u0': 0.005, 'tE': 70.0, 'piEN': 0.3, 'piEE': 0.1} # Year 7
trf.fit_outputs()
plt.show()
breakpoint()

