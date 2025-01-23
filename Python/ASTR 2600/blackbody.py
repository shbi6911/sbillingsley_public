#not written by Shane Billingsley

'''
This module contains tools to generate random blackbody spectra.
'''

import matplotlib.pyplot as plt
import numpy as np
import thermal as thermal

def generateBlackbody(T='random'):
    '''
    Generate a blackbody spectrum for a random (or specific) temperature.

        Arguments:

            T = [Optional]If set, this temperature is used to generate the
		blackbody curve. If not, a random temperature is generated
		between 2000-9000K.

        Returns:

            A blackbody spectrum. The format of this spectrum
            is as a dictionary that contains several key-value pairs.
            The keys inside this spectrum are highlighted in
            the example below.

        Example:

            # call the function to get a random spectrum
            s = generateBlackbody()

            # the wavelengths are in units of nanometers (=10**-9 m)
            w_nanometers = s['wavelengths']

            # the fluxes are in units of "W/m**3"
            f = s['flux']

            # The temperature used for the Planck function is in Kelvin"
            t = s['temp']

            # A title for the object is given in a string named objectID
            name = s['objectID']


    '''
    try:
        temperature = float(T)
    except (TypeError, ValueError):
        temperature = np.round(np.random.uniform(2000,9000))

    w = np.sort(np.random.uniform(300,1000, 1000))
    f = thermal.planck(w, T=temperature)

    # populate a fake spectrum
    s = {}

    # store the seed temperature
    s['temp'] = temperature

    # store wavelengths
    s['wavelengths'] = w

    # store fluxes
    s['flux'] = f

    # store the objectID for this fake spectrum
    s['objectID'] = 'A {:4.0f}K blackbody'.format(temperature)

    return s

def testBlackbody():
    '''
    This function tests the generateBlackbody() function and plots the result.
        Arguments:
            None

        Returns:
            Nothing

        Example:
            testBlackbody()

    '''
    s = generateBlackbody()

    # plot a spectrum
    plt.plot(s["wavelengths"], s["flux"])
    plt.title(s['objectID'])
    plt.xlabel("Wavelength (nanometers)")
    plt.ylabel("Flux (W/m$^3$)")
    plt.ylim(0, None)
    plt.show()
