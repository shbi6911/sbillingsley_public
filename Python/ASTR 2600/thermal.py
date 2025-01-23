'''This module contains functions related to the Planck function.'''

from constants import h, k, c
import numpy as np

def planck(w_nm, T=5800):
    """
    Compute the thermal emission spectrum.

    Arguments:
     w_nm = wavelength in nm, float or numpy array
     T = temperature in K (set to 5800K by default)

    Returns:
     intensity at the given wavelength(s), in W / m**3 / sr

    Example:
     >> import thermal
     >> print thermal.planck(550, T=5800)
     2.63106560099e+13
    """

    # convert nm to m
    w = w_nm/1.0e9

    # return Planck function
    num = 2 * h *c**2
    den = w**5 * (np.exp(h * c / (w * k * T)) - 1.0)
    return num/den

def wien(T):
    """
    For a given temperature, calculate the wavelength
    (in meters) where the Planck function is maximized.
    
    Arguments:
     T = temperature in K

    Returns:
     peak wavelength, in nm

    Example:

    """

    # Wien's Law, with T in K
    return  2.898e6 / T
