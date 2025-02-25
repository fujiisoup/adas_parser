import numpy as np
from scipy import integrate
import xarray as xr


def calculate_rate(
        energy, crosssection, temperature, mass,
        interp_method='trapz'):
    '''Calculate the rate coefficients from the crosssection, 
    assuming the Maxwell energy distribution.

    Parameters
    ----------
    energy: array-like
        Energy points in eV.
    crosssection: array-like
        Crosssection values in m^2.
    temperature: float
        Temperature in eV.
    mass: float
        Mass of the particles in amu.
    interp_method: str
        Interpolation method for the crosssection. Default is 'trapz'.
    kwargs: dict
        Additional keyword arguments for the interpolation function.

    Returns
    -------
    rate: array-like
        Rate coefficients in m^3/s.
    '''
    if isinstance(crosssection, xr.DataArray):
        crosssection = crosssection.transpose(
            *([d for d in crosssection.dims if d != 'energy'] + ['energy'])
        )
    
    # Convert temperature to joule
    kT = np.array(temperature) * 1.602176565e-19
    if kT.ndim > 0:
        kT = kT[(slice(None),) + (np.newaxis, ) * crosssection.ndim]
    # Convert mass to kg
    mass = mass * 1.66053906660e-27
    # Convert energy to Joules
    energy = np.array(energy) * 1.602176634e-19

    # integrand
    maxwell = 2 * np.sqrt(energy / (np.pi * kT)**3) * np.exp(-energy / kT)
    velocity = np.sqrt(2 * energy / mass)
    integrand = np.array(crosssection) * velocity * maxwell

    rate = integrate.trapezoid(y=integrand, x=energy, axis=-1)
    if isinstance(crosssection, xr.DataArray):
        dims = [d for d in crosssection.dims if d != 'energy']
        if np.array(temperature).ndim == 1:
            dims = ['temperature'] + dims
        rate = xr.DataArray(
            rate, dims=dims,
            coords={'temperature': temperature}
        )
        rate['temperature'].attrs['units'] = 'eV'
        for k, c in crosssection.coords.items():
            if k != "energy":
                rate.coords[k] = c
    return rate

