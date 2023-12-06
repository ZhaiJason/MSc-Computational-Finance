import numpy as np
from numpy.fft import fft, fftshift, ifftshift
from charfunction import charfunction

def kernel(ngrid, xmin, xmax, parameters, alpha=0, disc=1, flag=0):
    """
    Set up grids, compute the characteristic function, and apply discounting.

    Parameters:
    ngrid (int): Number of grid points.
    xmin (float): Minimum value of the grid in real space.
    xmax (float): Maximum value of the grid in real space.
    parameters (dict): Dictionary of parameters for the characteristic function.
    alpha (float): Shift parameter, especially for Feng-Linetsky and convolution.
    disc (int): Discount factor in the density (0 for no, 1 for yes).
    flag (int): Type of characteristic function (0 for backward, 1 for forward).

    Returns:
    tuple: Grid in real space (x), density (h), grid in Fourier space (xi), Fourier transform of density (H).
    """

    N = ngrid // 2
    dx = (xmax - xmin) / ngrid
    x = dx * np.arange(-N, N)
    dxi = 2 * np.pi / (xmax - xmin)
    xi = dxi * np.arange(-N, N)

    # Assuming charfunction is defined elsewhere
    H = charfunction(xi + 1j * alpha, parameters, flag)  # characteristic function

    if disc == 1:
        H *= np.exp(-parameters.rf * parameters.dt)  # apply discount

    h = np.real(fftshift(fft(ifftshift(H)))) / (xmax - xmin)  # discounted kernel

    return x, h, xi, H