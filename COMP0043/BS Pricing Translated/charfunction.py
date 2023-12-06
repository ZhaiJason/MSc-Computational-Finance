import numpy as np

from scipy.special import gamma
def charfunction(xi, parameters, flag=0):
    """
    Compute the characteristic function for different pricing problems.
    :param xi: Fourier space grid.
    :param parameters: Parameters for the distributions.
    :param flag: Flag for backward (0) or forward (1) characteristic function.
    :return: Characteristic function.
    """
    meancorrection = (parameters['rf'] - parameters['q']) * parameters['dt'] - \
                     np.log(charfunction0(-1j, parameters))
    F = np.exp(1j * meancorrection * xi) * charfunction0(xi, parameters)
    if flag == 0:
        F = np.conj(F)
    return F

def charfunction0(xi, parameters):
    """
    Compute the characteristic function for various distributions.
    :param xi: Fourier space grid.
    :param parameters: Parameters for the distributions.
    :return: Characteristic function.
    """
    distr_type = parameters['distr']
    
    if distr_type == 1:  # Normal
        m = parameters['m']
        s = parameters['s']
        F = np.exp(1j * xi * m - 0.5 * (s * xi)**2)

    elif distr_type == 2:  # NIG
        alpha = parameters['alpha']
        beta = parameters['beta']
        delta = parameters['delta']
        F = np.exp(-delta * (np.sqrt(alpha**2 - (beta + 1j * xi)**2) - np.sqrt(alpha**2 - beta**2)))

    elif distr_type == 3:  # VG
        theta = parameters['theta']
        s = parameters['s']
        nu = parameters['nu']
        F = (1 - 1j * xi * theta * nu + 0.5 * nu * (s * xi)**2)**(-1 / nu)

    elif distr_type == 4:  # Meixner
        alpha = parameters['alpha']
        beta = parameters['beta']
        delta = parameters['delta']
        F = (np.cos(beta / 2) / np.cosh((alpha * xi - 1j * beta) / 2))**(2 * delta)

    elif distr_type == 5:  # CGMY
        C = parameters['C']
        G = parameters['G']
        M = parameters['M']
        Y = parameters['Y']
        F = np.exp(C * gamma(-Y) * ((M - 1j * xi)**Y - M**Y + (G + 1j * xi)**Y - G**Y))

    elif distr_type == 6:  # Kou
        s = parameters['s']
        lambd = parameters['lambda']
        pigr = parameters['pigr']
        eta1 = parameters['eta1']
        eta2 = parameters['eta2']
        F = np.exp(-0.5 * (s * xi)**2 + lambd * ((1 - pigr) * eta2 / (eta2 + 1j * xi) + pigr * eta1 / (eta1 - 1j * xi) - 1))

    elif distr_type == 7:  # Merton
        s = parameters['s']
        alpha = parameters['alpha']
        lambd = parameters['lambda']
        delta = parameters['delta']
        F = np.exp(-0.5 * (s * xi)**2 + lambd * (np.exp(1j * xi * alpha - 0.5 * (delta * xi)**2) - 1))

    elif distr_type == 8:  # Levy alpha-stable
        alpha = parameters['alpha']
        beta = parameters['beta']
        gamm = parameters['gamm']
        m = parameters['m']
        c = parameters['c']
        F = np.exp(1j * xi * m - c * np.abs(gamm * xi)**alpha * (1 - 1j * beta * np.sign(xi) * np.tan(np.pi * alpha / 2)))

    else:
        raise ValueError('Invalid distribution number')

    return F

# def charfunction(xi, parameters, flag=0):
#     if len(parameters) == 2:
#         flag = 0
    
#     meancorrection = (parameters['rf'] - parameters['q']) * parameters['dt'] - np.log(charfunction0(-1j, parameters))
#     F = np.exp(1j * meancorrection * xi) * charfunction0(xi, parameters)
    
#     if flag == 0:
#         F = np.conj(F)
    
#     return F

# def charfunction0(xi, parameters):
#     distr = parameters.distr

#     if distr == 1:
#         m = parameters.m
#         s = parameters.s
#         F = np.exp(1j * xi * m - 0.5 * (s * xi)**2)

#     elif distr == 2:
#         alpha = parameters.alpha
#         beta = parameters.beta
#         delta = parameters.delta
#         F = np.exp(-delta * (np.sqrt(alpha**2 - (beta + 1j * xi)**2) - np.sqrt(alpha**2 - beta**2)))

#     elif distr == 3:
#         theta = parameters.theta
#         s = parameters.s
#         nu = parameters.nu
#         F = (1 - 1j * xi * theta * nu + 0.5 * nu * (s * xi)**2)**(-1/nu)

#     elif distr == 4:
#         alpha = parameters.alpha
#         beta = parameters.beta
#         delta = parameters.delta
#         F = (np.cos(beta / 2) / np.cosh((alpha * xi - 1j * beta) / 2))**(2 * delta)

#     elif distr == 5:
#         C = parameters.C
#         G = parameters.G
#         M = parameters.M
#         Y = parameters.Y
#         F = np.exp(C * np.math.gamma(-Y) * ((M - 1j * xi)**Y - M**Y + (G + 1j * xi)**Y - G**Y))
    
#     elif distr == 6:
#         s = parameters.s
#         lam = parameters.lam
#         pigr = parameters.pigr
#         eta1 = parameters.eta1
#         eta2 = parameters.eta2
#         F = np.exp(-0.5 * (s * xi)**2 + lam * ((1 - pigr) * eta2 / (eta2 + 1j * xi) + pigr * eta1 / (eta1 - 1j * xi) - 1))

#     elif distr == 7:
#         s = parameters.s
#         alpha = parameters.alpha
#         lam = parameters.lam
#         delta = parameters.delta
#         F = np.exp(-0.5 * (s * xi)**2 + lam * (np.exp(1j * xi * alpha - 0.5 * (delta * xi)**2) - 1))

#     elif distr == 8:
#         alpha = parameters.alpha
#         beta = parameters.beta
#         gamm = parameters.gamm
#         m = parameters.m
#         c = parameters.c
#         F = np.exp(1j * xi * m - c * np.abs(gamm * xi)**alpha * (1 - 1j * beta * np.sign(xi) * np.tan(alpha/2 * np.pi)))
    
#     else:
#         raise ValueError('Invalid distribution number')

#     return F