import numpy as np
def parameters(distr, T, dt, rf, q=None):
    param = {'distr': distr, 'T': T, 'dt': dt, 'rf': rf, 'q': 0 if q is None else q}

    if distr == 1:  # Normal
        m = 0  # mean
        s = 0.4  # standard deviation

        param.update({'m': m, 's': s * (dt ** 0.5), 'lambdam': 0, 'lambdap': 0})

    elif distr == 2:  # NIG
        alpha, beta, delta = 15, -5, 0.5
        param.update({'alpha': alpha, 'beta': beta, 'delta': delta * dt, 'lambdam': beta - alpha, 'lambdap': beta + alpha, 'FLc': delta, 'FLnu': 1})

    elif distr == 3:  # VG
        C, G, M = 4, 12, 18
        nu = 1 / C
        theta = (1 / M - 1 / G) * C
        s = ((2 * C) / (M * G)) ** 0.5

        param.update({'nu': nu / dt, 'theta': theta * dt, 's': s * (dt ** 0.5), 'lambdam': -M, 'lambdap': G})

    elif distr == 4:  # Meixner
        alpha, beta, delta = 0.3977, -1.4940, 0.3462
        param.update({'alpha': alpha, 'beta': beta, 'delta': delta * dt})

    elif distr == 5:  # CGMY
        C, G, M, Y = 4, 50, 60, 0.7
        param.update({'C': C * dt, 'G': G, 'M': M, 'Y': Y, 'lambdam': -M, 'lambdap': G, 'FLc': 2 * C * abs(np.gamma(-Y) * np.cos(np.pi * Y / 2)), 'FLnu': Y})

    elif distr == 6:  # Kou double exponential
        s, lambda_, pigr, eta1, eta2 = 0.1, 3, 0.3, 40, 12
        param.update({'s': s * (dt ** 0.5), 'lambda': lambda_ * dt, 'pigr': pigr, 'eta1': eta1, 'eta2': eta2, 'lambdam': -eta1, 'lambdap': eta2, 'FLc': (s ** 2) / 2, 'FLnu': 2})

    elif distr == 7:  # Merton jump-diffusion
        s, alpha, lambda_, delta = 0.4, 0.1, 0.5, 0.15
        param.update({'s': s * (dt ** 0.5), 'alpha': alpha, 'lambda': lambda_ * dt, 'delta': delta, 'lambdam': 0, 'lambdap': 0, 'FLc': (s ** 2) / 2, 'FLnu': 2})

    elif distr == 8:  # Stable
        alpha, beta, gamm, m, c = 2, 0, 0.3 / (2 ** 0.5), 0, 1 / ((1 + (beta * np.tan(alpha / 2 * np.pi)) ** 2) ** 0.5)
        param.update({'alpha': alpha, 'beta': beta, 'gamm': gamm * dt, 'm': m, 'c': c})

    else:
        raise Exception('Invalid distribution type')

    return param
