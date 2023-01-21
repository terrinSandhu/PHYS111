# Falk Feddersen (c) 2001
#(transaltion: Terrindeep Sandhu)
#
# function that takes the radian wave frequency and
# a vector of depths and returns the wavenumber at that
# depth by solving the dispersion relationship
#
# function k = get_wavenum(omega,h)

import numpy as np
import mpmath as mpm
def getWaveNum(omega,h):

# returns the wavenumber of the gravity wave
# dispersion relation, by using newtons method

# the initial guess will be the shallow water wavenumber

    g = 9.81

    k = omega/np.sqrt(g*h)


    f = g*k*np.tanh((k*h)) - omega**2

    while abs(f)>1e-10:
        dfdk = g*k*h*(mpm.sech(k*h))**2 + g*mpm.tanh(k*h)
        k = k - f/dfdk
        f = g*k*mpm.tanh(k*h) - omega**2
        print(k)
    return float(k)
