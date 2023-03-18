"""
UCSD PHYS/SIO 111:
Ch3 Example of Gaussian Packet IC with d\eta/dt = 0

transaltion:Terrindeep Sandhu WI23
"""

import numpy as np
import matplotlib.pyplot as plt

grav = 9.81
Lx = 2e5        # total domain size here 40 km
nx = 2e6        # total number of grid points
dx = Lx / nx    # grid
xx = np.arange(0, Lx-dx, dx) - Lx/2   # x domain
k1 = 2 * np.pi / Lx
kk = k1 * np.concatenate((np.arange(0, nx/2), np.arange(-nx/2+1, 0)))   # wavenumbers

h = 10   # try h=1000, h=100, h=10 m depth
gaussian_width = 50

# dispersion relationship
omega = np.sqrt(grav * np.abs(kk) * np.tanh(np.abs(kk) * h))

# IC in physical space
Fx = (gaussian_width * np.sqrt(2 * np.pi))**(-1) * np.exp(-xx**2 / (2 * gaussian_width**2))

# IC in Fourier space because d\eta/dt = 0
Ak = np.fft.fft(Fx)

# make the eta IC plot
fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10, 15))
axs[0].plot(xx, Fx)
axs[0].set_xlabel('x (m)')
axs[0].set_ylabel(r'$\eta$ (m)')
axs[0].set_title(f'Initial Condition of a Gaussian hump, d$\eta$/dt=0: h={h:.2f} (m)')
axs[0].set_xlim([-5*gaussian_width, 5*gaussian_width])
ylims = [-0.5 * np.max(Fx), np.max(Fx)]

print('hit enter to start')
input()
for j in range(100):
    t1 = (j - 1) * 10
    Feta = Ak * np.exp(1j * omega * t1) + Ak * np.exp(-1j * omega * t1)
    etaa = np.fft.ifft(Feta)
    eta = np.real(etaa)

    axs[1].cla()
    axs[1].plot(xx, eta), axs[1].grid()
    axs[1].set_xlabel('x (m)')
    axs[1].set_ylabel(r'$\eta$ (m)')
    axs[1].set_title(f'Gaussian bump IC in deep water: t={t1:04d} (s)')
    axs[1].set_xlim([-Lx/10, Lx/10])
    axs[1].set_ylim(ylims)

    axs[2].cla()
    axs[2].plot(xx, eta), axs[2].grid()
    axs[2].set_xlabel('x (m)')
    axs[2].set_ylabel(r'$\eta$ (m)')
    axs[2].set_title(f'Blowup of middle panel: Gaussian bump IC in deep water: t={t1:04d} (s)')
    axs[2].set_xlim([-Lx/30, Lx/30])
    axs[2].set_ylim(ylims)

    plt.pause(0.25)

plt.show()
