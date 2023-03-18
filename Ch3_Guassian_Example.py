"""
UCSD PHYS/SIO 111:
Ch3 Example of Gaussian Packet IC with d\eta/dt = 0

transaltion:Terrindeep Sandhu WI23
"""

import numpy as np
import matplotlib.pyplot as plt

# define constants and parameters
grav = 9.81
Lx = 2e5
nx = int(2e6)
dx = Lx / nx
xx = np.arange(0, Lx - dx, dx) - Lx / 2
k1 = 2 * np.pi / Lx
kk = np.concatenate((np.arange(0, nx / 2), np.arange(-nx / 2 + 1, 0))) * k1
gaussian_width = 20

# deep water dispersion relationship
omega = np.sqrt(grav * np.abs(kk))

# define initial condition in physical space
Fx = (gaussian_width * np.sqrt(2 * np.pi)) ** (-1) * np.exp(-xx ** 2 / (2 * gaussian_width ** 2)) * np.cos(2 * xx)

# define initial condition in Fourier space
Ak = np.fft.fft(Fx)

# plot initial condition
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
ax1.plot(xx, Fx)
ax1.set_xlabel('x (m)')
ax1.set_ylabel(r'$\eta$ (m)')
ax1.set_title('Initial Condition of a Gaussian hump, d$\eta$/dt = 0')
ax1.set_xlim(-5 * gaussian_width, 5 * gaussian_width)
ax1.set_ylim(-1*max(Fx), max(Fx))

# evolve the solution and plot it
for j in range(1, 401, 10):
    t1 = (j - 1) * 10
    Feta = Ak * np.exp(1j * omega * t1) + np.conj(Ak) * np.exp(-1j * omega * t1)
    etaa = np.fft.ifft(Feta)
    eta = np.real(etaa)
    
    ax2.cla()
    ax2.plot(xx, eta)
    ax2.set_xlabel('x (m)')
    ax2.set_ylabel(r'$\eta$ (m)')
    ax2.set_title(f'Gaussian bump IC in deep water: t = {t1:04d} (s)')
    ax2.set_xlim(0, Lx / 20)
    ax2.set_ylim(min(Fx), max(Fx))

    ax3.cla()
    ax3.plot(xx, eta)
    ax3.set_xlabel('x (m)')
    ax3.set_ylabel(r'$\eta$ (m)')
    ax3.set_title(f'Blowup of middle panel: Gaussian bump IC in deep water: t = {t1:04d} (s)')
    ax3.set_xlim(0, Lx / 60)
    ax3.set_ylim(min(Fx), max(Fx))
    
    plt.pause(0.25)

plt.show()
