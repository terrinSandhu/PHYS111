"""
UCSD PHYS/SIO 111:
Ch4 Distant Storm

transaltion:Terrindeep Sandhu WI23
"""

import numpy as np
import matplotlib.pyplot as plt

grav = 9.81
Lx = 4e+6
nx = 14e+6
dx = Lx / nx
xx = np.arange(0, Lx - dx, dx) - Lx / 2

k1 = 2 * np.pi / Lx
kk = k1 * np.concatenate((np.arange(0, nx / 2), np.arange(-nx / 2 + 1, 0)))

gaussian_width = 2

depth = 4000
omega = np.sqrt(grav * np.abs(kk))

Fx = (gaussian_width * np.sqrt(2 * np.pi)) ** (-1) * np.exp(-xx ** 2 / (2 * gaussian_width ** 2))
Ak = np.fft.fft(Fx)

time1 = 4 * 3600
time2 = 6 * 3600
time = time2
Feta = Ak * np.exp(1j * (omega * time))
eta = 2 * np.real(np.fft.ifft(Feta))

fig = plt.figure(figsize=(10, 15))
ax1 = fig.add_subplot(3, 1, 1)
ax1.plot(xx, Fx)
ax1.set_xlabel('x (m)')
ax1.set_ylabel('$\eta$ (m)')
ax1.set_title('Initial Condition of a $\eta$ = Gaussian hump, $d\eta/dt=0:$')
ax1.set_xlim([-5 * gaussian_width, 5 * gaussian_width])

xxkm = xx / 1000
Lxkm = Lx / 1000

ax2 = fig.add_subplot(3, 1, 2)
ax2.plot(xxkm, eta)
ax2.set_xlabel('x (km)')
ax2.set_ylabel('$\eta$ (m)')
ax2.set_title(f'Distant storm Case: t={time:04d} (s)')
ax2.set_xlim([0, Lxkm / 2])
ax2.grid()

ax3 = fig.add_subplot(3, 1, 3)
ax3.plot(xxkm, eta)
ax3.set_xlabel('x (km)')
ax3.set_ylabel('$\eta$ (m)')
ax3.set_title(f'Blowup of middle panel: t={time:04d} (s)')
ax3.set_xlim([0, Lxkm / 10])
ax3.grid()

plt.show()

x1 = 100 * 1000
cg1 = x1 / time1
k1 = 9.81 / (4 * cg1 ** 2)
lambda1 = 2 * np.pi / k1
lambda1A = 2000 / 16
print([lambda1, lambda1A])

A1_t1 = 7.9e-4

x2 = 200 * 1000
cg2 = x2 / time1
k2 = 9.81 / (4 * cg2 ** 2)
lambda2 = 2 * np.pi / k2
lambda2A = 4000 / 8
print([lambda2, lambda2A])

A2_t1 = 3e-4

x1_t2 = x1 + cg1 * (time2 - time1)
ax2.set_xlim([x1_t2 - 1000, x1_t2 + 1000] / 1000)
A1_t2 = 6.6e-4
A1_t2A = A1_t1 * np.sqrt(time1 / time2)

print("A1_t2 =", np.concatenate((A1_t1, A1_t2A)))

x2_t2 = x1+ cg2*(time2-time1)

L80 = 310; k80 = 2*np.pi/L80  
omega80 = np.sqrt(grav*k80) 
omega80 = np.sqrt(grav*k80*np.tanh(k80*depth)) 
cg80 = 0.5*np.sqrt(grav/k80)

L50 = 130; k50 = 2*np.pi/L50  
omega50 = np.sqrt(grav*k50)   
cg50 = 0.5*np.sqrt(grav/k50)
