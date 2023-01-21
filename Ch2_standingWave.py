"""
UCSD PHYS/SIO 111:
Shallow water standing wave example

transaltion:Terrindeep Sandhu WI23
"""
import numpy as np
import matplotlib.pyplot as plt
import time

grav = 9.81
h=1
a = 0.05  # individual wave amplitude
L=60  # length of basiin

x = np.arange(0, L, 0.1)

mode_number = 2

k = np.pi*mode_number/L

omega = np.sqrt(grav*h)*k

u0 = a*omega/(k*h)  #shallow water velocity amplitude

t = np.arange(0, 80, 0.2)
nt = len(t)


plt.ion()
fig,axs = plt.subplots(2)

for i in range(1,nt):

    axs[0].cla()
    axs[1].cla()

    eta1 = a*np.cos(k*x - omega*t[i])   # right going wave
    eta2 = a*np.cos(k*x + omega*t[i])     # left going wave
    eta = eta1+eta2                   #add them together
    u1 = u0*np.sin(k*x - omega*t[i])   # same for velocity
    u2 = u0*np.sin(k*x + omega*t[i])
    u = u1+u2

    #subplot(2,1,1)  
    axs[0].plot(x,eta)
    axs[0].set_xlabel('x (m)')
    axs[0].set_ylabel('\eta (m)')
    axs[0].set_title('standing wave example: sea surface L='+ "{:.1f}".format(L) +' (m), mode=' +str(mode_number) + ': t='+ "{:.1f}".format(t[i])+ '(s)')
    axs[0].axis([0, L, -2*a, 2*a])

    #subplot(2,1,2)  
    axs[1].plot(x,u)
    axs[1].set_xlabel('x (m)')
    axs[1].set_ylabel('u  (m/s)')
    axs[1].set_title('velocity ')
    axs[1].axis([0, L, -2*u0, 2*u0])

    fig.canvas.draw()
    fig.canvas.flush_events()
    time.sleep(0.05)
