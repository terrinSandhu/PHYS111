"""
UCSD PHYS/SIO 111:
Chapter 1 Demo for wave and wave velocity visualization

transaltion:Terrindeep Sandhu WI23
"""
import numpy as np
import time
import matplotlib.pyplot as plt
import get_waveNum

a = 0.5  #wave amplitude, m
f = 0.05  #wave freq, Hz
h = 100  #water depth, m
grav = 9.81  # gravity  m/s^2

omega = 2* np.pi* f  # radian wave frequencey  rad/s

k =  get_waveNum.getWaveNum(omega,h)  #wavenumber -- getWaveNum() ???


ak = a*k # nond parameters
kh = k*h # non-d parameter

Lambda = 2*np.pi/k  # wavelength

Lx = 600  # domain size, m
dx = 0.5  # grid spacing, m
x = np.arange(0, Lx, dx)

t=0  # time, s
dt = 0.5 # time step, s

tEnd = 100 # end time to stop

xmid = Lx/2
zA = -10  # first current meter vertical location, meters
zB = -70 # 2nd current meter vertical location, meters

if (abs(zA)>h or abs(zB)>h):
    print('ERROR: current meter below the bottom of the ocean')
    


plt.ion()
fig,axs = plt.subplots(2)



while t < tEnd:
    axs[0].cla()
    axs[1].cla()

    eta = a*np.cos( k*x - omega*t)

  
    axs[0].plot(x,eta)
    titleStr = "t=" + str(t)+ "(s): a=" +"{:.1f}".format(a) + "(m), f=" + "{:.1f}".format(f) +"(Hz), h=" + "{:.1f}".format(h) + "(m): ak=" + "{:.1f}".format(ak) +",  kh="+ "{:.1f}".format(kh)
    axs[0].set_title(titleStr)



    axs[0].plot( [xmid, xmid],[-2*a, 2*a],'r--' )
  
    axs[0].axis([0, Lx, -2*a,  2*a ])
    axs[0].set_xlabel('x (m)')
    axs[0].set_ylabel('z (m)')
    

    uA = a*omega*( np.cosh(k*(zA+h))/np.sinh(k*h) )*np.cos(k*xmid-omega*t)
    wA = a*omega*( np.sinh(k*(zA+h))/np.sinh(k*h) )*np.sin(k*xmid-omega*t)

    uB = a*omega*( np.cosh(k*(zB+h))/np.sinh(k*h) )*np.cos(k*xmid-omega*t)
    wB = a*omega*( np.sinh(k*(zB+h))/np.sinh(k*h) )*np.sin(k*xmid-omega*t)

    axs[1].quiver([xmid, xmid],[zA, zB],[uA, uB],[wA, wB],0.5)

    axs[1].plot([ xmid+0.75*zB, xmid-0.75*zB],[-h, -h],'k--')

    axs[1].set_xlabel('x (m)')
    axs[1].set_ylabel('z (m)')

    axs[1].axis([xmid+0.75*zB, xmid-0.75*zB,  -1.2*h, -1.5*zA ])



    

    fig.canvas.draw()
    fig.canvas.flush_events()
    time.sleep(0.1)


    t = t+dt
    #print(t)
