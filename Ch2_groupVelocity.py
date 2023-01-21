"""
UCSD PHYS/SIO 111:
Chapter 2 grooup velocity visualization

transaltion:Terrindeep Sandhu WI23
"""
import time
import numpy as np
import matplotlib.pyplot as plt


grav = 9.81

x = (np.arange(0, 60, 0.1)).tolist()

lambda1=3;     # first wavelength
lambda2 = 3.3  # second wavelength
k1 = 2*np.pi/lambda1
k2 = 2*np.pi/lambda2

omega1 = np.sqrt(grav*k1)
omega2 = np.sqrt(grav*k2)

t = np.arange(0, 40, 0.1)
nt = len(t)


plt.ion()
fig, ax = plt.subplots() 


for i in range(1,nt):
  ax.cla()

  k1N = [( (k1 * j)- omega1*t[i]) for j in x]
  eta1 = np.cos(k1N )
  k2N = [( (k2 * j)- omega2*t[i]) for j in x]
  eta2 = np.cos(k2N)

  eta =  [sum(z) for z in zip(eta1, eta2 )]  

  ax.plot(x,eta)
  ax.set_xlabel('x (m)')
  ax.set_ylabel('\eta (m)')
  tStr = 'wave packet example: lambda1=' + "{:.1f}".format(lambda1)  + '(m), lambda2=' + "{:.1f}".format(lambda2) +' (m): t=' + "{:.1f}".format(t[i]) +'(s)'
  ax.set_title(tStr)
  
  
  fig.canvas.draw()
  fig.canvas.flush_events()
  time.sleep(0.1)

