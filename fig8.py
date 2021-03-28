import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
from polarization import sky

plt.figure(figsize=(5,5))

theta0 = np.linspace(0, np.pi/2, 64, endpoint=False) # sun's zenith angle

s = sky(0.2) # optical depth of atmospheric layer
theta = np.r_[-s.theta, s.theta[::-1]]
N = []
for th0 in theta0:
    I = np.c_[s.transmission(th0, np.pi),
              s.transmission(th0, 0)[:,::-1]]
    P = (I[1]**2 + I[2]**2)/I[0]**2
    P_interp = interp1d(theta, P)
    N1 = [minimize_scalar(P_interp, theta[i:i+3]).x
          for i in np.flatnonzero((P[1:-1] < P[:-2])
                                  & (P[1:-1] < P[2:])
                                  & (P[1:-1] < 1e-4))]
    if len(N1)==1: N1 = [np.nan, N1[0], np.nan]
    elif len(N1)==2:
        if N1[0] < 0: N1.append(np.nan)
        else: N1.insert(0, np.nan)
    N.append(N1)

N = np.asarray(N).T/np.pi*180
plt.axis([0, 90, -90, 90])
plt.plot(theta0/np.pi*180, N[0], 'b', label=r'Arago ($\tau_0=0.2$)')
plt.plot(theta0/np.pi*180, N[1], 'g', label=r'Babinet ($\tau_0=0.2$)')
plt.plot(theta0/np.pi*180, N[2], 'r', label=r'Brewster ($\tau_0=0.2$)')

######################################################################
s = sky(0.05) # optical depth of atmospheric layer
N = []
for th0 in theta0:
    I = np.c_[s.transmission(th0, np.pi),
              s.transmission(th0, 0)[:,::-1]]
    P = (I[1]**2 + I[2]**2)/I[0]**2
    P_interp = interp1d(theta, P)
    N1 = [minimize_scalar(P_interp, theta[i:i+3]).x
          for i in np.flatnonzero((P[1:-1] < P[:-2])
                                  & (P[1:-1] < P[2:])
                                  & (P[1:-1] < 1e-4))]
    if len(N1)==1: N1 = [np.nan, N1[0], np.nan]
    elif len(N1)==2:
        if N1[0] < 0: N1.append(np.nan)
        else: N1.insert(0, np.nan)
    N.append(N1)

N = np.asarray(N).T/np.pi*180
plt.plot(theta0/np.pi*180, N[0], 'b--', label=r'Arago ($\tau_0=0.05$)')
plt.plot(theta0/np.pi*180, N[1], 'g--', label=r'Babinet ($\tau_0=0.05$)')
plt.plot(theta0/np.pi*180, N[2], 'r--', label=r'Brewster ($\tau_0=0.05$)')

plt.xlabel(r'$\theta_0$ = zenith angle of sun / deg', fontsize=14)
plt.ylabel('zenith angle of neutral point / deg', fontsize=14)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('fig8.eps')
plt.show()
