import numpy as np
import matplotlib.pyplot as plt
from sky import sky
from polarization import sky as p_sky

tau0 = 0.2 # optical depth of atmospheric layer
theta0 = np.pi/3 # sun's zenith angle

s = sky(tau0, theta0)
p = p_sky(tau0)

plt.subplot2grid((1,3),(0,0),1,2)
plt.axis([-90, 90, -0.05, 0.35])
plt.plot(np.r_[-p.theta, p.theta[::-1]]/np.pi*180,
         np.c_[p.transmission(theta0, np.pi),
               p.transmission(theta0, 0)[:,::-1]].T)
plt.plot(np.r_[-s.theta, s.theta[::-1]]/np.pi*180,
         np.r_[s.transmission(np.pi, tau0),
               s.transmission(0, tau0)[::-1]], 'c--')
plt.legend(['I','Q','U','I (no polarization)'], fontsize=14)
plt.xlabel(r'$\theta$ = zenith angle / deg', fontsize=14)
plt.ylabel(r'$I$ = stokes parameters / $F_{\odot}/\pi$', fontsize=14)
plt.text(50, -0.04, r'$\phi = \phi_0$', fontsize=14)

plt.subplot2grid((1,3),(0,2))
plt.axis([0, 90, -0.05, 0.35])
plt.plot(p.theta/np.pi*180,
         p.transmission(theta0, np.pi/2).T)
plt.plot(s.theta/np.pi*180,
         s.transmission(np.pi/2, tau0), 'c--')
plt.tick_params('y', labelleft=False)
plt.xlabel(r'$\theta$', fontsize=14)
plt.text(3, -0.04, r'$\phi = \phi_0 + 90^{\circ}$', fontsize=14)
plt.text(3, 0.325, r'$\tau_0 = 0.2$', fontsize=14)
plt.text(3, 0.3, r'$\theta_0 = 60^{\circ}$', fontsize=14)

plt.tight_layout()
plt.savefig('fig1.eps')
plt.show()

plt.close()
# fig2 ###########################################################

tau0 = 0.2 # optical depth of atmospheric layer
theta0 = 80/180*np.pi # sun's zenith angle

s = sky(tau0, theta0)
p = p_sky(tau0)

plt.subplot2grid((1,3),(0,0),1,2)
plt.axis([-90, 90, -0.05, 0.35])
plt.plot(np.r_[-p.theta, p.theta[::-1]]/np.pi*180,
         np.c_[p.transmission(theta0, np.pi),
               p.transmission(theta0, 0)[:,::-1]].T)
plt.plot(np.r_[-s.theta, s.theta[::-1]]/np.pi*180,
         np.r_[s.transmission(np.pi, tau0),
               s.transmission(0, tau0)[::-1]], 'c--')
plt.legend(['I','Q','U','I (no polarization'], fontsize=14)
plt.xlabel(r'$\theta$ = zenith angle / deg', fontsize=14)
plt.ylabel(r'$I$ = stokes parameters / $F_{\odot}/\pi$', fontsize=14)
plt.text(50, -0.04, r'$\phi = \phi_0$', fontsize=14)

plt.subplot2grid((1,3),(0,2))
plt.axis([0, 90, -0.05, 0.35])
plt.plot(p.theta/np.pi*180,
         p.transmission(theta0, np.pi/2).T)
plt.plot(s.theta/np.pi*180,
         s.transmission(np.pi/2, tau0), 'c--')
plt.tick_params('y', labelleft=False)
plt.xlabel(r'$\theta$', fontsize=14)
plt.text(3, -0.04, r'$\phi = \phi_0 + 90^{\circ}$', fontsize=14)
plt.text(3, 0.325, r'$\tau_0 = 0.2$', fontsize=14)
plt.text(3, 0.3, r'$\theta_0 = 80^{\circ}$', fontsize=14)

plt.tight_layout()
plt.savefig('fig2.eps')
plt.show()
