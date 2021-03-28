import numpy as np
import matplotlib.pyplot as plt
from polarization import sky

theta0 = np.pi/3 # sun's zenith angle

albedo = [0, 0.2] # reflectivity of ground
label = [r'$\tau_0=0.2$, $a=0$', r'$\tau_0=0.2$, $a=0.2$',
         r'$\tau_0=0.05$, $a=0$', r'$\tau_0=0.05$, $a=0.2$']

plt.axis([-90, 90, 0, 1])

s = sky(0.2) # optical depth of atmospheric layer
theta = np.r_[-s.theta, s.theta[::-1]]
for a in albedo:
    I = np.c_[s.transmission(theta0, np.pi, a),
              s.transmission(theta0, 0, a)[:,::-1]]
    P = np.sqrt(I[1]**2 + I[2]**2)/I[0] # degree of polarization
    plt.plot(theta/np.pi*180, P)

s = sky(0.05) # optical depth of atmospheric layer
for a in albedo:
    I = np.c_[s.transmission(theta0, np.pi, a),
              s.transmission(theta0, 0, a)[:,::-1]]
    P = np.sqrt(I[1]**2 + I[2]**2)/I[0] # degree of polarization
    plt.plot(theta/np.pi*180, P)

plt.xlabel(r'$\theta$ = zenith angle / deg', fontsize=14)
plt.ylabel(r'$P$ = degree of polarization', fontsize=14)
plt.text(-85, 0.93, r'$\theta_0 = 60^{\circ}$', fontsize=14)
plt.legend(label, fontsize=12)
plt.tight_layout()
plt.savefig('fig6.eps')
plt.show()

plt.close()
# fig7 ##############################################################

theta0 = 80/180*np.pi # sun's zenith angle

albedo = [0, 0.2] # reflectivity of ground
label = [r'$\tau_0=0.2$, $a=0$', r'$\tau_0=0.2$, $a=0.2$',
         r'$\tau_0=0.05$, $a=0$', r'$\tau_0=0.05$, $a=0.2$']

plt.axis([-90, 90, 0, 1])

s = sky(0.2) # optical depth of atmospheric layer
theta = np.r_[-s.theta, s.theta[::-1]]
for a in albedo:
    I = np.c_[s.transmission(theta0, np.pi, a),
              s.transmission(theta0, 0, a)[:,::-1]]
    P = np.sqrt(I[1]**2 + I[2]**2)/I[0] # degree of polarization
    plt.plot(theta/np.pi*180, P)

s = sky(0.05) # optical depth of atmospheric layer
for a in albedo:
    I = np.c_[s.transmission(theta0, np.pi, a),
              s.transmission(theta0, 0, a)[:,::-1]]
    P = np.sqrt(I[1]**2 + I[2]**2)/I[0] # degree of polarization
    plt.plot(theta/np.pi*180, P)

plt.xlabel(r'$\theta$ = zenith angle / deg', fontsize=14)
plt.ylabel(r'$P$ = degree of polarization', fontsize=14)
plt.text(-85, 0.93, r'$\theta_0 = 80^{\circ}$', fontsize=14)
plt.legend(label, fontsize=12)
plt.tight_layout()
plt.savefig('fig7.eps')
plt.show()
