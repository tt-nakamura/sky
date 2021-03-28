import numpy as np
import matplotlib.pyplot as plt
from polarization import sky

# fig3.1 ##########################################################
s = sky(0.05) # optical depth of atmospheric layer
theta0 = np.pi/3 # sun's zenith angle

x,y = np.meshgrid(np.linspace(-1,0,41), np.linspace(-1,1,81))
th,phi = np.pi/2*np.sqrt(x**2+y**2), np.arctan2(x,y) # clockwise
P = np.full_like(x, np.nan)

for i in range(phi.shape[0]):
    for j in range(phi.shape[1]):
        if th[i,j]>np.pi/2: continue
        I = s.transmission(theta0, phi[i,j])
        I = [np.interp(th[i,j], s.theta[::-1], I[::-1]) for I in I] 
        P[i,j] = np.sqrt(I[1]**2 + I[2]**2)/I[0] # degree of polarization

plt.figure(figsize=(5,5))
plt.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0)
plt.axis('equal')
plt.axis('off')

z = np.exp(1j*np.linspace(np.pi/2, 3*np.pi/2, 32))
levels = np.linspace(0,1,26)

cs = plt.contour(x,y,P,levels=levels)
plt.plot(np.real(z), np.imag(z), 'k', lw=1)
plt.vlines(0,-1,1,lw=1)
plt.text(-1, 0.95, r'$\tau_0 = 0.05$', fontsize=14)
plt.text(-1, 0.85, r'$\theta_0 = 60^{\circ}$', fontsize=14)
plt.clabel(cs, fmt='%1.2f')

# fig3.2 ##########################################################
theta0 = 80/180*np.pi # sun's zenith angle

x,y = np.meshgrid(np.linspace(0,1,41), np.linspace(-1,1,81))
th,phi = np.pi/2*np.sqrt(x**2+y**2), np.arctan2(x,y) # clockwise
P = np.full_like(x, np.nan)

for i in range(phi.shape[0]):
    for j in range(phi.shape[1]):
        if th[i,j]>np.pi/2: continue
        I = s.transmission(theta0, phi[i,j])
        I = [np.interp(th[i,j], s.theta[::-1], I[::-1]) for I in I] 
        P[i,j] = np.sqrt(I[1]**2 + I[2]**2)/I[0] # degree of polarization

z = np.exp(1j*np.linspace(-np.pi/2, np.pi/2, 32))
levels = np.linspace(0,1,26)

cs = plt.contour(x,y,P,levels=levels)
plt.plot(np.real(z), np.imag(z), 'k', lw=1)
plt.text(1, 0.95, r'$\tau_0 = 0.05$', ha='right', fontsize=14)
plt.text(1, 0.85, r'$\theta_0 = 80^{\circ}$', ha='right', fontsize=14)
plt.clabel(cs, fmt='%1.2f')

plt.tight_layout()
plt.savefig('fig3.eps')
plt.show()

plt.close()
# fig4.1 ##########################################################
s = sky(0.2) # optical depth of atmospheric layer
theta0 = np.pi/3 # sun's zenith angle

x,y = np.meshgrid(np.linspace(-1,0,41), np.linspace(-1,1,81))
th,phi = np.pi/2*np.sqrt(x**2+y**2), np.arctan2(x,y) # clockwise
P = np.full_like(x, np.nan)

for i in range(phi.shape[0]):
    for j in range(phi.shape[1]):
        if th[i,j]>np.pi/2: continue
        I = s.transmission(theta0, phi[i,j])
        I = [np.interp(th[i,j], s.theta[::-1], I[::-1]) for I in I] 
        P[i,j] = np.sqrt(I[1]**2 + I[2]**2)/I[0] # degree of polarization

plt.figure(figsize=(5,5))
plt.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0)
plt.axis('equal')
plt.axis('off')

z = np.exp(1j*np.linspace(np.pi/2, 3*np.pi/2, 32))
levels = np.linspace(0,1,26)

cs = plt.contour(x,y,P,levels=levels)
plt.plot(np.real(z), np.imag(z), 'k', lw=1)
plt.vlines(0,-1,1,lw=1)
plt.text(-1, 0.95, r'$\tau_0 = 0.2$', fontsize=14)
plt.text(-1, 0.85, r'$\theta_0 = 60^{\circ}$', fontsize=14)
plt.clabel(cs, fmt='%1.2f')

# fig4.2 ##########################################################
theta0 = 80/180*np.pi # sun's zenith angle

x,y = np.meshgrid(np.linspace(0,1,41), np.linspace(-1,1,81))
th,phi = np.pi/2*np.sqrt(x**2+y**2), np.arctan2(x,y) # clockwise
P = np.full_like(x, np.nan)

for i in range(phi.shape[0]):
    for j in range(phi.shape[1]):
        if th[i,j]>np.pi/2: continue
        I = s.transmission(theta0, phi[i,j])
        I = [np.interp(th[i,j], s.theta[::-1], I[::-1]) for I in I] 
        P[i,j] = np.sqrt(I[1]**2 + I[2]**2)/I[0] # degree of polarization

z = np.exp(1j*np.linspace(-np.pi/2, np.pi/2, 32))
levels = np.linspace(0,1,26)

cs = plt.contour(x,y,P,levels=levels)
plt.plot(np.real(z), np.imag(z), 'k', lw=1)
plt.text(1, 0.95, r'$\tau_0 = 0.2$', ha='right', fontsize=14)
plt.text(1, 0.85, r'$\theta_0 = 80^{\circ}$', ha='right', fontsize=14)
plt.clabel(cs, fmt='%1.2f')

plt.tight_layout()
plt.savefig('fig4.eps')
plt.show()

plt.close()
# fig5.1 ##########################################################
s = sky(0.5) # optical depth of atmospheric layer
theta0 = np.pi/3 # sun's zenith angle

x,y = np.meshgrid(np.linspace(-1,0,41), np.linspace(-1,1,81))
th,phi = np.pi/2*np.sqrt(x**2+y**2), np.arctan2(x,y) # clockwise
P = np.full_like(x, np.nan)

for i in range(phi.shape[0]):
    for j in range(phi.shape[1]):
        if th[i,j]>np.pi/2: continue
        I = s.transmission(theta0, phi[i,j])
        I = [np.interp(th[i,j], s.theta[::-1], I[::-1]) for I in I] 
        P[i,j] = np.sqrt(I[1]**2 + I[2]**2)/I[0] # degree of polarization

plt.figure(figsize=(5,5))
plt.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0)
plt.axis('equal')
plt.axis('off')

z = np.exp(1j*np.linspace(np.pi/2, 3*np.pi/2, 32))
levels = np.linspace(0,1,26)

cs = plt.contour(x,y,P,levels=levels)
plt.plot(np.real(z), np.imag(z), 'k', lw=1)
plt.vlines(0,-1,1,lw=1)
plt.text(-1, 0.95, r'$\tau_0 = 0.5$', fontsize=14)
plt.text(-1, 0.85, r'$\theta_0 = 60^{\circ}$', fontsize=14)
plt.clabel(cs, fmt='%1.2f')

# fig5.2 ##########################################################
theta0 = 80/180*np.pi # sun's zenith angle

x,y = np.meshgrid(np.linspace(0,1,41), np.linspace(-1,1,81))
th,phi = np.pi/2*np.sqrt(x**2+y**2), np.arctan2(x,y) # clockwise
P = np.full_like(x, np.nan)

for i in range(phi.shape[0]):
    for j in range(phi.shape[1]):
        if th[i,j]>np.pi/2: continue
        I = s.transmission(theta0, phi[i,j])
        I = [np.interp(th[i,j], s.theta[::-1], I[::-1]) for I in I] 
        P[i,j] = np.sqrt(I[1]**2 + I[2]**2)/I[0] # degree of polarization

z = np.exp(1j*np.linspace(-np.pi/2, np.pi/2, 32))
levels = np.linspace(0,1,26)

cs = plt.contour(x,y,P,levels=levels)
plt.plot(np.real(z), np.imag(z), 'k', lw=1)
plt.text(1, 0.95, r'$\tau_0 = 0.5$', ha='right', fontsize=14)
plt.text(1, 0.85, r'$\theta_0 = 80^{\circ}$', ha='right', fontsize=14)
plt.clabel(cs, fmt='%1.2f')

plt.tight_layout()
plt.savefig('fig5.eps')
plt.show()
