import numpy as np
from numpy.linalg import solve,eig
from scipy.special import roots_legendre
from rayleigh import phase_func

class sky:
    """ solve radiative transfer equation for the
    sun-lit sky by the method of discrete ordinates
    reference:
      S. Chandrasekhar
        "Radiative Transfer" sections 13, 48
      R. M. Goody and Y. L. Yung
        "Atmospheric Radiation" 2nd edition, section 8.2
    """
    def __init__(self, tau0, theta0, N=256):
        """
        tau0 = optical depth of atmospheric layer
        theta0 = zenith angle of sun (radian)
        N = number of discrete ordinates / 2
        """
        mu0 = np.cos(theta0)
        mu,w = roots_legendre(2*N)
        mu,w = mu[N:],w[N:]
        mu_ = mu[:,np.newaxis]

        M = (np.eye(N) - w*phase_func(mu_,mu))/(mu_*mu)
        h,g = eig(M)
        # swap zero eigenvalue to head
        j = np.argmin(h[0])
        h[0,j],h[0,0] = h[0,0],1
        g[0,:,j],g[0,:,0] = g[0,:,0],g[0,:,j]

        k = np.sqrt(h)
        f = g/mu_/k[:,np.newaxis,:]
        g = np.asarray([f+g, f-g])/2
        M -= np.eye(N)/mu0**2
        s = -phase_func(mu,-mu0)/mu/mu0/2 # / solar_const/pi
        h = solve(M,s)
        f = -h/mu*mu0
        h = np.asarray([f+h, f-h])/2

        k[0,0] = 0
        g[1,1] = -g[1,1]
        h[1,1] = -h[1,1]

        f = g[1]*np.exp(-k[:,np.newaxis,:]*tau0)
        M = np.asarray([g[0]+f, g[0]-f])
        M[0,0,:,0] = 1
        M[1,0,:,0] = mu + tau0/2
        f = -h[0]*np.exp(-tau0/mu0)
        v = np.asarray([f-h[1], f+h[1]])
        s = solve(M,v)
        L = np.asarray([s[0]+s[1], s[0]-s[1]])/2

        L[1,0,0] = s[1,0,0]/2
        L[0,0,0] = (s[0,0,0] - L[1,0,0]*tau0)/2

        self.theta = np.arccos(mu)
        self.tau0 = tau0
        self.mu0 = mu0
        self.mu = mu
        self.k = k
        self.g = g
        self.h = h
        self.L = L

    def intensity(self, phi, tau, refl):
        """
        phi = azimuthal angle from sun's direction
        tau = optical depth at observer (increase downward)
        refl = 1 or 0 if rays go up or down, resp.
        return I = ray's intensity / solar_const/pi
        I.shape = (N,); I[i] corresponds to self.theta[i]
        """
        t = self.L*np.exp([self.k*(tau - self.tau0), -self.k*tau])
        f = self.g if refl else np.asarray([self.g[1], self.g[0]])
        f[0,0,:,0] = 1
        f[1,0,:,0] = (tau + self.mu if refl else tau - self.mu)
        I = self.h[0 if refl else 1]*np.exp(-tau/self.mu0)
        I += np.einsum('kmij,kmj->mi', f, t)
        I[1:] *= 2*np.cos([phi, 2*phi])[:,np.newaxis]
        return np.sum(I,0)

    def transmission(self, phi, tau):
        return self.intensity(phi, tau, 0)

    def reflection(self, phi, tau):
        return self.intensity(phi, tau, 1)


if __name__ == '__main__':# example usage
    import matplotlib.pyplot as plt
    s = sky(0.2, np.pi/3)
    plt.plot(s.theta, s.transmission(0, 0.2), 'b')
    plt.plot(-s.theta, s.transmission(np.pi, 0.2), 'b')
    plt.show()
