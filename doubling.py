import numpy as np
from numpy.linalg import inv
from scipy.special import roots_legendre
from scipy.interpolate import interp1d
from rayleigh import phase_func

def doubling(R,T,M):
    """ double optical depth of atmospheric layer
    R,T = reflection and transmission matrix, resp.
    M = number of douling times
    reference:
      R. M. Goody and Y. L. Yung
        "Atmospheric Radiation" 2nd edition, section 8.3
    """
    for _ in range(M):
        A = np.eye(R.shape[-1]) - np.matmul(R,R)
        B = inv(A)
        A = np.matmul(T,B)
        B = np.matmul(A,R)
        R += np.matmul(B,T)
        T = np.matmul(A,T)
    return R,T

class sky:
    """ solve radiative transfer equation for the
    sun-lit sky by doubling method
    """
    def __init__(self, tau, N=128, M=20):
        """
        tau = optical depth of atmospheric layer
        N = number of discrete ordinates / 2
        M = number of doubling times
        """
        mu,w = roots_legendre(2*N)
        mu,w = mu[N:],w[N:]
        mu_ = mu[:,np.newaxis]
        tau1 = np.ldexp(tau, -M)

        R = tau1*w*phase_func(mu_,-mu)/2/mu_
        T = np.asarray([R[0], -R[1], R[2]])
        T += np.diag(np.exp(-tau1/mu))

        R,T = doubling(R,T,M)# doubling M times

        T -= np.diag(np.exp(-tau/mu))
        R = interp1d(mu, R/w, 'cubic', fill_value='extrapolate')
        T = interp1d(mu, T/w, 'cubic', fill_value='extrapolate')

        self.theta = np.arccos(mu)
        self.R = R
        self.T = T

    def intensity(self, theta0, phi, refl):
        """
        theta0 = zenith angle of sun (radian)
        phi = azimuthal angle from sun's direction
        return I = ray's intensity / solar_const/pi
        I.shape = (N,); I[i] corresponds to self.theta[i]
        if refl==1 (or 0), rays go up (or down)
        from layer's top (or bottom, resp.)
        """
        I = (self.R if refl else self.T)(np.cos(theta0))/2
        I[1:] *= 2*np.cos([phi, 2*phi])[:,np.newaxis]
        return np.sum(I,0)

    def transmission(self, theta0, phi):
        return self.intensity(theta0, phi, 0)

    def reflection(self, theta0, phi):
        return self.intensity(theta0, phi, 1)


if __name__ == '__main__':# example usage
    import matplotlib.pyplot as plt
    s = sky(0.2)
    plt.plot(s.theta, s.transmission(np.pi/3, 0), 'b')
    plt.plot(-s.theta, s.transmission(np.pi/3, np.pi), 'b')
    plt.show()
