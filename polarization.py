import numpy as np
from scipy.special import roots_legendre
from scipy.interpolate import interp1d
from rayleigh import phase_matrix
from doubling import doubling

class sky:
    """ solve radiative transfer equation for the
    polarization of sun-lit sky by doubling method
    reference:
      R. M. Goody and Y. L. Yung
        "Atmospheric Radiation" 2nd edition, section 8.3
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

        R,T = [],[]
        # for Fourier series in azimuthal angle
        for m,P in enumerate(phase_matrix(mu_,mu)):
            T1 = tau1*w/2/mu_*P
            R1 = -T1 if m==1 else T1#.copy()
            L,K = P.shape[0], N*P.shape[0]
            R1 = np.swapaxes(R1,1,2).reshape(K,K)
            T1 = np.swapaxes(T1,1,2).reshape(K,K)
            T1 += np.diag(np.tile(np.exp(-tau1/mu), L))

            R1,T1 = doubling(R1,T1,M)# doubling M times

            R1 = R1[:,:N].reshape(L,N,N)
            T1 = T1[:,:N].reshape(L,N,N)
            T1[0] -= np.diag(np.exp(-tau/mu))
            R.append(interp1d(mu, R1/w, 'cubic', fill_value='extrapolate'))
            T.append(interp1d(mu, T1/w, 'cubic', fill_value='extrapolate'))

        self.theta = np.arccos(mu)
        self.tau = tau
        self.mu = mu
        self.w = w
        self.R = R
        self.T = T

    def ground(self, albedo, mu0, refl):
        """ reflection from the ground
        assuming isotropic and unpolarized
        albedo = reflection coefficient
        mu0 = cos(zenith angle of sun)
        return correction to Stokes parameters (I,Q)
          due to reflection (shape(2,N)) / solar_const/pi
        reference:
          S. Chandrasekhar
            "Radiative Transfer" section 72
        """
        wmu = self.w * self.mu
        F = np.dot(np.dot(self.R[0](self.mu)[0], self.w), wmu)
        I = np.dot(self.T[0](mu0)[0], wmu)
        I += mu0*np.exp(-self.tau/mu0)
        I *= albedo/(1 - 2*albedo*F)
        A = self.T if refl else self.R
        J = I*np.dot(A[0](self.mu), self.w)
        if refl: J += I*np.exp(-self.tau/mu)
        return J

    def intensity(self, theta0, phi, refl, albedo=0):
        """
        theta0 = zenith angle of sun (radian)
        phi = azimuthal angle from sun's direction
        albedo = reflection coefficient of the ground
        return I = Stokes parameters (I,Q,U) / solar_const/pi
        I.shape = (3,N); I[i] corresponds to (I,Q,U)
        I[:,j] corresponds to self.theta[j]
        if refl==1 (or 0), rays go up (or down)
        from layer's top (or bottom, resp.)
        assume incident sun light is unpolarized
        """
        mu0 = np.cos(theta0)
        I = np.zeros((3, len(self.theta)))
        for m,A in enumerate(self.R if refl else self.T):
            I1 = A(mu0)/2
            if m:
                I1[:-1] *= 2*np.cos(m*phi)
                I1[-1]  *= 2*np.sin(m*phi)
            if m==0:   I[:-1] += I1
            elif m==1: I += [I1[0], I1[0], I1[1]]
            else:      I += I1

        if albedo:
            I[:-1] += self.ground(albedo, mu0, refl)

        return I

    def transmission(self, theta0, phi, albedo=0):
        return self.intensity(theta0, phi, 0, albedo)

    def reflection(self, theta0, phi, albedo=0):
        return self.intensity(theta0, phi, 1, albedo)

if __name__ == '__main__': # example usage
    import matplotlib.pyplot as plt
    s = sky(0.2)
    plt.plot(s.theta, s.transmission(np.pi/3, 0).T)
    plt.plot(-s.theta, s.transmission(np.pi/3, np.pi).T)
    plt.show()
