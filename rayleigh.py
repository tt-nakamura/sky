import numpy as np

def phase_func(mu1,mu2):
    """ angular distrubution of Rayleigh scattering
    mu1 = direction cosine of incident ray
    mu2 = direction cosine of scattered ray
    return [p_1,p_2,p_3] where p_m (m=0,1,2) are
    coeff of Fourier cosine series in azimuthal angle
    """
    m11,m22 = mu1*mu1, mu2*mu2
    a,b = m11*m22, m11+m22
    p = a-b+1
    return np.asarray([
            1 + 1.125*(a-b/3+1/9),
            0.75*mu1*mu2*np.sqrt(p),
            0.1875*p])

def phase_matrix(mu1,mu2):
    """ angular distrubution of Rayleigh scattering
    for polarized light; Stokes parameters are
    (I,Q) for m=0, (I,U) for m=1, (I,Q,U) for m=2
    U=0 (m=0), Q=I (m=1), V=0 (m=0,1,2)
    return P_1,P_2,P_3 where P_m (m=0,1,2)
    are coeff matrices of Fourier
    cosine series for I,Q, sine series for U
    P_1,P_2 are 2x2 matrices; P_3 is 3x3 matrix
    reference:
      S. Chandrasekhar
        "Radiative Transfer" sections 17, 69
    """
    m11,m12,m22 = mu1*mu1, mu1*mu2, mu2*mu2
    m112,m122 = m11*mu2, mu1*m22
    a,b,c = m11*m22, m11+m22, m11-m22
    p = a-b+1
    s = 0.75*np.sqrt(p)
    return np.asarray([
            [1 + 1.125*(a-b/3+1/9),
             0.375*(1-2*m11-b+3*a)],
            [0.375*(1-2*m22-b+3*a),
             1.125*p]
            ]),  np.asarray([
            [s*m12, s*mu1], [s*mu2, s]
            ]),  np.asarray([
            [0.1875*p,
             0.1875*(a+c-1),
             0.375*(m112-mu2)],
            [0.1875*(a-c-1),
             0.1875*(a+b+1),
             0.375*(m112+mu2)],
            [0.375*(m122-mu1),
             0.375*(m122+mu1),
             0.75*m12]])
