import numpy as np
import matplotlib.pyplot as plt

class Source:
    def __init__(self, ra, dec, jansky=100., index=-1., mfreq=150.):
        '''ra: right ascension (in radians)
        dec: declination (in radians)
        jansky: source flux density at specified frequency (mfreq)
        index: spectral index of flux density
        mfreq: specified frequency of jansky in MHz'''
        self.jansky = jansky
        self.mfreq = mfreq
        self.index = index
        self.ra = ra
        self.dec = dec
    def get_source_vector(self, lst, lat):
        '''Return topocentric vector pointing toward source given local 
        sidereal time (lst) and latitude (lat) in radians.'''
        ha = lst - self.ra
        x_eq = np.cos(self.dec) * np.cos(ha)
        y_eq = -np.sin(ha)
        z_eq = np.sin(self.dec)
        xyz_eq = np.array([x_eq, y_eq, z_eq])
        M_eq2top = np.array([[0., 1., 0],
                             [-np.sin(lat), 0, np.cos(lat)],
                             [ np.cos(lat), 0, np.sin(lat)]])
        xyz_top = np.dot(M_eq2top, xyz_eq)
        return xyz_top
        
