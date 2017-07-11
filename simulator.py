"""
Created on Tue Jun 27 13:47:35 2017

@author: champ4
"""
import numpy as np
import matplotlib.pyplot as plt

C = 3.0e8 # speed of light in m/s

def measurement_eq(A, I, s, bls, fqs):
    '''A: (s,fq) 2-dim array
    I: (s,fq) 2-dim array
    bls: (baselines,xyz)
    s: (sky position,xyz)'''
    b_s = np.dot(bls, s.T).T
    b_s.shape = (b_s.shape[0],1,b_s.shape[1]) # shape is (s,fq,bl)
    fqs.shape = (1,fqs.size,1)
    A.shape = A.shape + (1,)
    I.shape = I.shape + (1,)
    V = np.sum(A * I * np.exp(-2j*np.pi / C * fqs * b_s), axis=0).T
    return V # shape is (bls,fqs)

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
    
    def measurement_eq(self, antenna1 ,antenna2, lst, lat):
        t= np.array(antenna1) - np.array(antenna2)   
        x = self.jansky*((self.mfreq)/(150.0))**(self.index)*np.exp(-2j*
                np.pi*np.dot(t,self.get_source_vector(lst,lat))*self.mfreq*1e6/(3e8))
        return x