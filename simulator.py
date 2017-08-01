"""
Created on Tue Jun 27 13:47:35 2017

@author: champ4
"""
import numpy as np
#import matplotlib.pyplot as plt

C = 3.0e8 # speed of light in m/s

def measurement_eq(A, I, s, bls, fqs):
    '''measurment_eq(A, I, s, *M, *Hz)
    
    ''''''
    A: (s,fq) 2-dim array
    I: (s,fq) 2-dim array
    bls: (baselines,xyz)
    s: (sky position,xyz)
    '''
    b_s = np.dot(bls, s.T).T
    b_s.shape = (b_s.shape[0],1,b_s.shape[1]) # shape is (s,fq,bl)
    fqs.shape = (1,fqs.size,1)
    if A.ndim < 3:  
        A.shape = A.shape + (1,)
    if I.ndim < 3:
        I.shape = I.shape + (1,)
    V = np.sum(A * I * np.exp(-2j*np.pi / C * fqs * b_s), axis=0).T
    return V # shape is (bls,fqs)

def measurement_eq_split(A, I, s, ant_pos, fqs):
    '''A: (s,fq) 2-dim array
    I: (s,fq) 2-dim array
    ant_pos: (antenna's position,xyz)
    s: (sky position,xyz)
    fqs: (frequency) 1-dim array'''
#selection of points of interest
    s_valid = s[:,2] > 0
    s_above = s[s_valid,:]
    A_above = A[s_valid,:]
    I_above = I[s_valid,:]
#dot product
    a_s = np.dot(ant_pos, s_above.T).T         
#Matrix set up
    a_s.shape = (a_s.shape[0],1,a_s.shape[1]) # shape is (s,fq,ant_pos)
    fqs.shape = (1,fqs.size,1)
    if A_above.ndim < 3:  
        A_above.shape = A_above.shape + (1,)
    if I_above.ndim < 3:
        I_above.shape = I_above.shape + (1,)
#Equation for just an antenna
    V_a = (np.sqrt(A_above * I_above) * np.exp(-2j*np.pi / C * fqs * a_s))

    V_a_conj = V_a.conj()
#procduct of two antenna-baseline
    V = [np.sum(V_a[...,i]*V_a_conj[...,j], axis=0 ) 
            for i in range(V_a.shape[-1]) for j in range(i,V_a.shape[-1])]
    
    return np.array(V) #shape (bls, fqs)
   
    
def HealPix_Res_Map(A, I, s, fqs, Err = 1.0e-6, Nside = 512, bls=14):   

    #Things to Ask Aaron Parsons
# -Should Error be not complex?
# -When Finding Error, Should the Nside be dependent on what resolution we are working with. 
    #i.e.(2*pi/Nside) '$Sigma$' (i is an element of a) A_sub_i * I_sub_i - A_sub_a*I_sub_a / 4 + $Sigma$ Err_sub_i
    # -- If A_sub_i is of Nside 256, then should Nside be 256 as well, or should it be based on the highest Nside?
# -So far I feel like I should invest into making this its own class, from how it is running it might be cleaner/
    #After thoughts: work first on just setting it up then consider the option. but for class I think Error,
    #/n Nside and Bls would be the Init
    
    '''HealPix_Res_Map(A , I , s , fqs , Err = 1.0e-6, Nside = 512, bls = 14)
A(s,fqs) numpy.array, sources direction orginzed by Healpix nested 
I(s,fqs) numpy.array, sources strength  orginzed by Healpix nested
s(sky postion, xyz) np.array Healpix nested
error
Nside, Preferable if entering a values that are powers of 2 and divisible by 4
https://lambda.gsfc.nasa.gov/toolbox/tb_pixelcoords.cfm
!!!NESTED!!!
'''

#Creating an empty Array for future use    
    A_Nside_256 = np.zeros((A.shape[0]/4, A.shape[1]))
    A_Nside_128 = np.zeros((A_Nside_256.shape[0]/4,A_Nside_256.shape[1]))
#----------
    length1 =A_Nside_256.shape[0]
    length2 =A_Nside_128.shape[0]
#----------
    I_Nside_256 = np.zeros((I.shape[0]/4, I.shape[1]))
    I_Nside_128 = np.zeros((I_Nside_256.shape[0]/4,I_Nside_256.shape[1]))
#----------
    Error_256 = np.zeros(A_Nside_256.shape)
    Error_128 = np.zeros(A_Nside_128.shape)
#Populating Arrays    
    I_Nside_256 = np.array([np.sum(I[4*i:4*(i+1),...], axis=0 ) 
            for i in range(length1)])
    I_Nside_128 =np.array([np.sum(I_Nside_256[4*i:4*(i+1),...], axis=0 )
            for i in range(length2)])
#----------
    A_Nside_256 =np.array([np.sum(A[4*i:4*(i+1),...], axis=0 ) /4.0
            for i in range(length1)])
    A_Nside_128 =np.array([np.sum(A_Nside_256[4*i:4*(i+1),...], axis=0 ) /4.0
            for i in range(length2)])
#Populating Error Array for first step
    Error_256 = np.array([np.sum(A[4*a:4*(a+1),...]*I[4*a:4*(a+1),...] - (A_Nside_256[a,...]*
                 I_Nside_256[a,...]/4), axis = 0) for a in range(Error_256.shape[0])])
    Error_256=Error_256*(2*np.pi/Nside)
#----------
    Error_128 = np.array([np.sum(A_Nside_256[4*a:4*(a+1),...]*I_Nside_256[4*a:4*(a+1),...] - (A_Nside_128[a,...]*
                 I_Nside_128[a,...]/4), axis = 0) for a in range(Error_128.shape[0])])
    Error_128=Error_128*(4*np.pi/Nside)
#Creating of a boolean Array that will be the basis for the final array


#Boolean Arrays (TBC)

#System TBC

#Selection of Chunks from Boolean Arrays
   
    
    """
    a_s = np.dot(ant_pos, s_above.T).T         
    #Matrix set up
    a_s.shape = (a_s.shape[0],1,a_s.shape[1]) # shape is (s,fq,ant_pos)
    fqs.shape = (1,fqs.size,1)
    #Checks dimension of arrays
    if A_above.ndim < 3:  
        A_above.shape = A_above.shape + (1,)
    if I_above.ndim < 3:
        I_above.shape = I_above.shape + (1,)
    #Equation for just an antenna
   
    V_a = (np.sqrt(A_above) * np.sqrt(I_above) * np.exp(-2j*np.pi / C * fqs * a_s))
    V_a = V_a.astype(np.complex64)
    #Repeat V_a but for "less" sources.
    V_a_conj = V_a.conj()
    
    
    print(V_a.dtype)
    #procduct of two antenna to give a baseline
    
    V = [np.sum(V_a[...,i]*V_a_conj[...,j], axis=0 ) 
            for i in range(V_a.shape[-1]) for j in range(i,V_a.shape[-1])]
    """
    #return np.array(V) #shape (bls, fqs)
    

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