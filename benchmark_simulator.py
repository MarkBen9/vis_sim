import simulator
import aipy
import numpy as np

#NANTS = 350
NANTS = 16
NBLS = NANTS * (NANTS+1) / 2 # for HERA
#NSIDE = 512
NSIDE = 64
#NFREQS = 1024
NFREQS = 8

h = aipy.healpix.HealpixMap(nside=NSIDE)
s = np.array(h.px2crd(np.arange(h.npix()))).T
I = np.random.normal(size=(h.npix(),NFREQS)) # sky is real-valued
A = np.random.normal(size=(h.npix(),NFREQS)) + \
    1j*np.random.normal(size=(h.npix(),NFREQS)) # antenna response is complex
fqs = np.linspace(100e6,200e6,NFREQS)
antpos = np.random.normal(size=(NANTS,3))
bls = np.array([ai-aj for i,ai in enumerate(antpos) for aj in antpos[i:]])

vis = simulator.measurement_eq(A, I, s, bls, fqs)
