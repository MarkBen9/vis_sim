"""
Created on Thu Jun 29 21:56:12 2017

@author: champ4
"""

import unittest
import simulator
import numpy as np

class TestSource(unittest.TestCase):

    
    def test_init(self):
    
        s = simulator.Source(ra=np.pi/2, dec=0., jansky=1., index=-1.5, mfreq=100.)
        self.assertEqual(s.ra, np.pi/2)
        self.assertEqual(s.dec, 0.)
        self.assertEqual(s.jansky, 1.)
        self.assertEqual(s.index, -1.5)
        self.assertEqual(s.mfreq, 100.)
    
    
    def test_get_source_vector(self):
        
        
        s = simulator.Source(ra=0., dec=0.)
        # Equator, source overhead
        lst, lat = 0., 0.
        xyz_top = s.get_source_vector(lst, lat)
        self.assertEqual(xyz_top[0], 0.)
        self.assertEqual(xyz_top[1], 0.)
        self.assertEqual(xyz_top[2], 1.)
        
        s = simulator.Source(ra=0., dec=np.pi/2)
        # Test declination change
        lst, lat = 0., 0.
        xyz_top = s.get_source_vector(lst, lat)
        self.assertAlmostEqual(xyz_top[0], 0.)
        self.assertAlmostEqual(xyz_top[1], 1.)
        self.assertAlmostEqual(xyz_top[2], 0.)
        #XXX test hour angle change
        
        s = simulator.Source(ra=0.0, dec=0.0)
        # Test declination change
        lst, lat = np.pi/-8.0, 0.
        xyz_top = s.get_source_vector(lst, lat)
        self.assertAlmostEqual(xyz_top[0], -np.sin(np.pi/-8.0))
        self.assertAlmostEqual(xyz_top[1], 0.)
        self.assertAlmostEqual(xyz_top[2], np.cos(np.pi/-8.0))
        #PROBABLY change lst slowly
  
    def test_measurement_eq(self):
        s=simulator.Source(ra = 0.0, dec = 0.0)
#        


if __name__ == '__main__':
    unittest.main()
