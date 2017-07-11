"""
Created on Thu Jun 29 21:56:12 2017

@author: champ4
"""
import unittest
import simulator
import numpy as np

class TestFunctions(unittest.TestCase):
    def test_measurement_eq(self):
        s = np.array([[1.,0,0],[0.,1,0],[0.,0,1]]) # (3 pointings, xyz)
        A = np.array([[1.,2],[.5,1.], [.25,.5]]) # (3 pointings, 2 fqs)
        I = np.array([[1.,1],[1,1], [1,1]]) # (3 pointings, 2 fqs)
        bls = np.array([[10.,0,0],[0.,10,0],[0.,0,10]]) # (3 bls, xyz)
        fqs = np.array([.1e9,.2e9]) # (2 fqs) in Hz
        vis = simulator.measurement_eq(A, I, s, bls, fqs)
        self.assertEqual(vis.shape, (3,2)) # (3 bls, 2 fqs)
        self.assertTrue(np.all(np.iscomplex(vis))) # (3 bls, 2 fqs)
        
        # XXX test actual values (
        """Done, but only if I set the Almost Equal when the place of
        significance to about 5 due to imaginary numbers not fully matching, 
        this is most noticeable for imaginary numbers almost 0 *which the AlmostEqual
        takes into account and does neglect this difference* but for '4th' test
        it was off by the 5th place value"""
        
        A = np.array([[0.75,1.00,0.75]])
        I = np.array([[150.0,100.0,75.0]])
        s = np.array([[0,1.0,0],[1.0,0,1.0]])
        
        fqs = np.array([100.0e6,150.0e6,200.0e6])
        bls = np.array([[0.0,10.0,0.0],[1.0,7.0,0.0]])

        vis = simulator.measurement_eq(A, I, s, bls, fqs)
        
        
        self.assertAlmostEqual(np.real(vis[0,0]), (56.250  ))
        self.assertAlmostEqual(np.real(vis[0,1]), (200.000 ))
        self.assertAlmostEqual(np.real(vis[0,2]), (28.125  ))
        self.assertAlmostEqual(np.real(vis[1,0]), (-112.500)) #4th test
        self.assertAlmostEqual(np.real(vis[1,1]), (-200.000))
        self.assertAlmostEqual(np.real(vis[1,2]), (-56.250 ))
        
        self.assertAlmostEqual(np.absolute(vis[0,0]), (112.50000000000031))
        self.assertAlmostEqual(np.absolute(vis[0,1]), (200.0             ))
        self.assertAlmostEqual(np.absolute(vis[0,2]), (56.249999999999687))
        self.assertAlmostEqual(np.absolute(vis[1,0]), (224.99999999999997)) #4th test
        self.assertAlmostEqual(np.absolute(vis[1,1]), (200.000           ))
        self.assertAlmostEqual(np.absolute(vis[1,2]), (112.5             )) 
        
        self.assertEqual(vis.shape, (2,3)) # (2 bls, 3 fqs)
        
#        self.assertAlmostEqual/
        #First run had issues, it seems like there is no i(imaginary values)
        

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
