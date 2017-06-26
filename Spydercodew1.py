#! /usr/bin/env python
"""
Spyder Editor
Marcus U. B. 
Prj week 1:2, Measurement EQ Time started ()???? 6/19/2017
Prj week 1:3, cont. (9:30) 6/20/2017
"""
import numpy as np
import matplotlib.pyplot as plt

#[x,y,z].
# a,b,c are position vectors.
#b; a and c antenna's.

def sme(v,a,b,c):
    t=[]
    for i in range(len(a)):
        t.append(a[i]-c[i])    
    x = 100.0*(150.0/v)**(1)*np.exp(-2*np.pi*np.dot(t,b)*v/(3e8))
    return x

emp=[]
horz=[]

for i in range(181):
    horz.append(i)
    emp.append(sme(100,[0,0,0],[np.cos(i/180.0*np.pi),0,np.sin(i/180*np.pi)],[50,0,0]))

plt.figure(1)
plt.title('test')
plt.xlabel('Theta (Degrees)')
plt.ylabel('Jy')
plt.plot(horz,emp)
#plt.axis([0,180,149.98,150.02])

plt.show()
#To expand it to 2 or 3 repeaded the forloop and combined the array
