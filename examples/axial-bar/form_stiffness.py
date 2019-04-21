#!/usr/bin/env python
import numpy as np
import math
import numpy as np
np.set_printoptions(precision=4 , suppress=True)
from basis import N, Nprime

class KnotParameter:    
    """
    Class for defining a parameter that has certain range based on
    element.

    Author: Komahan Boopathy
    """
    def __init__(self, enum, bounds):
        """
        Construction logic for the object
        """
        self.enum = enum
        self.bounds = bounds
        return
    
    def __str__(self):
        """
        Define behavior within print() call
        """
        return str("enum : " + str(self.enum)) + \
               str("\nbounds : " + " ".join([str(i) for i in self.bounds]))

    def get_quadrature_points_weights(self, npoints):
        """
        Returns an arrays of quadrature points and weights in the
        original element parameter domain. The conversion from
        original parameter domain to the element parameter domain is
        done here.
        """
        
        # This is based on  weight 1.0 on interval [-1,1]
        xihat, wxihat = np.polynomial.legendre.leggauss(npoints)
        a = self.bounds[0]
        b = self.bounds[1]
        
        # scale weights
        wxi = (b-a)*wxihat/2.
        
        # transformation of variables into element space
        xi = (b-a)*xihat/2. + (b+a)/2.

        return xi, wxi
    
if __name__ == "__main__":   

    xi1 = KnotParameter(1, [0, 0.5])
    ## xi1pts, wxi1 = xi1.get_quadrature_points_weights(3)
    ## print xi1pts, sum(wxi1)

    e1jac = np.zeros((3,3))
    for i in range(0,3):
        for j in range(0,3):
            xi1pts, wxi1 = xi1.get_quadrature_points_weights(3)
            for xie, we in zip(xi1pts, wxi1):
                e1jac[i,j] += Nprime(xie,i+1,xi1.enum)*Nprime(xie,j+1,xi1.enum)*we/2  
    print e1jac

    xi2 = KnotParameter(2, [0.5, 1.0])
    ## xi2pts, wxi2 = xi2.get_quadrature_points_weights(3)
    ## print xi2pts, sum(wxi2)

    e2jac = np.zeros((3,3))
    for i in range(0,3):
        for j in range(0,3):
            xi2pts, wxi2 = xi2.get_quadrature_points_weights(3)            
            for xie, we in zip(xi2pts, wxi2):
                e2jac[i,j] += Nprime(xie,i+1,xi2.enum)*Nprime(xie,j+1,xi2.enum)*we/2      
    print e2jac

    print "global stiffness matrix"

    # Assemble Global Stiffnes matrix
    K = np.zeros((4,4))
    K[0:3,0:3] += e1jac[:,:]
    K[1:4,1:4] += e2jac[:,:]    
    print K
    print ""

    K = np.delete(K, 0, 1)
    K = np.delete(K, 0, 0)
    print K

    F = np.array([1.0/3.0, 2.0/3.0, 2.0/3.0, 1.0/3.0])
    print F    
    F = np.delete(F,0,0)
    print F

    U = np.linalg.solve(K, F)
    print U
    
    fe1 = np.zeros((3))    
    for i in range(0,3):
        xi2pts, wxi2 = xi1.get_quadrature_points_weights(3)
        for xie, we in zip(xi2pts, wxi2):
            fe1[i] += N(xie, i+1, xi1.enum)*we*2
    print "fe1", fe1, np.sum(fe1)

    fe2 = np.zeros((3))    
    for i in range(0,3):
        xi2pts, wxi2 = xi2.get_quadrature_points_weights(3)
        for xie, we in zip(xi2pts, wxi2):
            fe2[i] += N(xie, i+1, xi2.enum)*we*2
    print "fe2", fe2, np.sum(fe2)

    F = np.zeros((4))
    F[0] = fe1[0]
    F[1] = fe1[1] + fe2[0]
    F[2] = fe1[2] + fe2[1]
    F[3] = fe2[2]
    print F
