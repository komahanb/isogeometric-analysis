#!/usr/bin/env python
import numpy as np
import math
import numpy as np
np.set_printoptions(precision=4 , suppress=True)

from basis import N, Nprime
from basis import H, Hprime

import matplotlib.pyplot as plt

def assemble_stiffness(num_nodes, num_disps, elem_nodes, kset, bc_nodes):
    """
    Assemble and return the global stiffness matrix
    """
    K = np.zeros((num_nodes*num_disps,num_nodes*num_disps))

    # Assmble matrix
    for key in kset.keys():
        conn = elem_nodes[key]
        K[conn[0]*num_disps:conn[1]*num_disps+1, conn[0]*num_disps:conn[1]*num_disps+1] += kset[key]        

    # apply boundary conditions (wont work for two vars)
    K = np.delete(K, bc_nodes, 0)
    K = np.delete(K, bc_nodes, 1)
    
    return K

def assemble_force(num_nodes, num_disps, elem_nodes, fset, bc_nodes):
    F = np.zeros((num_nodes*num_disps))
    for key in fset.keys():
        conn = elem_nodes[key]
        F[conn[0]*num_disps:conn[1]*num_disps+1] += fset[key]
    return np.delete(F, bc_nodes, 0)

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
        return str(self.__dict__)
    
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

    """
    """
    a = 0.0
    b = 2.0
    
    num_nodes      = 4
    num_elems      = num_nodes - 1
    num_elem_nodes = 2
    num_disps      = 1
    ndofs          = num_nodes*num_disps
    le             = (b-a)/num_elems
    xvals = np.linspace(a,b,num_nodes)
    nodes = []
    for n in range(num_nodes):
        nodes.append(n)
    print "nodes", nodes

    bc_nodes = [nodes[0]]
    
    elem_nodes = {}
    for e in range(num_elems):        
        elem_nodes[e] = [nodes[e], nodes[e+1]]
    print elem_nodes
    
    xi = {}
    for e in range(num_elems):
        xi[e] = KnotParameter(e, [le*e, le*(e+1)])
    
    N = num_elem_nodes*num_disps
    
    #  
    det = 1

    # Gather local element stiffness matrices
    K = {}
    for key in xi.keys():
        xipts, wpts = xi[key].get_quadrature_points_weights(3)
        le = xi[key].bounds[1]-xi[key].bounds[0]
        K[key] = np.zeros((N,N))
        for i in range(0,N):
            for j in range(0,N):
                for xinode, wnode in zip(xipts, wpts):
                    K[key][i,j] += Hprime(xinode, le, i+1)*Hprime(xinode, le, j+1)*wnode/det

    # Gather local element force vectors
    F = {}
    for key in xi.keys():
        xipts, wpts = xi[key].get_quadrature_points_weights(3)
        le = xi[key].bounds[1]-xi[key].bounds[0]
        F[key] = np.zeros((N))
        for i in range(0,N):
            for xinode, wnode in zip(xipts, wpts):
                F[key][i] += H(xinode, le, i+1)*wnode*det
    
    kmat = assemble_stiffness(num_nodes, num_disps, elem_nodes, K, bc_nodes)
    fmat = assemble_force(num_nodes, num_disps, elem_nodes, F, bc_nodes)

    print kmat
    print fmat
    
    # Solve the linear system
    u = np.linalg.solve(kmat,fmat)
    
    print u

    plt.figure()
    plt.plot(xvals[1:], u , '-o', label='')
    plt.show()

    #TODO
    ## 1. Find det of jacobian
    ## 2. Use quadratic elemnet
    ## 3. Use same setup for isogeometric analysis
