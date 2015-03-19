#!/usr/bin/env python
"""
Created on Fri Jun 28 16:43:16 2013

Author: Oren Freifeld
Email: freifeld@dam.brown.edu
"""

import numpy as np
multiply = np.dot
_lsqrt = np.sqrt
_lcos = np.cos
_lsin = np.sin
#from numpy.linalg.linalg import norm



# If there are problems with the imports or functions from the cy.<module_name>
# modules, you may need to step  into these directories, and run
# "make clean" followed by "make". 



#import Cython
#cython_version = Cython.__version__

from cy.rodrigues import SO3so3
from cy.rot_btwn_two_vecs.rot_btwn_two_vecs import rmat_btwn_two_vecs



my_dtype = np.float32

class Aux:
    x_unit_vec = np.array([1,0,0],np.float)
    edges_tmp = np.zeros((3,2),my_dtype)
    _vec_a = np.zeros(3,my_dtype)
    _vec_b = np.zeros(3,my_dtype)
    _rotation_axis= np.zeros(3,my_dtype)
    _lperp= np.zeros(3,my_dtype)
    _rvec =  np.zeros(3,my_dtype)
    _rmat_step1 =  np.zeros((3,3),my_dtype)
    _rmat_step2 =  np.zeros((3,3),my_dtype)
    
    _AtimesS = np.zeros((3,3),my_dtype)

    
    @classmethod
    def rotation_matrix_between_two_vectors(cls,a,b):
        """ 
        Returns mat                    
        """            
        # Copy the data to local buffers
        cls._vec_a[:]=a
        cls._vec_b[:]=b    
        # Find a rotation rmat that sends a to b. 
        rmat_btwn_two_vecs(cls._vec_a,
                           cls._vec_b,
                           cls._rvec, 
                           cls._rmat_step1,
                           cls._rotation_axis,
                           cls._lperp) 
        return Aux._rmat_step1  
        

class Simple(object):
    _linear_system_A = np.zeros((2,2),my_dtype)
    _linear_system_b = np.zeros((2),my_dtype)
    
    @staticmethod
    def verify_is_3_by_2(X):       
        if X.shape != (3,2):
            err_msg = 'Expected shape = (3, 2); got {0} instead.'.format(X.shape)
            raise ValueError(err_msg)

    @staticmethod
    def verify_is_canonical(X): 
        Simple.verify_is_3_by_2(X)        
        if not np.allclose(X[-1],0,atol=1e-03):
            err_msg = 'Expected triangle to lie in the xy plane; got\n{0}\n instead.'.format(X)            
            raise ValueError(err_msg)
        if X[0,0] <= 0:
            err_msg = 'Expected 1st edge to point in the positive x direction; '
            err_msg += 'got\n{0}\n instead.'.format(X[:,0]) 
            raise ValueError(err_msg)
        if X[1,1] <= 0:
            err_msg = 'Expected 2nd edge to have a positive y coordinate; '
            err_msg += 'got\n{0}\n instead.'.format(X[:,1]) 
            raise ValueError(err_msg)            

    @classmethod 
    def compute_deformation(cls,X,Y,Q,A,S,Rx,Ry,tf_verify_checks=False):
        """
        Return Q.
        This assumes Rx and Ry have been pre-computed. 01/11/2014
        """
#        Rx[:] = Simple.compute_canonical_rotation(X)  
#        Ry[:] = Simple.compute_canonical_rotation(Y)    
        Simple.compute_canonical_rotation(X,out=Rx)  
        Simple.compute_canonical_rotation(Y,out=Ry)   

        if tf_verify_checks:
            Simple.verify_is_canonical(Rx.dot(X))
            Simple.verify_is_canonical(Ry.dot(Y))           
        
        AtimesS = Aux._AtimesS
   
        Simple.compute_planar_deformation(Rx.dot(X),Ry.dot(Y),A,S)
       
        if tf_verify_checks:
            Simple.verify_ASX_equals_Y(Rx.dot(X),Ry.dot(Y),A,S)
     
         # Technically, using eye*S should have 
         # worked too, as Rx time X should lie in the xy plane. 
         # However, due to accumlation of numerical errors
         # it is better to enforce 1 in the last entry. 
        AtimesS[:2,:2]=A[:2,:2]*S           
        AtimesS[2,2]=1       
                 
        Q[:]= Ry.T.dot(AtimesS).dot(Rx)
        return Q
            
    @classmethod
    def compute_scale(cls,X,Y,tf_verify_canonical=False):
        """
        Assumes X and Y are canonical.
        """
        if tf_verify_canonical:
            cls.verify_is_canonical(X)  
            cls.verify_is_canonical(Y) 
        
        S = Y[0,0]/X[0,0]
        return S
    
     
    @classmethod    
    def compute_planar_deformation(cls,X,Y,A,S,tf_verify_canonical=False):
        """
        Return A and S.
        After scaling is taken into account,
        the first column must be [1 0]^T since we want [1 0]^T to be 
        an eigenvector with 1 as an eigenvalue. 
        We have  
                    [1 a  [x1  = x2                         
                     0 b]  y1]   y2
        with a and b unknowns. We know that y1>0 and y2 >0.
        ==>    b = y2/y1 and a = (x2-x1)/y1  

         
        """  

        if tf_verify_canonical:
            cls.verify_is_canonical(X)  
            cls.verify_is_canonical(Y)
        
        S[0] = Simple.compute_scale(X,Y)

        x1,y1 =  S[0] * X[:2,1]
        x2,y2 =  Y[:2,1]
 
        a = (x2-x1)/y1
        b = y2 / y1
  
        A[:] = np.eye(3)
                     
        A[:2,1]=a,b
        A[2,2]=0
        return A,S

    @staticmethod
    def verify_ASX_equals_Y(X,Y,A,S):
        if not np.allclose(A.dot(S * X) ,  Y,atol=1e-03):
            raise ValueError('ASX =/= Y\n {0}\n\n{1}'.format(A.dot(S * X) ,Y))
            
    @staticmethod
    def verify_QX_equals_Y(X,Y,Q):
        if not np.allclose(Q.dot(X) ,  Y, atol=1e-03):
            raise ValueError('QX =/= Y')
    
    
    @classmethod
    def compute_canonical_rotations(cls,arr_of_edges,out):
        if arr_of_edges.shape[0]!=out.shape[0]:
            raise ValueError
        if arr_of_edges.shape[1:]!=(3,2):
            raise ValueError
        if out.shape[1:]!=(3,3):
            raise ValueError           
        
        for X,R in zip(arr_of_edges,out):
#            R[:]=cls.compute_canonical_rotation(X)
            cls.compute_canonical_rotation(X,out=R)
    @classmethod
    def compute_canonical_rotation(cls,edges,out):
        """
        Returns: rmat
        
        Remark:         
            There is another way to do it (using cross product)
            which is simpler code and understand. 
            However, that other method turned out to be slower (at least in 
            python). I suspect that in cython or Julia the other method would 
            in fact be faster. But this is still in the TODO list. 
            
            Oren Freifeld, 01/12/2014 
        
        """
        if edges.dtype != my_dtype:
            raise ValueError('Expected type {0}; got {1} instead.'.format(my_dtype,edges.dtype))
        if edges.shape != (3,2):
            raise ValueError('Expected (3,2); got {0} instead.'.format(edges.shape))
        
        first_edge,second_edge = edges[:,0],edges[:,1]                    
        rmat_step1 = Aux.rotation_matrix_between_two_vectors(first_edge,Aux.x_unit_vec)
        
        multiply(rmat_step1,edges,Aux.edges_tmp)                                             
        
        
        vec = Aux.edges_tmp[:,1]
        """
        [1 0 0     [ x1         [x2 = x1]                   # 
         0 c -s  *   y1    =    [y2 = +sqrt(y1^2+z1^2)      # isometry & y2>0
         0 s c]      z1]         z2 = 0                     #
         
        Unknowns: c and s 
        """
         
        A = cls._linear_system_A
        b = cls._linear_system_b
        
        A[0,0]=vec[1]
        A[0,1]=-vec[2]
        A[1,0]=vec[2]
        A[1,1]=vec[1]

        b[0]= _lsqrt(vec[1]**2+vec[2]**2)       
        detA = A[0,0]* A[1,1]-A[1,0]*A[0,1]        
        cos_theta =  A[1,1]*b[0] /detA
        sin_theta = -A[1,0]*b[0] /detA

        
        # "Proof": sin^2(theta)+cos^2=1
        if not np.abs(cos_theta **2 + sin_theta ** 2 - 1) < 1e-03:
            raise ValueError(cos_theta **2 + sin_theta ** 2)

        theta = np.arctan2(sin_theta,cos_theta)
        c = _lcos(theta)
        s = _lsin(theta)
        
        rmat_step2 = Aux._rmat_step2
        rmat_step2[:] = ([[1, 0 ,0 ],
                          [0,+c,-s ],
                          [0,+s,+c ]])                    
        rmat_step2.dot(rmat_step1,out=out)
#        return multiply(rmat_step2,rmat_step1) 
        
    @classmethod
    def compute_deformations(cls,Xs,Ys,Rxs,Rys,As,Ss,Qs,verbose=False):    
        for i in xrange(nTriangles):
            if verbose==True:
                if i % 1000 == 0:
                    print i
            # input
            X = Xs[i]
            Y = Ys[i]
            
            # Canonical rotations
            Rx = Rxs[i]
            Ry = Rys[i]        
            
            # output
            Q = Qs[i]
            A = As[i]
            S = Ss[i]                              
             
            Simple.compute_deformation(X,Y,Q,A,S,Rx,Ry,tf_verify_checks=0)


                

            
if __name__ == '__main__':


    edge_pair_shape = [3,2]
    nTriangles = 10000    
    
    # Edge-pairs for the template
    Xs = np.zeros((nTriangles,3,2),my_dtype)
    # Canonical rotations for the template
    Rxs = np.zeros((nTriangles,3,3),my_dtype)
    # Edge-pairs for the other mesh
    Ys = np.zeros_like(Xs)       
    # Canonical rotations for the other mesh
    Rys = np.zeros((nTriangles,3,3),my_dtype)

    # The 3x3 mesh deformations.     
    Qs = np.zeros((nTriangles,3,3),my_dtype)
    
    # The planar deformations
    As = np.zeros_like(Qs)
    As[:,0,0]= As[:,1,1]=As[:,1,1]=1 # initialization
    # Scales
    Ss = np.ones((nTriangles,1),dtype=my_dtype)
    
    # fake some data HERE YOU SHOULD USE YOUR OWN DATA = BUT NOTE THE DIMENSIONS.
    Xs[:] = np.random.standard_normal(Xs.shape)
    Ys[:] = np.random.standard_normal(Ys.shape)

    
    Simple.compute_canonical_rotations(Xs,out=Rxs)
    Simple.compute_canonical_rotations(Ys,out=Rys)
    

    print "Compute deformations:"
    Simple.compute_deformations(Xs,Ys,Rxs,Rys,As,Ss,Qs)      
    print 'Done.' 
    
    print 
    
    # Verify
    for i in xrange(nTriangles):
        if i % 5000 == 0 or i == nTriangles-1:
            print i
        # input
        X = Xs[i]
        Y = Ys[i]                     
        # output
        Q = Qs[i]
                         
        Simple.verify_QX_equals_Y(X,Y,Q)    
        
    print 'Verification done.'   
    
        
    
    
    