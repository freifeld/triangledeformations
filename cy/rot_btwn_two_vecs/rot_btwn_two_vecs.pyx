 
import numpy as np
cimport numpy as np

from cpython cimport bool
 
    

cdef extern from "math.h":
    float sqrt(float) 
    float asin(float)
    float acos(float)
    float sin(float)  
    float cos(float) 
    
    
    

def rvec_btwn_two_vecs(np.ndarray[np.float32_t,ndim=1] a not None,
                       np.ndarray[np.float32_t,ndim=1] b not None,
                       np.ndarray[np.float32_t,ndim=1] rvec not None,
                       np.ndarray[np.float32_t,ndim=1] rotation_axis not None,
                       np.ndarray[np.float32_t,ndim=1] perp not None):  
    """
    Input: a,b 
    Output: rvec
    Note: a and b will be modified
    """
  

    a /= sqrt(a[0]*a[0]+a[1]*a[1] +a[2]*a[2])
    b /= sqrt(b[0]*b[0]+b[1]*b[1] +b[2]*b[2])

#    rotation_angle = np.arccos( a[0] * b[0] + a[1] * b[1] + a[2] * b[2] )
    cdef:
        float rotation_angle
        float tmp
        float pi = 3.141592653589793
        
        float a0 = a[0], a1 = a[1], a2 = a[2]
        float b0 = b[0], b1 = b[1], b2 = b[2]
        bool tf_sanity_checks = False
     
    rotation_angle = acos( a0 * b0 + a1 * b1 + a2 * b2 )
            
    rotation_axis[0] = a1*b2-b1*a2
    rotation_axis[1] = -a0*b2+b0*a2
    rotation_axis[2] = a0*b1-b0*a1
    

    tmp = sqrt(rotation_axis[0]*rotation_axis[0]+
               rotation_axis[1]*rotation_axis[1]+
               rotation_axis[2]*rotation_axis[2])
               
    if tmp == 0:
        if (a * b).sum() == 1:
            # These vectors are equal. No need to rotate.
            rvec[:]=0
            return rvec
        if np.allclose(a , -b):
            # HACK. TODO: FIX THIS.
            if  np.allclose(a , np.array([1.0,0,0])) or np.allclose(a , np.array([-1.0,0,0])):
                rvec[:] = [ 0.        ,  0.        ,  pi]
                return rvec                
            raise NotImplementedError
        # Let's find some vector perpendicular to a.

        if a0:
            perp[0]= -1.0 * a0
            perp[1]= +1.0 * a1                
        elif a1:               
            perp[1]= +1.0 * a1  
        elif a2:                              
            perp[:]= (1.0,0,0)
        else:
            raise NotImplementedError
            
        if  tf_sanity_checks:
            if not np.allclose(np.inner(a , perp) , 0):
                print 'a and b'
                print a,b
                print 'np.inner(a , perp)'
                print np.inner(a , perp)
                raise ValueError('Failed to find perpendicular vector')
                
            if (perp[0] * perp[0] + perp[1] * perp[1] + perp[2] * perp[2] ) == 0:
                raise ValueError
                
        tmp = sqrt((perp ** 2).sum())
        rotation_axis[:] = perp
      
    rotation_axis /= tmp

    rvec[0] = rotation_angle * rotation_axis[0]
    rvec[1] = rotation_angle * rotation_axis[1]
    rvec[2] = rotation_angle * rotation_axis[2]
    
    return rvec     
    
    
    
    
def rmat_btwn_two_vecs(np.ndarray[np.float32_t,ndim=1] a not None,
                       np.ndarray[np.float32_t,ndim=1] b not None,
                       np.ndarray[np.float32_t,ndim=1] rvec not None,
                       np.ndarray[np.float32_t,ndim=2] rmat not None,
                       np.ndarray[np.float32_t,ndim=1] rotation_axis not None,
                       np.ndarray[np.float32_t,ndim=1] perp not None):
    """
    Input: a,b 
    Output: rmat
    Note: a and b will be modified
    """                           
    rvec_btwn_two_vecs(a,b,rvec,rotation_axis,perp)                                                      
    rvec2rmat_float32(rvec,rmat)                       
    return rmat                         
    
                           
def rvec2rmat_float32(np.ndarray[np.float32_t,ndim=1] rvec not None,
                       np.ndarray[np.float32_t,ndim=2] rmat not None):
    """
        rvec2rmat_float32(rvec,rmat)
        rvec.shape should be (3,)
        rmat.shape should be (3,3)        
        Operation is done inplace (in rmat).
    """                    
    
    cdef:
        float t  # Theta
        float s  # Sin
        float c  # Cos
        
        float r0
        float r1
        float r2 
               
        # r is going to be the normalized version of rvec
#        np.ndarray r = np.zeros(3, dtype= np.float32)  
    
    t  =  sqrt(rvec[0] * rvec[0] +
                  rvec[1] * rvec[1] + 
                  rvec[2] * rvec[2])
    if t:              
        r0 = rvec[0]/t 
        r1 = rvec[1]/t 
        r2 = rvec[2]/t        
        
        s = sin(t)
        c = cos(t)

#        # diagonal    
        rmat[0,0] = c + (1-c)*r0*r0 
        rmat[1,1] = c + (1-c)*r1*r1 
        rmat[2,2] = c + (1-c)*r2*r2
        
        # the rest
        rmat[0,1] = (1-c)*r0*r1-s * r2 
        rmat[0,2] = (1-c)*r0*r2+s * r1
        rmat[1,2] = (1-c)*r1*r2-s * r0 
        
        rmat[1,0] = (1-c)*r0*r1+s * r2 
        rmat[2,0] = (1-c)*r0*r2-s * r1  
        rmat[2,1] = (1-c)*r1*r2+s * r0

    else:
        # diagonal    
        rmat[0,0] = 1 
        rmat[1,1] = 1 
        rmat[2,2] = 1
        
        # the rest
        rmat[0,1] = 0
        rmat[0,2] = 0
        rmat[1,2] = 0 
        
        rmat[1,0] = 0  
        rmat[2,0] = 0  
        rmat[2,1] = 0
                                       
                           
                           
                           
