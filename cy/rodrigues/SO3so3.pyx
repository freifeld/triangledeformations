 
import numpy as np
cimport numpy as np


import  cv2    
#def rmat2rvec_1(np.ndarray[np.float32_t,ndim=2] rmat not None,
#                       np.ndarray[np.float32_t,ndim=1] rvec not None):            
#    rvec[:]=cv2.Rodrigues(src=rmat,dst=rvec)[0].flatten()
    

cdef extern from "math.h":
    float sqrt(float) 
    float asin(float)
    float acos(float)
    float sin(float)  
    float cos(float) 
 
def rmat2rvec_float32(np.ndarray[np.float32_t,ndim=2] rmat not None,
                       np.ndarray[np.float32_t,ndim=1] rvec not None):            
    """
        rmat2rvec_float32(rmat,rvec)
        rmats.shape should be 3,3
        rvecs.shape should be 3
        Operation is done inplace (in rvec).
    """
     
    
    cdef: 
        float theta
        float traceR_minus1_over2   #   (trace(R)-1) / 2         
    
    traceR_minus1_over2 = ((rmat[0,0] + rmat[1,1] + rmat[2,2]) -1 ) / 2
    
    # Numerics might make C-function acos unhappy. Make sure we are in [-1,+1]
    if traceR_minus1_over2 <= -1:
#        print 'A'
        traceR_minus1_over2 = -1
        # arcos(-1) is -pi so the sin will be zero. 
        # We can't divide by the sin.
        # I have no time to fix it now. So cowardly resort to opencv 
        cv2.Rodrigues(src=rmat,dst=rvec)
        return
    elif traceR_minus1_over2 >= 1:
#        print 'B'
        traceR_minus1_over2 = 1
        # acos(+1) is pi so the sin will be zero. 
        # We can't divide by the sin.
        cv2.Rodrigues(src=rmat,dst=rvec) 
        return
#    else:
#        print 'C' ,  traceR_minus1_over2
    theta = acos(traceR_minus1_over2)
    print 'theta' ,  theta
    if theta:
        rvec[0] = -(rmat[1,2]-rmat[2,1])/2     
        rvec[1] =  (rmat[0,2]-rmat[2,0])/2
        rvec[2] = -(rmat[0,1]-rmat[1,0])/2   
        rvec *= theta / sin(theta)
    else:
        rvec[:]=0 


 
        
#def rmats2rvecs_1_float32(np.ndarray[np.float32_t,ndim=3] rmats not None,
#                         np.ndarray[np.float32_t,ndim=2] rvecs not None):
#    for rvec,rmat in zip(rvecs,rmats):
#        rvec[:]=cv2.Rodrigues(src=rmat,dst=rvec)[0].flatten()

def rmats2rvecs_2_float32(np.ndarray[np.float32_t,ndim=3] rmats not None,
                         np.ndarray[np.float32_t,ndim=2] rvecs not None):
    for rvec,rmat in zip(rvecs,rmats):
        rmat2rvec_float32(rmat,rvec)


def rmats2rvecs_3_float32(np.ndarray[np.float32_t,ndim=3] rmats not None,
                         np.ndarray[np.float32_t,ndim=2] rvecs not None):
    cdef:
        float sinTheta 
        np.ndarray rvec = np.zeros(3, dtype= np.float32)
        np.ndarray rmat = np.zeros((3,3), dtype= np.float32)
    for rvec,rmat in zip(rvecs,rmats):
        rvec[0] = -(rmat[1,2]-rmat[2,1])/2     
        rvec[1] =  (rmat[0,2]-rmat[2,0])/2
        rvec[2] = -(rmat[0,1]-rmat[1,0])/2
           
        sinTheta = sqrt(rvec[0] * rvec[0] +
                        rvec[1] * rvec[1] + 
                        rvec[2] * rvec[2])
        
        if sinTheta:
            rvec  /= sinTheta
            rvec  *= asin(sinTheta)       
        else:
            rvec[:]=0                 

def rmats2rvecs_float32(np.ndarray[np.float32_t,ndim=3] rmats not None,
                         np.ndarray[np.float32_t,ndim=2] rvecs not None):
    """
        rmats2rvecs_float32(rmats,rvecs)
        rmats.shape should be N,3,3
        rvecs.shape should be N,3
        Operation is done inplace (in rvecs).
    """
 
    
#    cdef np.ndarray rvec = np.zeros(3, dtype= np.float32)
#    cdef np.ndarray rmat = np.zeros((3,3), dtype= np.float32)
    
    cdef:
        Py_ssize_t i
        float r0
        float r1
        float r2

   
        float theta
        float traceR_minus1_over2   #   (trace(R)-1) / 2   
        float theta_over_sinTheta        
        
    for i in range(rvecs.shape[0]):
      

        traceR_minus1_over2 = ((rmats[i,0,0] + rmats[i,1,1] + rmats[i,2,2]) -1 ) / 2
        
        # Numerics might make C-function acos unhappy. Make sure we are in [-1,+1]
        if traceR_minus1_over2 <= -1:
            traceR_minus1_over2 = -1
            # arcos(-1) is -pi so the sin will be zero. 
            # We can't divide by the sin.
            # I have no time to fix it now. So cowardly resort to opencv 
            cv2.Rodrigues(src=rmats[i],dst=rvecs[i])
            continue
        elif traceR_minus1_over2 >= 1:
            traceR_minus1_over2 = 1
            # acos(+1) is pi so the sin will be zero. 
            # We can't divide by the sin.
            cv2.Rodrigues(src=rmats[i],dst=rvecs[i]) 
            continue
        theta = acos(traceR_minus1_over2)

     
        if theta:
            rvecs[i,0] = -(rmats[i,1,2]-rmats[i,2,1])/2      
            rvecs[i,1] =  (rmats[i,0,2]-rmats[i,2,0])/2
            rvecs[i,2] = -(rmats[i,0,1]-rmats[i,1,0])/2 
            
            theta_over_sinTheta = theta / sin(theta)
            rvecs[i,0] *= theta_over_sinTheta
            rvecs[i,1] *= theta_over_sinTheta
            rvecs[i,2] *= theta_over_sinTheta
        else: # probably wrong!
            rvecs[i,0]=0
            rvecs[i,1]=0
            rvecs[i,2]=0 

#    print 'YES'
         
        
        
        
 ##########################
#def rvec2rmat_1_float32(np.ndarray[np.float32_t,ndim=1] rvec not None,
#                       np.ndarray[np.float32_t,ndim=2] rmat not None):
#    # cv2 has a bug: dst is not affected. Workaround: use LHS.                       
#    rmat[:]=cv2.Rodrigues(src=rvec,dst=rmat)[0]   

#r is going to be the normalized version of rvec
#cdef:

#    np.ndarray r = np.zeros(3, dtype= np.float32)     
#    list r = [0.0,0.0,0.0]
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
#        r[0] = rvec[0]/t 
#        r[1] = rvec[1]/t 
#        r[2] = rvec[2]/t
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
        
#        # diagonal    
#        rmat[0,0] = c + (1-c)*r[0]*r[0] 
#        rmat[1,1] = c + (1-c)*r[1]*r[1] 
#        rmat[2,2] = c + (1-c)*r[2]*r[2]
#        
#        # the rest
#        rmat[0,1] = (1-c)*r[0]*r[1]-s * r[2] 
#        rmat[0,2] = (1-c)*r[0]*r[2]+s * r[1]
#        rmat[1,2] = (1-c)*r[1]*r[2]-s * r[0] 
#        
#        rmat[1,0] = (1-c)*r[0]*r[1]+s * r[2] 
#        rmat[2,0] = (1-c)*r[0]*r[2]-s * r[1]  
#        rmat[2,1] = (1-c)*r[1]*r[2]+s * r[0]

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
            
    
def rvecs2rmats_2_float32(np.ndarray[np.float32_t,ndim=2] rvecs not None,
                        np.ndarray[np.float32_t,ndim=3] rmats not None):
    cdef:
        float t  # Theta
        float s  # Sin
        float c  # Cos
        # r is going to be the normalized version of rvec
        # np.ndarray r = np.zeros(3, dtype= np.float32)     
        float r0
        float r1
        float r2

        Py_ssize_t i


        
    for i in range(rvecs.shape[0]):   
        
        t  =  sqrt(rvecs[i,0] * rvecs[i,0] +
                   rvecs[i,1] * rvecs[i,1] + 
                   rvecs[i,2] * rvecs[i,2])
        if t:
            s = sin(t)
            c = cos(t)
            
            r0 = rvecs[i,0]/t 
            r1 = rvecs[i,1]/t 
            r2 = rvecs[i,2]/t 
            
            # diagonal    
            rmats[i,0,0] = c + (1-c)*r0*r0 
            rmats[i,1,1] = c + (1-c)*r1*r1 
            rmats[i,2,2] = c + (1-c)*r2*r2
            
            # the rest
            rmats[i,0,1] = (1-c)*r0*r1-s * r2 
            rmats[i,0,2] = (1-c)*r0*r2+s * r1
            rmats[i,1,2] = (1-c)*r1*r2-s * r0 
            
            rmats[i,1,0] = (1-c)*r0*r1+s * r2 
            rmats[i,2,0] = (1-c)*r0*r2-s * r1  
            rmats[i,2,1] = (1-c)*r1*r2+s * r0 
            
#            r[0] = rvecs[i,0]/t 
#            r[1] = rvecs[i,1]/t 
#            r[2] = rvecs[i,2]/t 
#            
#            # diagonal    
#            rmats[i,0,0] = c + (1-c)*r[0]*r[0] 
#            rmats[i,1,1] = c + (1-c)*r[1]*r[1] 
#            rmats[i,2,2] = c + (1-c)*r[2]*r[2]
#            
#            # the rest
#            rmats[i,0,1] = (1-c)*r[0]*r[1]-s * r[2] 
#            rmats[i,0,2] = (1-c)*r[0]*r[2]+s * r[1]
#            rmats[i,1,2] = (1-c)*r[1]*r[2]-s * r[0] 
#            
#            rmats[i,1,0] = (1-c)*r[0]*r[1]+s * r[2] 
#            rmats[i,2,0] = (1-c)*r[0]*r[2]-s * r[1]  
#            rmats[i,2,1] = (1-c)*r[1]*r[2]+s * r[0]   
        else:
            # diagonal    
            rmats[i,0,0] = 1
            rmats[i,1,1] = 1 
            rmats[i,2,2] = 1
            
            # the rest
            rmats[i,0,1] = 0 
            rmats[i,0,2] = 0
            rmats[i,1,2] = 0
            
            rmats[i,1,0] = 0
            rmats[i,2,0] = 0
            rmats[i,2,1] = 0               

def rvecs2rmats_float32(np.ndarray[np.float32_t,ndim=2] rvecs not None,
                        np.ndarray[np.float32_t,ndim=3] rmats not None ,
                        np.ndarray[np.float32_t,ndim=2] rmats_twoD not None ):
    """
        rvecs2rmats_float32(rmats,rvecs,rmats_twoD)
        rvecs.shape should be N,3
        rmats.shape should be N,3,3
        rmats_twoD.shape should be N,9
        Note it is assumed that rmats_twoD is in fact rmats.reshape((rmats.shape[0],9)).
        In effect, they both point to the same data. Of course, we could have done it
        inside cython, but (to my surprise), this clumsy way is in fact slightly faster.
        
        Operation is done inplace (in rmats_twoD, and this also in rmats).
    """
    cdef:
        float t  # Theta
        float s  # Sin
        float c  # Cos
        float one_minus_c # 1 - Cos
        # r0,r1,r2, is going to be the normalized version of rvec   
        float r0
        float r1
        float r2
        Py_ssize_t i 
     
#    cdef np.ndarray rmats_twoD =  rmats.reshape(rmats.shape[0],9)
#    rmats_twoD =  rmats.reshape(rmats.shape[0],9)
    
    cdef:
        float rvecsi0
        float rvecsi1
        float rvecsi2
    for i in range(rvecs.shape[0]):   
      
#        t  =  sqrt(rvecs[i,0] * rvecs[i,0] +
#                   rvecs[i,1] * rvecs[i,1] + 
#                   rvecs[i,2] * rvecs[i,2])
#      
#        r0 = rvecs[i,0]/t 
#        r1 = rvecs[i,1]/t 
#        r2 = rvecs[i,2]/t

        rvecsi0 = rvecs[i,0]
        rvecsi1 = rvecs[i,1]
        rvecsi2 = rvecs[i,2]
        
        t  =  sqrt(rvecsi0 * rvecsi0 +
                   rvecsi1 * rvecsi1 + 
                   rvecsi2 * rvecsi2)
        if t:
            r0 = rvecsi0/t 
            r1 = rvecsi1/t 
            r2 = rvecsi2/t
            
            
            s = sin(t)
            c = cos(t)
            
            one_minus_c = 1 - c
#            # diagonal    
#            rmats_twoD[i,0] = c + one_minus_c*r0*r0 
#            rmats_twoD[i,4] = c + one_minus_c*r1*r1 
#            rmats_twoD[i,8] = c + one_minus_c*r2*r2
#            
#            # the rest
#            rmats_twoD[i,1] = one_minus_c*r0*r1-s * r2 
#            rmats_twoD[i,2] = one_minus_c*r0*r2+s * r1
#            rmats_twoD[i,5] = one_minus_c*r1*r2-s * r0 
#            
#            rmats_twoD[i,3] = one_minus_c*r0*r1+s * r2 
#            rmats_twoD[i,6] = one_minus_c*r0*r2-s * r1  
#            rmats_twoD[i,7] = one_minus_c*r1*r2+s * r0


               
            rmats_twoD[i,0] = c + one_minus_c*r0*r0 
            rmats_twoD[i,1] = one_minus_c*r0*r1-s * r2
            rmats_twoD[i,2] = one_minus_c*r0*r2+s * r1
            rmats_twoD[i,3] = one_minus_c*r0*r1+s * r2
            rmats_twoD[i,4] = c + one_minus_c*r1*r1 
            rmats_twoD[i,5] = one_minus_c*r1*r2-s * r0 
            rmats_twoD[i,6] = one_minus_c*r0*r2-s * r1
            rmats_twoD[i,7] = one_minus_c*r1*r2+s * r0
            rmats_twoD[i,8] = c + one_minus_c*r2*r2 

        else:
            # diagonal    
            rmats_twoD[i,0] = 1 
            rmats_twoD[i,4] = 1 
            rmats_twoD[i,8] = 1
            
            # the rest
            rmats_twoD[i,1] = 0
            rmats_twoD[i,2] = 0
            rmats_twoD[i,5] = 0 
            
            rmats_twoD[i,3] = 0 
            rmats_twoD[i,6] = 0  
            rmats_twoD[i,7] = 0            

  



def rmats2rvecs_float64(np.ndarray[np.float64_t,ndim=3] rmats not None,
                         np.ndarray[np.float64_t,ndim=2] rvecs not None):
    """
        rmats2rvecs_float64(rmats,rvecs)
        rmats.shape should be N,3,3
        rvecs.shape should be N,3
        Operation is done inplace (in rvecs).
    """
    raise NotImplementedError('Need to fix a bug like in the 32bit')
    cdef float sinTheta 
    cdef float theta_over_sinTheta
 
    
    cdef:
        Py_ssize_t i
        float r0
        float r1
        float r2
        
    for i in range(rvecs.shape[0]): 
        r0 = -(rmats[i,1,2]-rmats[i,2,1])/2     
        r1 =  (rmats[i,0,2]-rmats[i,2,0])/2
        r2 = -(rmats[i,0,1]-rmats[i,1,0])/2           
        sinTheta = sqrt(r0 * r0 + r1 * r1 + r2 * r2)
       
        if sinTheta:

            # Now of course this SHOULD be in [0,1], but there are numerical issues...
            # As a result, asin will return NaN, and the game is over. Workaround:
            if sinTheta > 1.0:
                sinTheta = 1.0


            theta_over_sinTheta = asin(sinTheta) / sinTheta
             
 
            rvecs[i,0]  = r0 * theta_over_sinTheta 
            rvecs[i,1]  = r1 * theta_over_sinTheta
            rvecs[i,2]  = r2 * theta_over_sinTheta 
            
 

        else:
 
            rvecs[i,0]=0
            rvecs[i,1]=0
            rvecs[i,2]=0 

def rvecs2rmats_float64(np.ndarray[np.float64_t,ndim=2] rvecs not None,
                        np.ndarray[np.float64_t,ndim=3] rmats not None ,
                        np.ndarray[np.float64_t,ndim=2] rmats_twoD not None ):
    cdef:
        float t  # Theta
        float s  # Sin
        float c  # Cos
        float one_minus_c # 1 - Cos
        # r0,r1,r2, is going to be the normalized version of rvec   
        float r0
        float r1
        float r2
        Py_ssize_t i 
     
 
    
    for i in range(rvecs.shape[0]):   
       
        t  =  sqrt(rvecs[i,0] * rvecs[i,0] +
                   rvecs[i,1] * rvecs[i,1] + 
                   rvecs[i,2] * rvecs[i,2])
        if t:
            r0 = rvecs[i,0]/t 
            r1 = rvecs[i,1]/t 
            r2 = rvecs[i,2]/t
            
            
            s = sin(t)
            c = cos(t)
            
            one_minus_c = 1 - c
            # diagonal    
            rmats_twoD[i,0] = c + one_minus_c*r0*r0 
            rmats_twoD[i,4] = c + one_minus_c*r1*r1 
            rmats_twoD[i,8] = c + one_minus_c*r2*r2
            
            # the rest
            rmats_twoD[i,1] = one_minus_c*r0*r1-s * r2 
            rmats_twoD[i,2] = one_minus_c*r0*r2+s * r1
            rmats_twoD[i,5] = one_minus_c*r1*r2-s * r0 
            
            rmats_twoD[i,3] = one_minus_c*r0*r1+s * r2 
            rmats_twoD[i,6] = one_minus_c*r0*r2-s * r1  
            rmats_twoD[i,7] = one_minus_c*r1*r2+s * r0
        else:
            # diagonal    
            rmats_twoD[i,0] = 1 
            rmats_twoD[i,4] = 1
            rmats_twoD[i,8] = 1
            
            # the rest
            rmats_twoD[i,1] = 0
            rmats_twoD[i,2] = 0
            rmats_twoD[i,5] =0
            
            rmats_twoD[i,3] = 0
            rmats_twoD[i,6] = 0
            rmats_twoD[i,7] = 0      

def rvec2rmat_float64(np.ndarray[np.float64_t,ndim=1] rvec not None,
                       np.ndarray[np.float64_t,ndim=2] rmat not None): 
    
    cdef:
        float t  # Theta
        float s  # Sin
        float c  # Cos
        # r is going to be the normalized version of rvec
        np.ndarray r = np.zeros(3, dtype= np.float32)  
    
    t  =  sqrt(rvec[0] * rvec[0] +
                  rvec[1] * rvec[1] + 
                  rvec[2] * rvec[2])
    if t:
        r[0] = rvec[0]/t 
        r[1] = rvec[1]/t 
        r[2] = rvec[2]/t
        
        
        s = sin(t)
        c = cos(t)
        
        # diagonal    
        rmat[0,0] = c + (1-c)*r[0]*r[0] 
        rmat[1,1] = c + (1-c)*r[1]*r[1] 
        rmat[2,2] = c + (1-c)*r[2]*r[2]
        
        # the rest
        rmat[0,1] = (1-c)*r[0]*r[1]-s * r[2] 
        rmat[0,2] = (1-c)*r[0]*r[2]+s * r[1]
        rmat[1,2] = (1-c)*r[1]*r[2]-s * r[0] 
        
        rmat[1,0] = (1-c)*r[0]*r[1]+s * r[2] 
        rmat[2,0] = (1-c)*r[0]*r[2]-s * r[1]  
        rmat[2,1] = (1-c)*r[1]*r[2]+s * r[0]    
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