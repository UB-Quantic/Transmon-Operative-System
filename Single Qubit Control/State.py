# Author: Xabier Oyanguren Asua <oiangu9@gmail.com>
#TODOOOO! Introduce subclassing to implement mixed state manipulation

import numpy as np
import scipy.optimize as opt

class PureState:
    """An object to sustain the different representations of an state vector
    in the two dimensional Hilbert space (for pure states):
        
      - Rather as a superposition of the two states of the computational basis
    in the form |State> = a |0> + b |1> with a,b Complex scalars;
      - As a denisty matrix, given by:
          denisty_matrix = |State> <State|
      - Or as a unit vector in the surface of the Bloch Sphere: 
        |State> = cos(theta/2) |0> + e**(i*phi)*sin(theta/2) |1>
        for theta in (0,pi), phi in (0,2pi)
        thus, the point in the sphere is:
            (u,v,w) = (sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)) in R3
        
    Parameters
    ---------
    alpha : Complex double
        The projection of the state on |0>
    beta : Complex double
        The projection of the state on |1>
    
    Note the assertion that |alpha|**2 + |beta|**2 = 1
    
    Attributes
    ---------
    _amplitudes : NumpyArray, shape(2,1) Complex
        stores the stateVector in the Computational basis
    _denisty_matrix : Numpyarray, shape(2,2) Complex
        stores the density matrix of the system in the computational basis. It
        is obtained in the initialization by denisty_matrix = |State> <State|
    _blochVector : NumpyArray, shape(1,3) Real
        Unit vector in the surface of a R3 sphere representing the State.
        It is obtained in the initialization using the denisty matrix as
        suggested in the first reference
        
    Reference
    --------
    https://en.wikipedia.org/wiki/Bloch_sphere
    
    Michael A. Nielsen , Isaac L. Chuang, Quantum Computation and Quantum 
        Information, Cambridge University Press, 2001
    
    """
    
    def __init__(self,alpha, beta):
        assert (np.absolute(alpha)**2 + np.absoute(beta)**2 ==1), "WaveVector must be Normalized"
            
        self._amplitudes = np.array([[alpha],[beta]])
        self._density_matrix = np.dot(self._amplitudes,
                      self._amplitudes.transpose().conjugate())
        self._blochVector = np.array([2*self._density_matrix.real[0,1], 
                                      2*self._density_matrix.imag[1,0],
                                      self._density_matrix[0,0]-self._density_matrix[1,1]])
    
    def getAmplitudes(self):
        return self._amplitudes
    
    def getBlochVector(self):
        return self._blochVector
    
    def getDensityMatrix(self):
        return self._density_matrix
        
        