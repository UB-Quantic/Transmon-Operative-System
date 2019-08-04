# SingleSuperconductingQubitExperiment CLASS

This is a class that façades the translation of generic algorithmic orders
to the actual analogical transmon like superconducting qubit implementation.

The class employs the Labber API to communicate with the specific hardware.

The idea is that one day this class could be the way in which quantum programmers
will remotely control the sate of the qubit and recieve measurement
results. Just by producing the execution of this class in a local directory of 
the quantum computer's controlling server and without requiring them to know 
nothing about the particular superconducitng
qubit implementation, full control of the qubit drive will be possible. 
Programmers will just need to work in the abstract Bloch Sphere scheme.
The class then makes all the hardware communications and both logical and analogical
compilations required to translate a desired unitary operaiton into an effective
experimental realizaion.

Note that the class depends on some files that should be in the same directory
as the defining script. These files, descibed bellow in the section "Requirements
for correct function", are supossed to be generated/updated by the administrator of
the quantum computer whenever the qubit remote control server is initialized.
They contain the required transmon empyrical parameters and hardware model
variables necessary to execute all the procedure. They are an attempt to avoid
hard coding and abstract from the main code the most arbitrary/variable attributes
of the system. 

The script defining the class is fully documented, so feel free to open it
and read the headers of each method.

## Main Methods

### driveQubitTo(self, wishedAlpha, wishedBeta, initialAlpha=1, initialBeta=0)
-----------------------------------------------------------------
Given a desired final state vector, in three steps, this method
realizes the physical single qubit rotation using the virtual Z gate
strategy that follows: 
    
a) logical compilation:
 - generation of the required single rotation matrix, 
 - translation of this rotation to a series of two pi/2 X gates 
 alternated in between three Z gates with certain angles
 
b) analogical compilation:
 - using the fact that Z gates that are followed by single qubit gates
 can be implemented in superconducting qubits in the capacitive coupling
 as an additional relative phase applied to the coming driving pulses;
 the sequence of Z and X pi/2 gates is translated to two pi/2 pulses
 that have the additional necessary phases. So, the rotation is
 translated into a pair of microwave pulses
 
c) physical execution of the experiment
 - Communication is established with the Instrument Server of Labber, 
 and channels to the particular devices are enabled
 - Pulses are sent to the driving devices, and they are armed to wait
 until the trigger
 - Logging channels are armed as well to allow inmediate readout
 - The trigger is sent and the interaction with the physical qubit begins
 - The readout signal is anlayzed and the method outputs the relevant
 information
 
     #### Arguments
     
     wishedAlpha : double (Complex)
         The projection on the |0> basis state we wish the qubit to drive to
     wishedBeta : double (Complex)
         The projection on the |1> basis state we wish the state to drive to
    
     #### Outputs
     
     The measurement's result
 
### executeUnitaryOperation(self,operation)
-------------------------------------------
It is the same method as executeRotation, but here one can directly
input a desired operation on the initial ground state |0> of the qubit,
instead of having to give in the output state.

    #### Argument

    operation : numpyArray, shape(2,2) Complex
        The unitary operation one wishes to apply to the qubit, written as
        a matrix in the computational basis
    
    #### Outputs
    
    The measurement's result 

## Requirements for Correct funciton

The class takes for granted that the following files are in the same directory
as its defining script, for the reasons exposed in the introduction:

- CircuitParams.txt
    A txt file containing info about the empyrical variables of the qubit
    It should have the following layout:    

        QubitFrequency <float>
        Anharmonicity <flaot>
        
- DeviceParams.txt
    A txt file containing information about the employed particular devices.
    
        AWG_SampleRate <flaot>
        wSSB <flaot>

- State.py
    This is a small side-class used by the main one, it is exposed bellow. 

## TODOs

In the header of the class' .py file there is a list of genral TODO-s in order
to expand the usage  of the class.

## Example

A trivial example of usage of this class: The preparation of a |+> State and 
inmediate measurement of the sate:
        
    myQubit = SingleSuperconductingQubitExperiment()
    result = myQubit.driveQubitTo( wishedAlpha=1/np.sqrt(2), wishedBeta=1/np.sqrt(2))
    print('The qubit-s state is',result)

## References used

[1] - P. Krantz, M. Kjaergaard, F. Yan, T. P. Orlando, S. Gustavsson, and
    W. D. Oliver, A quantum engineer’s guide to superconducting qubits, 
    Appl. Phys. Rev. 6, 021318 (2019).
    
[2] - McKay, D. C., Wood, C. J., Sheldon, S., Chow, J. M. & Gambetta, J. M. 
Efficient Z-gates for quantum computing. Preprint at 
https://arxiv.org/abs/1612.00858 (2017)
    
[3] - Michael A. Nielsen , Isaac L. Chuang, Quantum Computation and Quantum 
Information, Cambridge University Press, 2001

## Contact

Xabier Oyanguren Asua 

E-mail: <oiangu9@gmail.com>

# State side-class

This is an object to sustain the different representations of an state vector
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
    
    #### Parameters
    alpha : Complex double
        The projection of the state on |0>
    beta : Complex double
        The projection of the state on |1>
    
    Note the assertion that |alpha|**2 + |beta|**2 = 1
    
    #### Attributes
    
    _amplitudes : NumpyArray, shape(2,1) Complex
        stores the stateVector in the Computational basis
    _denisty_matrix : Numpyarray, shape(2,2) Complex
        stores the density matrix of the system in the computational basis. It
        is obtained in the initialization by denisty_matrix = |State> <State|
    _blochVector : NumpyArray, shape(1,3) Real
        Unit vector in the surface of a R3 sphere representing the State.
        It is obtained in the initialization using the denisty matrix as
        suggested in the first reference
        
    #### Reference
    
    https://en.wikipedia.org/wiki/Bloch_sphere
    
    Michael A. Nielsen , Isaac L. Chuang, Quantum Computation and Quantum 
        Information, Cambridge University Press, 2001
    
    
