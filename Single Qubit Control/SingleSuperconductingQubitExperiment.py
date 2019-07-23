# Author: Xabier Oyanguren Asua <oiangu9@gmail.com>

import numpy as np
from State import State
from HZbasis import HZbasis
import scipy.optimize as opt
import Labber as Labber

class SingleSuperconductingQubitExperiment:
    """ An object that façades the translation of generic algorithmic orders
    to the actual analogical transmon like superconducting qubit implementation. 
    
    Parameters
    ---------
    drag : Boolean , default=True
        If True, DRAG pulsing will be used during capacitive coupling. Note
        that in this case, a drag beta callibration (to reduce leakage) will 
        be required, which might take an undesired time
    
    Attributes
    ---------
    circuitParameters : python dictionary
    
    deviceSpecs : python dictionary
    
    Xpi_2Pulse : python dictionary
    
    dragBeta : real double
    
    
    NOTE
    ----
        Throughout the class, it will always be assumed that operator matrixes
    are expressed in the computational basis [|0>,|1>] of the quantization axis
    depicted by the z axis in the Bloch Sphere
    
    TODO
    ----
        Implement dynamical pulsing to block decoherenig effects by longitudinal
    noise. DragZ pulsing is an alternative to this.
        Implement a Qubit simulator to see the theoretical effects on the Bloch
    sphere due to the sent pulses (Labber apparently has a module for this).
        Support for mixed states should be implemented, as well as a more general
    class for an N qubit system (for now, the Lab only has a single qubit, so
    this class covers all the possible operations on it).
        A calibration for the pulses should be implemented, as the pulse that
    leaves out the devices is not really the pulse arriving to the qubit!
        DRAGZ should be imnplemented in order to get both Stark shift and
    leakage correcitons working at once
        Perhaps another way of introducing the experimental parameters (that's
    capacitance values, impedances etc) for the correct pulse genertaion could
    be thought. The point is that the only pulse shapes used for driving will
    be X pi/2.
        Should there be an initial calibration funciton for getting the qubit
    resonace freqcy, T1,T2... etc, or rather this will be done already externally
    and these parameters will be updated in the txt file from where the class
    will read them? Ideally, if this is to be used by programers remotely, 
    maybe in each program (for each created singlequbit instance) a full 
    calibration should be performed?

    NOTES
    -----
    #SHOULD READOUT BE DONE SEPARATEDLY OF THE COUPLING? making them in the same
    function is required for the arm/trig way. If instead coupling and readout 
    are done in different Labber config files, perhaps dividing drive and readout
    in two methods is more convenient, such that the external user can directly
    choose the readout technique used (or not used) directly when called
    
    BIBLIOGRAPHY
    (for the Qubit Control scheme by Capacitive Coupling and Virtual Z gating)
    ---------
    [1] - P. Krantz, M. Kjaergaard, F. Yan, T. P. Orlando, S. Gustavsson, and
    W. D. Oliver, A quantum engineer’s guide to superconducting qubits, 
    Appl. Phys. Rev. 6, 021318 (2019).
    
    [2] - McKay, D. C., Wood, C. J., Sheldon, S., Chow, J. M. & Gambetta, J. M. 
    Efficient Z-gates for quantum computing. Preprint at 
    https://arxiv.org/abs/1612.00858 (2017)
    
    [3] - Michael A. Nielsen , Isaac L. Chuang, Quantum Computation and Quantum 
    Information, Cambridge University Press, 2001
    
    """
    
    
    def __init__(self, drag=True, simulation=False, dragz=False, dynamical=False):

        self.drag = drag
        self.simulation = simulation
        
        # at instance initialization, the parameters of the transmon should be
        # extracted from the local file generated in the calibration of the
        # system when the server was initialized by the QC admin
        self.circuitParameters = self._extractParamsFromFile("CircuitParams.txt")
        self.deviceSpecs = self._extractParamsFromFile("DeviceParams.txt")
        
        self.Xpi_2pulse = self._calibrateXPulse() #dictionary with T...
        if (drag):
            self.dragBeta = self._calibrateDragBeta()
        
    def driveQubitTo(self, wishedAlpha, wishedBeta, initialAlpha=0, initialBeta=0):
        """Given a desired final state vector, in three steps, this method
        realizes the physical singlequbit rotation using the virtual Z gate
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
         
         Arguments
         ---------
         wishedAlpha : double (Complex)
             The projection on the |0> basis state we wish the qubit to drive to
         wishedBeta : double (Complex)
             The projection on the |1> basis state we wish the state to drive to
        
         Outputs
         -------
         The measurement's result
        
        """
        # should it be established a timer to avoid calling this method too 
        # fast? I mean, in order to leave the qubit enough time to return to
        # the ground state
        ## while(timer no acabe):
            ## timer get time
        initialState = State(initialAlpha, initialBeta)
        finalState = State(alpha=wishedAlpha, beta=wishedBeta) 
        
        atOnceRotation = self._getRotationMatrix(initialState, finalState)
        
        Xpi_2ZgateSequence = self._logicalCompilationToXpi_2Z(atOnceRotation)
        
        pulseSeq = self._analogicalCompilationVZgates(Xpi_2ZgateSequence)
        
        measurementOutput = self.executeExperiment(pulseSeq)

        # initialize timer
        return measurementOutput
    
    def executeUnitaryOperation(self,operation):
        """It is the same method as executeRotation, but here one can directly
        input a desired operation on the initial ground state |0> of the qubit,
        instead of having to give in the output state
        
        Argument
        -------
        operation : numpyArray, shape(2,2) Complex
            The unitary operation one wishes to apply to the qubit, written as
            a matrix in the computational basis
        
        Outputs
        ------
        The measurement's result
        
        """
        Xpi_2ZgateSequence = self._logicalCompilationToXpi_2Z(operation)
        
        pulseSeq = self._analogicalCompilationVZgates(Xpi_2ZgateSequence)
        
        measurementOutput = self.executeExperiment(pulseSeq)  # if readout must be armed and if it is always interesting to make a readout
        
       # self.executeIQLOdriving(pulseSeq) # if instead arming of log channels is not necessary this should be done modularly
        
       # measurementOutput = self.executeReadout() 
        
        return measurementOutput
    
   
    def executeExperiment(self, pulseSeq):
        # Nota: en obras...
        #cuando ejecutas el exe del measurement, los dispositivos ya deben estar
        # de antemano en el Instrument server??? Aka, la comunicacion con los devices 
        # puede que si que sea necesaria hacerla por codigo
        # que es master channel?!
        Labber.ScriptTools.setExePath('C:\Program Files (x86)\Labber\Program')
        scriptPath = os.path.dirname(os.path.abspath(__file__))
        IQLOdriving = ScriptTools.MeasurementObject(os.path.join(scriptPath,\
                                                          'PulseSender.hdf5'))
        readout = ScriptTools.MeasurementObject(os.path.join(scriptPath,\
                                        str(self.readout_modality)+'.hdf5'))
        IQLOdriving.updateValue('',)
        
    
#    def executeExperiment(self,pulseSeq):
#        
#        devices = self._establishCommunication()
#        
#        self._armDevices(devices, pulseSeq)
#        
#        # self._armReadoutDevices(devices)
#        
#        triggerON
#        
#        return self._signalAnalyzer(devices)
#    
#    def _establishCommunication():
#        """Method that makes sure Labber instrument server is running and
#        establishes connection with the driving and readout devices.
#            Note that the exact model names and network interface adresses
#        are extracted from the external file "DeviceParams.txt".
#        Returns
#        ------
#        python dict with entries:
#            'AWG' : Labber.InstrumentClient object
#                An object associated to the specific AWG used
#            'LO' : Labber.InstrumentClient object
#                An object associated to the specific LO used
#            'client' : Labber.LabberBlockingClient object
#                The reference that controls the Instrument Server
#        
#        """
#        #TODO : Force Labber Instrument Server to initialize
#        client = Labber.connectToServer('localhost')
#        AWG = client.connectToInstrument(self.deviceSpecs['AWG'],
#                                         dict(interface=self.deviceSpecs['AWG_interf'],
#                                              adress=self.deviceSpecs['AWG_adress'], 
#                                              startup='Set config', lockbool=True),
#                                              bCreatenew=True)
#        LO = client.connectToInstrument(self.deviceSpecs['LO'],
#                                        dict(interface=self.deviceSpecs['LO_interf'],
#                                             adress=self.deviceSpecs['LO_adress'], 
#                                             startup='Set config', lockbool=True),
#                                             bCreatenew=True)
#        return {'AWG':AWG, 'LO':LO, 'client':client}
#    
#    def _armDevices(self,devices,pulseSeq):
#        """Method that sets the parameters of each device and orders them to 
#        wait until the trigger signal arrives in order to carry out them 
#        simultaneously.
#        
#        Arguments
#        --------
#        devices : python dictionary with entries:
#            'AWG' : Labber.InstrumentClient object
#                An object associated to the specific AWG used
#            'LO' : labber.InstrumentClient object
#                An object associated to the specific LO used
#            'client' : Labber.LabberBlockingClient object
#                The reference that controls the Instrument Server
#        
#        pulseSeq : python dictionary with entries:
#            'I' : NumpyArray, shape(1,) real doubles
#                The pulse amplitude samples to be sent to the AWG for the I
#                quadrature
#            'Q' : NumpyArray, shape(1,) real doubles
#                The pulse amplitude samples to be sent to the AWG for the Q
#                quadrature
#            't_AWG' : NumpyArray, shape(1,) real doubles
#                The time points for the pulses by the AWG
#            't_LO' : NumpyArray, shape(1,) real doubles
#                The time points for the pulses by the LO
#            'wLO' : double
#                The frequency that will be ordered to the local oscillator
#        
#        """
#        devices['AWG'].arm()
#        device['AWG'].setInstrConfig()
#        devices['AWG'].startInstrument()
#        # ein biko dana da armie AWG eta measurement devicek eta LOan triggerrari
#        # konekte, de forma que LO aktibetan dozunien se pone en marcha tol tinglao
#    
    
    def _getRotationMatrix(self, initialState, finalState):
        """
        Method that returns one of the possible unitary evolution matrixes, U 
        required to evolve the given initial state to the final state:
            
            U |initialState> = |finalState>
        
        Any unitary transormation on a Single qubit Bloch sphere can be seen 
        as a rotation of a certain angle theta (t) about a particular axis 
        defined by the real unit vector n=(nx,ny,nz) on the sphere. 
        
        Rn(t) = exp(-i t*n*sig/2) = cos(t/2)*I - i*sin(t/2)*(nx*X+ny*Y+nz*Z)
        
        where sig denotes the three component vector of Pauli matrices (X,Y,Z)
        Thus, the function solves the system of equations in the Comput. basis:
        Rn(t) |initialState> = |finalState>
        with unknowns t, n=(nx,ny,nz) , using numpy functions
        
        References
        ---------
        - Bibl. [3]
        
        Parameters
        ----------
        initialState : State
            The state of the wave-vector at time zero of the unitary evolution
        
        finalState : State
            The desired statevector inmediately post-evolution
            
        Returns
        -------
        U (singleUnitaryOperator) : NumpyArray, shape(2,2) (matrix on Complex)
            A matrix in the Computational Basis that expresses the single 
            transformation needed for the change |initialState> -> |finalState>
            
        """
        def Rnt(t,nx,ny,nz):
            return np.array([[np.cos(t/2)-1j*nz*np.sin(t/2),-np.sin(t/2)*(ny+1j*nx)]
                              ,[-np.sin(t/2)*(1j*nx-ny),np.cos(t/2)+1j*nz*np.sin(t/2)]])
    
        vect_ini = initialState.getAmplitudes().flatten()
        re_ini = vect_ini.real
        im_ini = vec_ini.imag
        vect_fin = finalState.getAmplitudes().flatten()
        re_fin = vect_fin.real
        im_fin = vect_fin.imag
        def Rnt(t,nx,ny,nz):
            return np.array([[np.cos(t/2)-1j*nz*np.sin(t/2),-np.sin(t/2)*(ny+1j*nx)]
                   ,[-np.sin(t/2)*(1j*nx-ny),np.cos(t/2)+1j*nz*np.sin(t/2)]])  
        def F(unknowns):
            (t,nx,ny,nz) = unknowns
            A = np.array([[re_ini[0],im_ini[0],-re_ini[1],im_ini[1]],
                         [im_ini[0],-re_ini[0],-im_ini[1],-re_ini[1]],
                         [re_ini[1],-im_ini[1],re_ini[0],im_ini[0]],
                         [im_ini[1],re_ini[1],im_ini[0],-re_ini[0]]])
            x = np.array([[np.cos(t/2)],[0],[0],[0]])+np.sin(t/2)*np.array([[0],[nz],[ny],[nx]])
            B = np.array([[re_fin[0]],[im_fin[0]],[re_fin[1]],[im_fin[1]]])
            return (np.dot(A,x)-B).flatten().tolist()
        
        solut = opt.fsolve(F,[np.pi/4,1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)])
        #norm = np.sqrt(solut[1]**2+solut[2]**2+solut[3]**2)
        #(solut[1],solut[2],solut[3])=(solut[1]/norm,solut[2]/norm,solut[3]/norm)
        return Rnt(solut[0],solut[1],solut[2],solut[3])
        
            
    def _logicalCompilationToXpi_2Z(self,unitary):
        """Given an input unitary transormation expressed as a matrix in the 
        computational basis, this rutine obtains the parameters t (theta),
        g (gamma), p (phi), l (lambda), that satisfy the following theorem:
            Any arbitrary SU(2) gate can be written as:
                U(t,g,p,l) = [[cos(t/2), -i*exp(i*l)*sin(t/2)],
                              [-i*exp(i*p)*sin(t/2), exp(i*(l+p))*cos(t/2)]]
            which up to a global phase is represented as:
                U(t,p,l) = Z(p) X (t) Z(l)
            where Z(j),X(j) represent rotations about Z and X of angle j
            Gate identies yield:
                U(t,p,l) = Z(p-pi/2) X (pi/2) Z(pi-t) X(pi/2) Z(l-pi/2)
            Thus, knowing t,p,l we are able to compile any single qubit
            operation into a sequence of two pi/2 X gates and Z gates (that
            ideally can be implemented as Virtual Z gates, s.t. effectively
            the only driving pulses required to perform the operation 
            will be two pi/2 X pulses)
        
        References
        -------
        - Bibl. [2]
        
        Arguments
        --------
        unitary : NumpyArray, shape(2,2) Complex
            It is the unitary operation we want to translate
            
        Output
        ------
        : tuple (a,b,c,d,e) real doubles
            The rotation angles for the respective (Z,X,Z,X,Z) gate tuple of
            the theorem. Note b,d are pi/2 always
        
        """
        def F(unknowns):
            (t, g, p, l, d1,d2,d3,d4) = unknowns # d1,..d4 are dummie variables to satisfy the numpy function
            f1=np.cos(t/2)*np.cos(g)-unitary.real[0,0]
            f2=np.cos(t/2)*np.sin(g)-unitary.imag[0,0]
            f3=np.sin(t/2)*(np.cos(p)*np.sin(g)+np.sin(p)*np.cos(g))-unitary.real[1,0]
            f4=np.sin(t/2)*(np.sin(p)*np.sin(g)-np.cos(p)*np.cos(g))-unitary.imag[1,0]
            f5=np.sin(t/2)*(np.cos(l)*np.sin(g)+np.sin(l)*np.cos(g))-unitary.real[0,1]
            f6=np.sin(t/2)*(np.sin(l)*np.sin(g)-np.cos(l)*np.cos(g))-unitary.imag[0,1]
            f7=np.cos(t/2)*(np.cos(l+p)*np.cos(g)-np.sin(l+p)*np.sin(g))-unitary.real[1,1]
            f8=np.cos(t/2)*(np.cos(l+p)*np.sin(g)+np.sin(l+p)*np.cos(g))-unitary.imag[1,1]
            return [f1,f2,f3,f4,f5,f6,f7,f8]
        solut = opt.fsolve(F,[np.pi/4,np.pi/4,np.pi/4,np.pi/4,1,1,1,1])
        (t,g,p,l)=solut[:5]
        return (p-np.pi/2, np.pi/2, np.pi-t, np.pi/2, l-np.pi/2)
        
    #(Zgate(p-np.pi/2),Xgate(np.pi/2),Zgate(np.pi-t),Xgate(np.pi/2),Zgate(l-np.pi/2))
    
        #decomposed=np.array([[np.cos(t/2), -1j*np.exp(1j*l)*np.sin(t/2)],
        #        [-1j*np.exp(1j*p)*np.sin(t/2), np.exp(1j*(l+p))*np.cos(t/2)]])
        #print("\ndecomposed matrix yields recomposed\n", decomposed,
        #               "\n with the global phase \n", decomposed*np.exp(1j*g))
    
    def _analogicalCompilation(self,ZXZXZdecompositionAngles):
        """Method that prepares the I,Q pulses for the AWG and the frequency
        of the LO, so they can be sent to the devices.
        
        It is called analogical compilation since it translates the sequence
        of Z(p-pi/2) X (pi/2) Z(pi-t) X(pi/2) Z(l-pi/2) gates outputed by the
        logical compiler into the actual microwave pulses that will drive the
        rotation in the real transmon qubit.
        In particular, note that it is here where virtual Z gating is applied
        by forgetting about the three Z gates, and instead adding their 
        rotation angles as relative phases to the pi/2 X pulses, as described
        in the bibliography paper [1].
        More in particular still, as before and after this gate sequence we
        will not apply any additional operation, the only Z gate that is
        relevant is the one in between the two finite time X gates, and thus
        the only angle we shall care about of the ZXZXZ decomposition is the
        theta described in the logical compilation function. Thus, for our
        purposes:
        Z(p-pi/2) X (pi/2) Z(pi-t) X(pi/2) Z(l-pi/2) ~ X (pi/2) Z(pi-t) X(pi/2)
        And thus, this can be implemented by just two pulses:
            X(pi/2 rot, with pi-t phase relative to) X(pi/2 rot)
        
        As DRAG pulsing may be used, the pulse envelope is set to be a Gaussian
        for convenience with the derivative, following its use in bibl. [2].
        The X pulses are designed following the indications in bibl. [1]
        by sending the drivng signal entirely in the I quadrature, while the
        drag part occupies the Q channel.
        
        References
        ---------
        - Bibl. [1] and [2]
        
        Arguments
        --------
        ZXZXZdecompositionAngles : tuple (a,b,c,d,e) real doubles
            The rotation angles for the respective (Z,X,Z,X,Z) gate tuple of
            the U(a,b,c,d,e) = Z(a) X (b) Z(c) X(d) Z(e) decomposition with
            b,d= pi/2 in all cases
        
        Outputs
        -------
        python dictionary 
            Containing the following data:
            'I' : NumpyArray, shape(1,) real doubles
                The pulse amplitude samples to be sent to the AWG for the I
                quadrature
            'Q' : NumpyArray, shape(1,) real doubles
                The pulse amplitude samples to be sent to the AWG for the Q
                quadrature
            't_AWG' : NumpyArray, shape(1,) real doubles
                The time points for the pulses by the AWG
            't_LO' : NumpyArray, shape(1,) real doubles
                The time points for the pulses by the LO
            'wLO' : double
                The frequency that will be ordered to the local oscillator
        
        """
        # en obras....
        # definir V0 necesaria, sigma necesaria y buffer between pulses necesaria!
        T_pi_2 = # NEED TO GET THE ADEQUATE PULSE LENGTH FOR A PI/2 PULSE EXPERIMTLY
        bufferBetweenPulses = # NEED THIS AS WELL!
        wLO = self.circuitParameters['QubitFrequency'] - self.deviceSpecs['AWGwssb'] # wssb+wLO = wqubit
        V0 =
        
        def GaussianEnvelope(t,V0,T,sigma=T/4):
            """Method that returns the amplitude of a Gaussian Envelope 
            with parameters V0,T,sigma for a time t where:
                -t is the time of interest
                -V0 is the scalar amplitude of the pulse
                -T is HALF the pulse length, usualy set to T = 4*sigma
                -sigma is the standard deviation of the pulse
            """
            return V0*(np.exp(-(t-T)**2/(2*sigma**2))-np.exp(-T**2/(2*sigma**2)))/(1-np.exp(-T**2/(2*sigma**2)))
        
        def DerivativeGaussianEnvelope(t,V0,T,sigma=T/4):
            return -V0*(T-t)*np.exp(-(T-t)**2/(2*sigma**2))/(sigma**2*(np.exp(-T**2/(2*sigma**2))-1))
        
        timeLO = np.linalg(0, 2*T_pi_2 + bufferBetweenPulses + 2*T_pi_2, 
                         self.machineSpecs['LO_rateSample'])
        timeAWG = np.linalg(0, 2*T_pi_2 + bufferBetweenPulses + 2*T_pi_2, 
                         self.machineSpecs['AWG_rateSample'])
        
        # I pulses drive rotations in X while Q pulses in Y
        # in our case, as we only need X pi/2 pulses, if drag is off, Q channel
        # will be null at all times. If drag is on, as the drive is done
        # through the I channel, the derivative pulse will go on Q
        Ipulses=[]
        Qpulses=[]
        for t in timeAWG:
            if t < 2*T_pi_2:
                Ipulses.append(GaussianEnvelope(t,V0,T)*np.sin(wssb*t))
            elif t < 2*T_pi_2 + bufferBetweenPulses:
                Ipulses.append(0.0)
            else:
                Ipulses.append(GaussianEnvelope(t,V0,T)*np.sin(wssb*t+
                               ZXZXZdecompositionAngles[2]))
            if (self.drag):
                if t < 2*T_pi_2:
                    Qpulses.append(self.dragBeta*DerivativeGaussianEnvelope(t,V0,T)*
                                   np.cos(wssb*t))
                elif t < 2*T_pi_2 + bufferBetweenPulses:
                    Qpulses.append(0.0)
                else:
                    Qpulses.append(self.dragbeta*DerivativeGaussianEnvelope(t,V0,T)*
                                   np.cos(wssb*t+ZXZXZdecompositionAngles[2]))
            else:
                Qpulses.append(0.0) 
                # si Q channel es full ceros para evitar errores podria
                # haerse que simplemente no se outputease nada a ese channel
        # DEPHASING para Virtual Z! alzu sartun AWGren parametrotzat edo 
        # I eta Qe malieu biher dozuz?
        return {'I':Ipulses, 'Q':Qpulses, 't_LO':timeLo, 't_AWG':timeAWG, 'wLO':wLO}
        
    def _extractParamsFromFile(self,file_name):
        #DEVICE NAMESENA IMPLEMENTE!
        """Method that reads a data file named "file_name" and extracts the 
        parameters and respective numeric values listed there in the following
        format:
            NameOfParameter numericValue
            NameOfParameter2 numericValue
        where the numeric value is expected to be given in any of the following
        formats:
            8
            2.71828
            8.23e-07
            8e-7
        Collected data are registered in a Python dictionary.
        
        Argument
        --------
        file_name : string
        
        Output
        -----
        params : python dictionary
            with key='NameOfParameter', value= NumericValue (float)
        
        """
        datafile = open(file_name,"r")
        params={}
        for line in datafile:
            name=''
            i=0
            while(line[i]!=' '):
                name+=line[i]
                i+=1
            i+=1
            number=''
            while(line[i]!='e' and line[i]!='\n'):
                number+=line[i]
                i+=1
            number=float(number)
            if(line[i]=='e'):
                i+=1
                exponent=''
                while(line[i]!='\n'):
                    exponent+=line[i]
                    i+=1
                exponent=float(exponent)
                params[name] = number*10**exponent
            else:
                params[name] = number
        datafile.close()    
        return params
        
    def _calibrateDragBeta(self):
        # Given omega, V0 and anharmonicity a, and the correct pi/2 pulse
        # leakage measurements swept with varying beta in [-2,2]
        # A test could be made by applying an identity roation, 
        # decomposed in two pi/2 pulses and see the resulting leakage
        # the best record will correspond to the optimized beta
        
    def _calibrateXPulse(self):
        # dados los parametros OMEGA y V0 calcula la T ideal para pi/2 y la buffer