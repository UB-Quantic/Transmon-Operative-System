import numpy as np
import matplotlib.pyplot as plt
import Labber as Labber 
"""
 if an error arises it is probably due to the fact that you haven't installed
the Labber api properly. The only use of the api in this script are the last three 
lines of code. Nevertheless, these are not really necessary, so you can just 
erase them, and execute the Instruments manually
"""

import PulseSequencetoTxtConversor

# the parameters of the two X pulses are set as:

T = 10e-9 # half of the duration of each Gaussian pulse; T=4sigma per default
sampleRate=5e9
bufferBetweenPulses = 2e-9
wawg= 500e6*2*np.pi
V0 = 1
phase_2_pulse = 3.1415/2 # phase of the second pulse respect to the first one
beta=0.5e-9 #drag beta
n_points=(2*T + bufferBetweenPulses + 2*T)*sampleRate #total length of the sequence aka #time divisions
drag=True

# we define the envelope functions
def GaussianEnvelope(t,V0,T,sigma='default'):
    """Function that returns the amplitude of a Gaussian Envelope 
    with parameters V0,T,sigma for a time t where:
        -t is the time of interest
        -V0 is the scalar amplitude of the pulse
        -T is HALF the pulse length, usualy set to T = 4*sigma
        -sigma is the standard deviation of the pulse, per default it's T/4
    """
    if sigma=='default':
        sigma=T/4
    return V0*(np.exp(-(t-T)**2/(2*sigma**2))-np.exp(-T**2/(2*sigma**2)))/(1-np.exp(-T**2/(2*sigma**2)))

def DerivativeGaussianEnvelope(t,V0,T,sigma='default'):
    if sigma=='default':
        sigma=T/4
    return -V0*(T-t)*np.exp(-(T-t)**2/(2*sigma**2))/(sigma**2*(np.exp(-T**2/(2*sigma**2))-1))

# the time array
timeAWG = np.linspace(0, 2*T + bufferBetweenPulses + 2*T, num=int(n_points))

Ipulses=[]
Qpulses=[]
# we fill the pulse "arrays" with voltage amplitudes
for t in timeAWG:
    if t < 2*T:
        Ipulses.append(GaussianEnvelope(t,V0,T)*np.sin(wawg*t)) # first pulse (I channel)
        Qpulses.append(-GaussianEnvelope(t,V0,T)*np.cos(wawg*t)) # first pulse (Q channel)

    elif t < 2*T + bufferBetweenPulses:
        Ipulses.append(0.0) # buffer in between
        Qpulses.append(0.0)
    else:
        Ipulses.append(GaussianEnvelope(t-(2*T + bufferBetweenPulses),V0,T)*np.sin(wawg*(t-(2*T + bufferBetweenPulses))+phase_2_pulse)) # second pulse (I channel)
        Qpulses.append(-GaussianEnvelope(t-(2*T + bufferBetweenPulses),V0,T)*np.cos(wawg*(t-(2*T + bufferBetweenPulses))+phase_2_pulse)) # second pulse (Q channel)
   
if(drag==True):
    for i,t in enumerate(timeAWG):
        if t < 2*T:
            Ipulses[i]+= beta*DerivativeGaussianEnvelope(t,V0,T)*np.sin(wawg*t+np.pi/2) # first pulse (I channel)
            Qpulses[i]+= -beta*DerivativeGaussianEnvelope(t,V0,T)*np.cos(wawg*t+np.pi/2) # first pulse (Q channel)

        elif t > 2*T + bufferBetweenPulses:
            Ipulses[i]+= beta*DerivativeGaussianEnvelope(t-(2*T + bufferBetweenPulses),V0,T)*np.sin(wawg*(t-(2*T + bufferBetweenPulses))+np.pi/2+phase_2_pulse) # second pulse (I channel)
            Qpulses[i]+= -beta*DerivativeGaussianEnvelope(t-(2*T + bufferBetweenPulses),V0,T)*np.cos(wawg*(t-(2*T + bufferBetweenPulses))+phase_2_pulse+np.pi/2) # second pulse (Q channel)
       

# plot the designed pulses to see if they are the desired ones
plt.plot(timeAWG, Ipulses)
plt.plot(timeAWG, Qpulses)

# we prepare the dictionary for the txt generator class
pulseSeq={}
pulseSeq['Q1'] = Qpulses
pulseSeq['I1'] = Ipulses
pulseSeq['SampleRate'] = sampleRate
pulseSeq['OutputNumber'] = 1

# generathe the txt data file
converter = PulseSequencetoTxtConversor.PulseSequencetoTxtConversor()
converter.convertDicArraysToTxt(pulseSeq)

"""
At this point if you click start in the Labber driver, the correct sequence
will be imported and shown.

For the following lines of code, they should work only if you have the Instument
Server open and the LAbber python api installed. Anyhow, they can be omitted and
manually start the drivers
"""

client = Labber.connectToServer('localhost')
pulseGenerator = client.connectToInstrument('SingleQubit Pulse Generator with Externaly Importable Pulse_Seqs',dict(name='bai',interface='None',address='',startup='Set config', lockbool=False),bCreateNew=True)
pulseGenerator.startInstrument()
# if you use these lines to start the instrument remember to click on get Trace if you want to visualize the generated signal       
