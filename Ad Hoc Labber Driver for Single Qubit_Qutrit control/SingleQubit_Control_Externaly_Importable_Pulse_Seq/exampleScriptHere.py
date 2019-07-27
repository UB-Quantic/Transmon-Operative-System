# PARA RAFA:

import numpy as np
import matplotlib.pyplot as plt
import Labber as Labber # si te da error es que no has instalado el api de 
                        # Labber en tu path, asi que lo que puedes hacer es 
                        # borrar esta y las ultimas tres lineas de codigo del script
                        # que precisamente usaban esta libreria
                        # no hace falta, puedes controlar labber manualmente
                        
# respecto al import siguiente, para poder hacerlo,
# primero ejecuta dentro de la carpeta del driver customizado, una vez puesto donde toca
# de Labber, el script AddPackageToPath.py tal y como se comenta en la seccion
# de instalacion del README
# Esto es solo para que puedas importar PulseSequencetoTxtConversor directamente
# como a continuacion (tengo que mejorar esto, lo se...)
# si no te sale bien simplemente hazlo todo en la misma carpeta de drivers customizados
# en donde toca de Labber y en vez de la linea siguiente pon esta que está comentada:
# import PulseSequencetoTxtConversor

from SingleQubit_Control_Externaly_Importable_Pulse_Seq import PulseSequencetoTxtConversor

#los parametros de los dos pulsos

T = 10e-9 # la mitad de la longitud de cada pulso gausiano tal que T=4sigma por default
sampleRate=5e9
bufferBetweenPulses = 2e-9
wssb= 500e6*2*np.pi
V0 = 1
phase_2_pulse = 3.1415/2 # phase of the second pulse w respect to the first one
beta=0.5e-9 #drag beta
n_points=(2*T + bufferBetweenPulses + 2*T)*sampleRate #total length of the sequence/time divisions
drag=True

# definimos las enevlopes que usaremos
def GaussianEnvelope(t,V0,T,sigma='default'):
    """Method that returns the amplitude of a Gaussian Envelope 
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

# los puntos temporales
timeAWG = np.linspace(0, 2*T + bufferBetweenPulses + 2*T, num=int(n_points))

Ipulses=[]
Qpulses=[]
# llenamos los vectores de puntos de los pulsos
for t in timeAWG:
    if t < 2*T:
        Ipulses.append(GaussianEnvelope(t,V0,T)*np.sin(wssb*t)) # first pulse (I channel)
    elif t < 2*T + bufferBetweenPulses:
        Ipulses.append(0.0) # buffer in between
    else:
        Ipulses.append(GaussianEnvelope(t-(2*T + bufferBetweenPulses),V0,T)*np.sin(wssb*(t-(2*T + bufferBetweenPulses))+phase_2_pulse)) # second pulse (I channel)
    
    if (drag==True): # depending of wheter we want drag or not
        if t < 2*T:
            Qpulses.append(beta*DerivativeGaussianEnvelope(t,V0,T)*np.cos(wssb*t)) # first pulse (Q channel)
        elif t < 2*T + bufferBetweenPulses:
            Qpulses.append(0.0) # buffer
        else:
            Qpulses.append(beta*DerivativeGaussianEnvelope(t-(2*T + bufferBetweenPulses),V0,T)*np.cos(wssb*t+phase_2_pulse)) # second pulse (Q channel)
    else:
        Qpulses.append(0.0) 

#ploteamos aqui mismo para comprobar que funciona
plt.plot(timeAWG, Ipulses)
plt.plot(timeAWG, Qpulses)

# preparamos el dictionary para el generador del txt
pulseSeq={}
pulseSeq['Q1'] = Qpulses
pulseSeq['I1'] = Ipulses
pulseSeq['SampleRate'] = sampleRate
pulseSeq['OutputNumber'] = 1

# genera el txt
converter = PulseSequencetoTxtConversor.PulseSequencetoTxtConversor()
converter.convertDicArraysToTxt(pulseSeq)

#ahora si le das a start en Labber al driver deberia salir la secuencia correcta
# en teoria los siguientes comandos deberian funcionar si tienes abierto el 
# instrument server y instalado el API de Python, si no da igual, a mano, que esto es pa probar

client = Labber.connectToServer('localhost')
pulseGenerator = client.connectToInstrument('SingleQubit Pulse Generator with Externaly Importable Pulse_Seqs',dict(name='bai',interface='None',address='',startup='Set config', lockbool=False),bCreateNew=True)
pulseGenerator.startInstrument()
# si lo haces asi para ver el plot tienes que darle get trace, si no solo lo ejecuta internamente

        

