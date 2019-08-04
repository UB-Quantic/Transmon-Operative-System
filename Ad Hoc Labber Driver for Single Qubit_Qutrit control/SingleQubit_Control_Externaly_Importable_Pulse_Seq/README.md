# INSTALLATION
------------
If you wish to use this driver pack in Labber, just copy the folder named 
"SingleQubit_Control_Externaly_Importable_Pulse_Seq", containing the files
listed in the section "Files in the Package", in the Labber Custom Driver directory.
If working on Windows, the path is probably "C:\Users\UserName\Labber\Drivers"
but anyhow, the easy way to know it is: in Labber Instrument Server, click
Edit > Perferences > Folders, then the path in "Local drivers" is what you are
looking for.

After that, in order to be able to import the PulseSequencetoTxtConversor class
to your scripts wherever these scripts are located, open your pyhton interpreter
with working directory: the path described in the earlier paragraph, then execute the
script "AddPackageToPath.py". This script just adds the package to path with these
commands:

	import os, sys
	CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
	sys.path.append(os.path.dirname(CURRENT_DIR))

Note that this addition to path is not permanent, so it should be executed whenever
one intends to use the class from an external directory.
If this seems rather complicated, the easy way is always to run your script in the
same driver folder, so you can directly import tha PulseSequencetoTxtConversor 
class without worrying about the path stuff. See the example script for an example
of this way.

# FUNCTIONALITIES AND PURPOSE
--------------------------
This driver pack is intended to give a Superconducting qubit engineer the
ability to generate in a python Script a customized set of I,Q pulses
for AWGs in the form of numpy arrays and send them correctly to the signal generating
hardware using Labber. Employing the pyhton class attached, the numpy arrays may 
be translated to a txt file that will then be interpreted by the Labber Signal
Generating driver, "SingleQubit_Control_Externaly_Importable_Pulse_Seq", which
using the Labber Signal Communication protocol will ultimately send the custom pulse
information to the actual AWG hardware.

In short, the main idea with this pack is to give the ability to send I Q signal pulses
designed by numpy scripting to the particular AWG used in the experiments.

The easy way of doing this would have probably been to simply make a custom Signal
Generator Labber Driver ad hoc, but what we have done here goes beyond that.
The default SingleQubit_Signal_Generator driver of Labber has been used as a template 
to build up the functionality, so it is an add-on to the default signal generator.

This modified version of SingleQubit_Signal_Generator allows to continue designing
sequences of pulses with the given options in the configuration window but with the
very relevant option to intersperse numpy script designed arbitrary waveforms among them.

In particular we have developed the following options:

- One may continue using this drivers as the good old SingleQubit_Signal_Generator by
unchecking the "Import pulse sequence from external file" option in the driver config window.

- If you instead wish to import pulses externaly, then there are three options:

  > Ignore the pulses outputed by the SingleQubit_Signal_Generator options in the config 
 window and just output the array info imported externally. For that, check the box in the
 window that says "Ignore the pulses designed by the driver and output the ones of the .txt"

  > Mix the pulses designed with the SingleQubit_Signal_Generator like interface with the
 imported ones. Here you might intersperse the imported pulses rather after the pre-pulses
 and before the rest of pulses, or after the pre-pulses and after the main pulses but before
 the readout and tomography pulses.

Note that these three options allow all the possible combinations of pulse sequences:
 	> Just driver pre-config pulses (pre-pulses, main pulses and readout/tomography pulses)
	> Just imported pulses
	> (Pre-pulses) + imported ones + (main pulses) + (readout,tomography pulses)
   	where by parenthesis we mean that they can be added or not at wish.
      Note that as pre-pulses can be chosen to overlap with the rest of the pulses, we even allow
   the option to superimpose external and driver given pulses. Options are unlimited!


# HOW TO IMPORT CORRECTLY EXTERNAL I,Q PULSES
-------------------------------------------
## TXT FILE
  
The drivers, when asked to import pulses from an external file, expect the pulse
information to be written in a txt file with the following exact layout:

"""
	SampleRate <double_in_s^-1>
	OutputNumber <integer>
	I1 <double_1> <double_2> .... <double_m1>
	Q1 <double_1> <double_2> .... <double_m1>
	I2 <double_1> <double_2> .... <double_m2>
	Q2 <double_1> <double_2> .... <double_m2>
    		...
	In <double_1> <double_2> .... <double_mn>
	Qn <double_1> <double_2> .... <double_mn>

"""
- where SampleRate is the inverse of the time period per defined
point in the arrays (ideally this should coincide with the sample
rate of the arbitrary waveform generator AWG to be used);

- mi is the maximum number of points in each array. If mi are 
different, every array will be filled until they reach a length
equal to max({mi})=m with zeros. In fact m*(SampleRate^-1) is equal
to the total time length of the train of pulses;

- n must be equal to OutputNumber, where  nOutputNumber refers to
the number of output IQ channels. If it is larger than the number of
outputs designed in the configuration window, then the latter will be
forcely changed.

Note that SampleRate must concide with the sample rate introduced in the
driver configuration. If not, the sampleRate in the txt will be forced.

## UTILITY TO CONVERT NUMPY ARRAYS INTO THE TXT

In order not to worry about the layout of the txt file, the package includes a python class
called PulseSequencetoTxtConversor that allows direct generation of a correctly formatted txt
out of an input python dictionary with the desired I,Q pulse numpy arrays.

Example of use:

Given d is a dictonary with the following entries:

	'SampleRate': double
        'OutputNumber': integer
        'I1': numpyArray shape(m1,)
        'Q1': numpyArray shape(m1,)
            ....
        'Qn': numpyArray shape(mn,)

- where SampleRate is the inverse of the time period per defined
point in the arrays (ideally this should coincide with the sample
rate of the arbitrary waveform generator AWG to be used);
- mi is the maximum number of points in each array;
- n must be equal to OutputNumber, where  nOutputNumber refers to
the number of output IQ channels.

Then:
 pyhton code:

	import PulseSequencetoTxtConversor
	conversor = PulseSequencetoTxtConversor()
	conversor.convertDicArraysToTxt( d )

should succesfully generate the necessary txt in the path where the pyhton
script of class PulseSequencetoTxtConversor is located. This txt must be placed
in the "SingleQubit_Control_Externaly_Importable_Pulse_Seq" driver directory as explained
in the installation section. Next the drivers will be able to generate correctly the desired
pulse sequence.


# FILES IN THE PACKAGE
-------------------
 - SingleQubit_Control_Externaly_Importable_Pulse_Seq.py

 - SingleQubit_Control_Externaly_Importable_Pulse_Seq.ini

	The above pair are the two driver files used by Labber

 - Pulses_for_the_Signal_Generator.txt
	This is the default file from where the pulses will be extracted. By default the package
	includes in this file an example sequence.

 - PulseSequencetoTxtConversor.py
	This is the conversor offered to authomatize the generation of a correct txt out of a
	python dictionoray with the desired numpy array pulses

 - exampleScript.py
	This is an example script that generates the I,Q pulses for a sequence of two X(theta) rotation
	pulses. It assumes that the package is correctly installed in the cusotm driver path.

 - AddPackageToPath.py
	This is a simple script that, adds the path where the script is to your Python path for
	the running interpreter.

 - This README.md