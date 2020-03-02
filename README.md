## Signal Processing
# ENE418 Project 1

# Part 1 Fault detection of a wind-turbine pitch drive based on measured current signals

A scaled wind turbine pitch drive is built by Motor – ABB 2.2 kW, 4 pole, rated at 1420 rpm, 7.2 Nm and Drive – ABB FOC ACS 355 (Operated in FOC mode) as shown below. Current signals from the system are collected by Sensors – LEM LTS 6NP (Hall effect Current sensors). The bearing - SKF 6205 is used in the system.

To detect localized faults on the system, the phase currents are measured in healthy and faulty conditions and stored in files (File contents: Three phase currents (ia, ib, ic) and time). Download measured data from files
Data_Healthy_1420RPM_48p8Hz_70Load.mat
fault 1.mat, fault 2.mat, fault 3.mat
• You can freely choose any signal processing method or suggest any algorithm to find a localized fault in the pitch motor.
• Calculate characteristic frequencies for faulty stator winding and bearings
• The method must be implemented in Matlab® environment. After implementing the algorithm, you should evaluate its performance.
• To compare with the healthy case, a common or comparison plot of healthy and faulty cases can be of interest.

# Part 2 Implementation and Evaluation of a Sinusoid Extraction Algorithm
• You can freely choose one of the signal processing algorithms (“current reference generator”) from any reference or suggest reference. The algorithm should be able to extract a fundamental sinusoid from a severely distorted waveform.
• The algorithm must be implemented in Matlab® environment (do not use Simulink).
• After implementing the algorithm, you should evaluate its performance with the following test cases:
o Your primary test signal contains six odd-order harmonics (3rd, 5th, 7th, 9th, 11th, and 13th), and their amplitudes are 15% of the fundamental’s amplitude.
o Three fundamental frequencies are evaluated in all experiments (49, 50, and 51 Hz).
ENE 418 - miniproject – signal processing
o Experiment 1: Let the algorithm run long enough to be in steady state. Calculate the total harmonic distortion, THD %, as well as individual distortions corresponding to each of the six harmonics at the output of your signal processing system. Compare those measures to the corresponding distortions at the input of your system. Study also the phase delay between your output and the fundamental input—is it appropriate for a current-reference-generator application?
o Experiment 2: Let the algorithm run long enough to be in steady state. Then increase the amplitudes of the fundamental and all harmonics by 50%; and after reaching steady state again, decrease the amplitudes back to their original values. Study the behavior and response times of your algorithm in those step experiments.
o Experiment 3: Use first an input signal that contains only the fundamental, and let the algorithm run long enough to be in steady state. Then switch on the six harmonics. Study the behavior and response time of your algorithm.
o You may create additional experiments if you like.


