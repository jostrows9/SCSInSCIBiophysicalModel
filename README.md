# SCSInSCIBiophysicalModel

### Overview 

Epidural electrical spinal cord stimulation (SCS) can improve motor control in animal models and humans with spinal cord injury. The current understanding is that the stimulation excites sensorimotor circuits below the injury by interacting with the natural activity of large afferent fibers entering into the spinal cord from the dorsal roots. Here, we built a computational framework to study the mechanisms of this interaction.

This repository contains the code underlying the simulations performed in Balaguer et al. (2025), and was built upon several previously-validated biophysical models (references provided below).

### Getting Started 

All simulations are run in NEURON. Instructions for setting up NEURON can be found [here](https://www.neuron.yale.edu/neuron/download). Once installed, NEURON must be configured by running the following from the home directory: 

```nrnivmodl mod_files```

The rest of the required packages can be found in `requirements.txt` and can be installed using:

```pip install -r requirements.txt``` 

### Scripts

This repository includes three scripts which are the generic versions of the main simulations included in Balaguer et al., 2025. A general description of these scripts, their use case, and the corresponding figure(s) can be found below. Information on the exact use of each script can be found in the function descriptions. 

- `record_mn_membrane_potential_supraspinal_scs_simulation`: Records a motoneuron membrane potential with (1) supraspinal input alone, (2) SCS input alone, and (3) the combination of supraspinal and SCS input. This script can be used to visualize how supraspinal and SCS inputs interact in the motoneuron membrane, as in Figure 2a-c, Figure 7a, and Figure 8b (Balaguer et al., 2025).
- `run_mn_pool_supraspinal_scs_simulation`: Simulates a motoneuron pool (N=100) with continuous SCS and supraspinal inputs. This script underlies the bulk of the simulations in Balaguer et al., 2025 (Figure 2e-h, Figure 3, Figure 4a-k, Figure 7b-c, and Figure 8c-f).
- `run_fine_motor_task_simulation`: Simulates a fine motor task, during which the brain input fluxuates as a sinusoidal wave. This script can be used to visualize how SCS improves performance in a fine motor task after lesion, as in Figure 4l-n and Figure 8d (Balaguer et al., 2025).

Each script can be run from the home directory with the following: 

```python modeling_scripts/[script_name].py``

This repo also contains functions to perform the analysis of EMG data seen in Balaguer et al., 2025. These functions are located in `analysis_scripts/emg_analysis.py`: 

- `calculate_p2p_amplitude`: Calculates the peak-to-peak amplitude of EMG responses to stimulation.
- `compute_isi_distribution`: Returns computed interspike-interval (ISI) histogram.

### References 

Booth, V., Rinzel, J. & Kiehn, O. Compartmental model of vertebrate motoneurons for Ca2+-dependent spiking and plateau potentials under pharmacological treatment. J. Neurophysiol. 78, 3371–3385 (1997).

Formento, E. et al. Electrical spinal cord stimulation must preserve proprioception to enable locomotion in humans with spinal cord injury. Nat. Neurosci. 21, 1728 (2018).

Greiner, N. et al. Recruitment of Upper-Limb Motoneurons with Epidural Electrical Stimulation of the Primate Cervical Spinal Cord. Nat. Commun. 12, (2021).

McIntyre, C. C. & Grill, W. M. Extracellular Stimulation of Central Neurons: Influence of Stimulus Waveform and Frequency on Neuronal Output. J. Neurophysiol. 88, 1592–1604 (2002).

