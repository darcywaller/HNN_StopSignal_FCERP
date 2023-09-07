# Biophysical microcircuit mechanisms of frontocentral evoked potentials support a race model interpretation of the Stop-Signal P3

Diesburg, D. A.<sup>1</sup>, Wessel, J. R. <sup>2, 3</sup> & Jones, S. R. <sup>1, 4</sup>  
   
<sup><sup>1</sup> Department of Neuroscience, Carney Institute for Brain Sciences, Brown University, Providence, United States  
<sup>2</sup> Department of Psychological and Brain Sciences, University of Iowa, Iowa City, Iowa, United States  
<sup>3</sup> Department of Neurology, Carver College of Medicine, University of Iowa Hospitals and Clinics, Iowa City, Iowa, United States  
<sup>4</sup> Center for Neurorestoration and Neurotechnology, Providence VAMC, Providence, United States</sup>  

***
This repository contains code to replicate the main findings associated with the manuscript [“Biophysical microcircuit mechanisms of frontocentral evoked potentials support a race model interpretation of the Stop-Signal P3”](biorxiv link goes here). Data used to compute empirical ERPs to which model results were compared were drawn from the OSF repository associated with [Wessel, 2020](https://osf.io/v3a78/). 

Abstract:
*Human frontocentral event-related potentials (FC-ERPs) are ubiquitous neural correlates of cognition and control, but their generating multiscale mechanisms remain mostly unknown. We used Human Neocortical Neurosolver(HNN)’s biophysical model of the canonical neocortical column under exogenous drive (Neymotin et al., 2020) to estimate the cell and circuit mechanisms underpinning the Stop-Signal-locked FC-ERP during action-stopping in the Stop-Signal task. We first demonstrated that a conserved (from primary sensory cortices, e.g., Jones et al., 2007) sequence of simulated external proximal and distal drives can produce the FC-ERP. We then used this model of the FC-ERP to examine mechanisms underlying condition differences in FC-ERP features between successful and failed action-stopping and in particular their adherence to the predictions of the behavioral “race model” (Logan and Cowan, 1984). These modeling results provide novel testable predictions of the thalamocortical mechanisms underlying action-stopping and demonstrate the utility of biophysical models such as HNN in the investigation of evoked responses associated with higher-order cognition.*  
  

***

**Please Note:** The specific version of HNN used in this study has not been released at the time of publication. The version used here differs from the current (May 2021) release of HNN GUI in the way layer V calcium dynamics are calculated, which can lead to slight differences in the shape of the dipole waveform.  
For an exact replication of the published simulations, please use the parameter files in the 'HNN Parameters' directory (described below) and replace the file  'L5_pyramidal.py' in your local HNN GUI directory with the [L5_pyramidal.py](https://github.com/kohl-carmen/HNN-AEF/blob/main/L5_pyramidal.py) file provided by Kohl et al. (2021).


***
## Data
The empirical EEG data used in this study to create the FC-ERPs compared to model output was published as part of a previous investigation ([Wessel (2020)](https://www.jneurosci.org/content/40/2/411)). The data is openly available on [OSF](https://osf.io/v3a78/) and therefore not included it in this repository. To download the data needed to replicate these analyses, simply download the 'EEG data' folder available in Wessel's OSF repository. (However, note that this directory contains 235 EEG datasets and is therefore somewhat large. If you simply wish to generate the manuscript figures from our model parameters and saved ERP data, this can be done using the plotting scripts described in the next section.)
*	<span>**SUCCSTOP.txt**</span> Is the stop-signal-locked FC-ERP on successful stop trials.
*	<span>**FAILSTOP.txt**</span> Is the stop-signal-locked FC-ERP on failed stop trials.

***
## Code
This repository contains .m files (we ran with MATLAB R2017b) which reproduce the main result figures in the manuscript, using the empirical data downloaded from OSF and model parameters and output in 'HNN_sims' and 'HNN_params' directories respectively. 
 *	**getSSTERPs.m**
    *	This script epochs preprocessed EEG data from Wessel (2020), plots frontocentral ERPs, saves .txt files for ERPs imported into HNN GUI, and performs statistical analysis of peaks/onsets of FC-ERP features.
    *   Note that installation of [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php) is required to successfully execute this script. In addition, paths to the OSF data and EEGLAB should be changed in the first section to match local user paths. Support functions from /functions are necessary as well.
     
*   **plotERPs.m**
    * This script plots the grand-average, stop-signal-locked FC-ERPs for successful and failed stop trials included in figure 2C and outputs results of statistical tests of FC-ERP features. Requires you have downloaded /stats or  have run the analysis code from getSSTERPs.m to output those variables.

*	**plot_SSmodel.m**
    *	This script plots Figure 3A.
    *   HNN model for successful stop trials.
    
*   **plot_FStimingmodel.m**  
    *	This script plots Figure 3B.
    *   HNN model for failed stop trials, optimized from SSmodel parameters with no change in prox2 drive's synpatic weights.

*   **plot_FSonespikemodel.m**  
    *	This script plots Figure 4A.
    *   Alternate HNN model for failed stop trials in which prox2 was reduced from one to two spikes and optimization conducted.

*   **plot_FSstrengthmodel.m**  
    *	This script plots Figure 4B.
    *   Alternate HNN model for failed stop trials in the FS model was optimized from SSmodel parameters with changes in prox2 synaptic weights, but not timing allowed.
    
*   **plot_FSinitsomaticinhibmodel.m**  
    *	This script plots Figure 5A.
    *   Alternate HNN model for failed stop trials in which somatic inhibition resulting from intial prox1 drive was simulated.

*   **plot_FSoptsomaticinhibmodel.m**  
    *	This script plots Figure 5B.
    *   Alternate HNN model for failed stop trials in which initial somatic inhibition parameters from Fig. 5A were optimized.

***
## HNN_params
The ‘HNN_params’ directory contains one .param file per for each simulation that generated a subfigure in the paper. Each file contains all parameters required for a simulation in the Human Neocortical Neurosolver (HNN) software (but again note that the resulting dipole will look different than in the manuscript figures if one uses them to run simulations with hnn-core or without the updated calcium dynamics for HNN GUI).
*	<span>**SSmodel**</span> contains parameters for the simulation of SUCCSTOP (Figure 3A).
*	<span>**FStimingmodel**</span> contains parameters for the simulation of FAILSTOP with no change in prox2 drive's synaptic weights from the SSmodel params (Figure 3B).
*	<span>**FSonespikemodel**</span> contains parameters for the alternate simulation of FAILSTOP with a prox2 drive of one instead of two spike(s) (Figure 4A).
*	<span>**FSstrengthmodel**</span> contains parameters for the alternate simulation of FAILSTOP with no change in prox2 drive's timing/spread from the SSmodel params (Figure 4B).
*	<span>**FSinitsomaticinhibmodel**</span> contains parameters for the alternate simulation of FAILSTOP changes made to prox1 to simulate somatic inhibition (Figure 5A).
*	<span>**FSoptsomaticinhibmodel**</span> contains parameters for the alternate simulation of FAILSTOP parameters simulating somatic inhibition in prox1 optimized (Figure 5B).

***
## HNN_sims
The ‘HNN_sims’ directory contains HNN output associated with each of the parameter files in 'HNN_params'.
*	Each subdirectory corresponds to a parameter file in 'HNN_Parameters' and contains the following .txt files:  
    <sub> **dpl**	contains the averaged dipole *(column 1: time steps in ms, column 2: aggregate dipole, column 3: Layer II/III dipole, column 4: Layer V dipole)*   
     **dpl_0 - dpl_49**	contains the dipole associated with each trial *(here, 50 trials, column structure as in dpl)*   
     **spk_0**	contains the spike timing of incoming proximal and distal inputs and outcoming spikes from model basket and pyramidal cells*        </sub> 

***

Further information, code, and data may be available upon request. 
Please refer to the manuscript or contact darcy.diesburg@gmail.com. 
For further information regarding the Human Neocortical Neurosolver, or to run simulations using the parameter files provided,
please refer to [jonescompneurolab/hnn](https://github.com/jonescompneurolab/hnn) or [hnn.brown.edu](https://hnn.brown.edu/).
