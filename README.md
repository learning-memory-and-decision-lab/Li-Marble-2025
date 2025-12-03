

## System Requirements
- **MATLAB (version R2023b)**
- Compatible with standard Windows/Mac/Linux systems running MATLAB R2021b, R2023b or newer.
- Dependencies: Psychtoolbox-3 package.

## Repository Structure

### Analysis Path
This folder includes MATLAB scripts and functions for analyzing EEG, eye-tracking, and behavioral data:
- eeg_analysisFunc
- eyeRegressionFunc
- eyeRegressionPredResp
- computeLearningRate
- getEEGCluster
- chanlocs and shared_variables (support files)

#### Setup Instructions (Analysis):
1. If accessing full dataset, prepare approximately 70 GB of space.
2. Create your directory structure:

```
basePath/
├── analysis/
├── behaveData/
│   ├── allSubCombined/
│   │   └── XXXX_allBlockData.mat
│   └── subCombined/
│       └── XXXX_3and4BlockData.mat
├── eegData/
│   └── XXXX_ALP_FILT_STIM.mat
├── eyeData/
│   └── XXXX.mat
└── figDir/ 
```
Downloading example data should take 1-2 minutes
3. Download EEG data (`XXXX_ALP_FILT_STIM.mat`) into eegData.
4. Download behavioral data files (XXXX_allBlockData.mat) into behaveData/allSubCombined.
5. Place eye-tracking data files (XXXX.mat) into the eyeData folder.
6. Download the subFunctions and helperFunctions folders and place in basePath.

Functions included:
- eeg_analysisFunc
- eyeRegressionFunc
- eyeRegressionPredResp
- computeLearningRate
- getEEG_clusterSize
- shared_variables
- chanlocs
- modelCP
- modelOB

Running the script on data included in this repository takes ~20 minutes. It should generate all parts of the figures that can be generated without significant EEG or pupil results. For the parts of the figures which require significant EEG or pupil results, it will generate results using a template frontal positive ERP and pupil dilation. There are only 5 trials in the sample EEG data, so a lot of the analyses won't work correctly.

### Task Path
- Contains MATLAB scripts implementing the behavioral and EEG/pupil recording tasks using Psychtoolbox-3.

Setup:
1. Install MATLAB (version R2023b or newer).
2. Install Psychtoolbox-3.
3. Run the MATLAB scripts to present tasks and collect data.

### Data Path
- Contains simulated data to run analysis script

### License Information
Copyright Harrison Marble and Tiantian Li, Brown University, 2025.
Licensed under the Academic Free License version 3.0.
