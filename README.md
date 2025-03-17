

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

basePath/
├── behaveData/
│   ├── allSubCombined/ (contains XXXX_allBlockData.mat)
│   └── subCombined/
├── eegData/ (contains XXXX_ALP_FILT_STIM.mat)
├── eyeData/
├── figDir/

3. Download EEG data (`XXXX_ALP_FILT_STIM.mat`) into eegData.
4. Download behavioral data files (XXXX_allBlockData.mat) into behaveData/allSubCombined.
5. Place eye-tracking data files into the eyeData folder.
6. Download the sharedMatlabUtilities and place in basePath.

Functions included:
- eeg_analysisFunc
- eyeRegressionFunc
- eyeRegressionPredResp
- computeLearningRate
- getEEGCluster
- shared_variables
- chanlocs

### Task Path
- Contains MATLAB scripts implementing the behavioral and EEG/pupil recording tasks using Psychtoolbox-3.

Setup:
1. Install MATLAB (version R2023b or newer).
2. Install Psychtoolbox-3.
3. Run the MATLAB scripts to present tasks and collect data.

### Data Path
- Contains simulated data to run analysis script
