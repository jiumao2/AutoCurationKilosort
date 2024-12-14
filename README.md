# AutoCurationKilosort: (Semi-) Automated Curation for Kilosort Output

A tool designed to assist researchers in automating the curation process of Kilosort outputs. It helps streamline the curation of sorting results, reducing manual effort while maintaining control over the final output for increased accuracy and trustworthiness.

## How to use it

- Add the path to the repository to the MATLAB path

```MATLAB
addpath(genpath('path_to_AutoCurationKilosort'))
```

- Edit the `settings.json` file to set the curation parameters

- Copy `AutoCurationKilosort.m` to the your data folder. Edit the installation path in the script.

```MATLAB
folder_data = './catgt_Exp_g0';
setting_filenames = 'path_to_AutoCurationKilosort/settings.json';
```

- Run the script `AutoCurationKilosort.m`

## The pipeline

- Remove the clearly bad units with defined SNR thresholds
- Remove the outlies inside each cluster in the PC-feature space
- Find the potential splits and merges and make the dicisions automatically / manually
- Compute the quality metrics for each unit
- Label the units as 'good' or 'mua' based on the quality metrics
- Center the troughs of the waveforms to the spike times

## Why should we use it

- Do not waste time on the clearly bad units (noise)
- It is tiring to look at thousands of units, after you have gained the clear understanding of the data
- It is laborious to curate the waveforms of each unit by remove the large noise in the waveforms
- It is difficult to find the potential splits and merges and make the dicisions
- It is error-prone and inconsistant to label the units as 'good' or 'mua' manually




## References for code

> Sean Bone (2024). JSON+C parsing for MATLAB (https://github.com/seanbone/matlab-json-c/releases/tag/v1.1), GitHub.
> Algorithm AS 217 APPL. STATIST. (1985) Vol. 34. No.3 pg 322-325
