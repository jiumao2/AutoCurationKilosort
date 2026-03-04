# AutoCurationKilosort: Automated Curation for Kilosort Output

[![View AutoCurationKilosort on GitHub](https://img.shields.io/badge/GitHub-AutoCurationKilosort-blue.svg)](https://github.com/jiumao2/AutoCurationKilosort)

**AutoCurationKilosort** is a MATLAB-based pipeline designed to streamline and automate the curation of Kilosort outputs. By minimizing manual intervention while maintaining strict quality control, this tool helps researchers process large-scale electrophysiology datasets with greater speed, accuracy, and reproducibility.

## 💡 Why Use AutoCurationKilosort?

Manual spike sorting curation is often laborious, error-prone, and inconsistent. This tool addresses these bottlenecks by offering:
- **Time Efficiency:** Automatically discards obvious noise units, saving you from wasting time on unviable data.
- **Reduced Fatigue:** Eliminates the need to manually review thousands of units once you have a clear understanding of your dataset.
- **Automated Refinement:** Cleans up unit waveforms by removing large noise artifacts without requiring manual intervention.
- **Objective Standardization:** Replaces inconsistent manual labeling with objective, metric-driven unit classification ('good' vs. 'mua').

## ⚙️ The Pipeline

The tool processes your Kilosort output through a sequential, automated pipeline:

1. **Noise Filtering:** Removes clearly bad units based on predefined Signal-to-Noise Ratio (SNR) thresholds.
2. **Outlier Rejection:** Identifies and removes outliers within each cluster in the Principal Component (PC) feature space.
3. **Split & Merge Detection:** Highlights potential splits and merges, allowing for either fully automatic or manual decision-making.
4. **Metric Computation:** Calculates standardized quality metrics for each remaining unit.
5. **Automated Labeling:** Classifies units as `'good'` or `'mua'` (multi-unit activity) based on the computed quality metrics.
6. **Waveform Alignment:** Centers the troughs of the waveforms to the precise spike times for accurate visualization and downstream analysis.

## 🚀 Getting Started

### 1. Installation
Add the repository path to your MATLAB environment:
```matlab
addpath(genpath('path_to_AutoCurationKilosort'))
```

### 2. Configure Settings  
Edit the settings.json file to set your specific curation parameters.

### 3. Run the Curation
Copy `AutoCurationKilosort.m` to your data folder. Edit the installation paths in the script to match your local setup:

```matlab
folder_data = './catgt_Exp_g0';
setting_filenames = 'path_to_AutoCurationKilosort/settings.json';
```

Finally, run the script `AutoCurationKilosort.m` in MATLAB.

## 📖 Documentation
See [Documentation.md](Documentation.md) for details on the pipeline algorithms and advanced configuration.

## 📚 References for code

* Siegle, J. H. et al. Survey of spiking in the mouse visual system reveals functional hierarchy. Nature 592, 86–92 (2021).
* [Quality metrics](https://allensdk.readthedocs.io/en/latest/_static/examples/nb/ecephys_quality_metrics.html)
* Sean Bone (2024). [JSON+C parsing for MATLAB](https://github.com/seanbone/matlab-json-c/releases/tag/v1.1), GitHub.
* Algorithm AS 217 APPL. STATIST. (1985) Vol. 34. No.3 pg 322-325
