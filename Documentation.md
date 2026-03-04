# 📖 Documentation: AutoCurationKilosort

This document details the pipeline, configuration options, and troubleshooting strategies for the **AutoCurationKilosort** tool. 

## ⚠️ Important Notes

* **Validation for First-Time Users:** It is highly recommended to run this tool and subsequently validate the results manually using `phy`. This is especially important for small datasets or when using the pipeline for the first time.
* **Parameter Tuning:** Default parameters were optimized using rat cortex data. You should adjust these parameters in `settings.json` based on the specific properties of your recordings (e.g., if recording from the striatum or other regions with different baseline activity).
* **Visualization (Verbose Mode):** While getting comfortable with the pipeline, set the `verbose` parameter to `true`. This saves diagnostic visualizations to a `Fig` folder, helping you verify the automated decisions. Once confident, set it to `false` to significantly speed up processing.

---

## ⚙️ Pipeline Overview

The curation process runs through the following sequential steps:

1. **Noise Filtering:** Discards units with a critically low firing rate (default: < 0.05 Hz) and low Signal-to-Noise Ratio (default: SNR < 3).
2. **Outlier Rejection:** Removes spike outliers within each cluster in the 2D Principal Component (PC) feature space. This is calculated using the Median Absolute Deviation (MAD, default threshold: 5). A valid spike should not register as an outlier in any dimension.
3. **Split & Merge Detection *(Optional)*:** Identifies potential cluster splits and merges. This step flags them for review but does not modify the data, allowing the user to finalize the split/merge manually in `phy`.
4. **Duplicate Removal:** Identifies and resolves duplicated units. Two units are flagged as duplicates if they share the same peak amplitude channel and have overlapping spike times (default: 10% of spikes occurring within 0.5 ms of each other). The unit with fewer overall spikes is removed.
5. **Quality Metrics Computation:** Calculates standardized [ECEPHYS quality metrics](https://allensdk.readthedocs.io/en/latest/_static/examples/nb/ecephys_quality_metrics.html), including: *ISI violations, Amplitude cutoffs, Presence ratio, Median Amplitude, Isolation distance, D prime, Nearest-neighbor miss rate, Nearest-neighbor hit rate,* and *L ratio*.
6. **Automated Labeling:** Classifies units as `'good'` or `'mua'` (multi-unit activity) based on user-defined thresholds. 
7. **Spike Time Realignment:** Adjusts the spike times of each cluster to perfectly center the troughs of the waveforms.

### Default Quality Criteria
| Metric | `'good'` Threshold | `'mua'` Threshold |
| :--- | :--- | :--- |
| **ISI Violations** | < 0.05 | < 0.50 |
| **Presence Ratio** | > 0.95 | > 0.90 |
| **Amplitude Cutoffs** | < 0.05 | < 0.10 |
| **NN Hit Rate** | > 0.80 | N/A |

---

## 🚀 Usage Guide

1. **Configure Settings:** Open `settings.json` to define your curation parameters (see the tables below).
2. **Setup Script:** Copy `AutoCurationKilosort.m` into your active data folder. 
3. **Update Paths:** Edit the installation and data paths inside the script to match your local directories:

```matlab
path_autocuration = 'path_to_AutoCurationKilosort\AutoCurationKilosort';
folder_data = './catgt_Exp_g0';
setting_filenames = 'path_to_settings/settings.json';
```

4. **Execute:** Run `AutoCurationKilosort.m` in MATLAB.

## 🎛️ Configuration Parameters (`settings.json`)

The `settings.json` file contains several modules you can tweak. Below are explanations of the parameters within each key module.

### `detectNoiseClusters`
| Parameter | Default | Description |
| :--- | :--- | :--- |
| `n_random_spikes` | 100 | Number of random spikes sampled to find the channel with the largest amplitude. |
| `min_firing_rate` | 0.05 | Minimum firing rate (Hz) required to keep a cluster. |
| `min_signal_to_noise_ratio` | 3 | Minimum SNR required to keep a cluster. |
| `waveform_window` | [-48, 31] | Sample window used to extract the waveform. |
| `baseline_window` | [-48, -32] | Sample window used to calculate the baseline noise variance. |

### `removeNoiseInsideCluster`
| Parameter | Default | Description |
| :--- | :--- | :--- |
| `verbose` | true | Toggles plotting of the removed spike outliers. |
| `n_pc_feature` | 2 | Number of Principal Components used for outlier detection. |
| `n_random_spikes` | 100 | Number of random spikes sampled to find the peak amplitude channel. |
| `threshold` | 5 | MAD threshold factor used to remove spikes via `isoutlier()`. |
| `waveform_window` | [-31, 32] | Sample window used to extract the waveform. |

### `duplicatedClusters`
| Parameter | Default | Description |
| :--- | :--- | :--- |
| `verbose` | true | Toggles plotting for clusters with highly overlapping spike times. |
| `dt` | 0.5 | Time window (ms) to consider two spikes as overlapping/duplicated. |
| `overlap_percentage` | 10 | Percentage of overlapping spikes required to flag a cluster as a duplicate. |
| `n_random_spikes` | 100 | Number of random spikes sampled to find the peak amplitude channel. |
| `waveform_window` | [-31, 32] | Sample window used to extract the waveform. |

### `realignSpikeTimes`
| Parameter | Default | Description |
| :--- | :--- | :--- |
| `n_random_spikes` | 100 | Number of random spikes sampled to find the peak amplitude channel. |
| `waveform_window` | [-64, 63] | Sample window used to extract the waveform. |
| `baseline_window` | [-64, -33] | Sample window used to calculate the baseline. |
| `n_channels_included` | 4 | Number of adjacent channels to include for alignment computation. |
| `verbose` | true | Toggles plotting of the before-and-after alignment. |

*(Note: Advanced PC-related metrics and explicit labeling criteria are also defined in this file under the `qualityMetrics` block).*

---

## 🛠️ Common Challenges & Solutions

### 1. Large Fraction of Noise Clusters
It is common for Kilosort to output a massive number of noise clusters, typically characterized by low firing rates and low SNRs. 
* **Solution:** The `detectNoiseClusters` function filters these out automatically. SNR is calculated using the peak-to-trough amplitude and the baseline variance:

$$
\text{SNR} = \frac{\text{amplitude}^2}{\text{Var(baseline)}}
$$

The baseline variance is calculated using the range defined in `baseline_window` (default is `[-48, -32]` samples).

### 2. Contamination by Large Noise Artifacts
Recordings can suffer from high-amplitude noise artifacts caused by poor grounding, loose connections, or animal movement. 
* **Solution:** The `removeNoiseInsideCluster` function cleans these up by projecting spikes into a PC-feature space and removing outliers using the Median Absolute Deviation (MAD) method (default threshold: 5). Any spike flagged as an outlier in any dimension is discarded.

### 3. Duplicate Units
Kilosort may detect the same unit multiple times across different channels. This presents as a shifted waveform shape and an abnormally high central peak in the cross-correlogram (distinct from a falsely split unit).
* **Solution:** The `duplicatedClusters` function resolves this by looking for substantial temporal overlap. If two units share the same peak channel and >10% of their spikes overlap within 0.5 ms, they are flagged as duplicates, and the unit with fewer spikes is removed.

### 4. Subjectivity in Quality Determination
Determining unit quality visually can be highly subjective and difficult to reproduce. 
* **Solution:** The pipeline uses strict, quantitative metrics to categorize units, ensuring consistency:
  * **`good`**: Well-isolated units. Characterized by a clear waveform, high SNR, clean ISI histogram, healthy autocorrelogram, and tight amplitude distribution.
  * **`mua`**: Multi-unit activity. The cluster contains viable neural data but is not cleanly isolated enough to be considered a single unit.
  * **`noise`**: Highly contaminated clusters (filtered out).
