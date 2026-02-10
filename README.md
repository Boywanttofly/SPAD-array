# SPAD-array
The experimental data and Matlab code for validation the model of the output PDF of SPAD array 
# SPAD Array Output PDF Model: Experimental Data and MATLAB Code

This repository contains the experimental data and MATLAB code for validating the analytical model of the output Probability Density Function (PDF) for a SPAD (Single-Photon Avalanche Diode) array.

## File Description

*   **`SPAD_Model_Numerical_Solving.m`**
    *   Obtains the PDF of the actual output voltage of the SPAD array.
    *   Calculates the predicted PDF using the analytical model.
    *   Computes the Kullback-Leibler (KL) divergence between experimental and model-predicted PDFs.
    *   Plots the comparative results.
    *   **Usage**: Modify the file paths inside the script to run directly.

*   **`SPAD_output_wave_display.m`**
    *   Displays the actual avalanche waveforms of the SPAD array.
    *   Performs parameter estimation from the waveforms.
    *   **Usage**: Modify the file paths inside the script to run directly.

*   **`crosstalk_model_vpp_local.m`**
    *   Displays the distribution of the maximum value (Vpp) of output pulses under weak light.
    *   Validates the geometric decay crosstalk model.
    *   Estimates the probability of a single crosstalk event.

*   **`createFit_for_exp_data_DC.m`**
    *   Distribution fitting script for the experimental dark count data.

*   **`noise_fit.m`**
    *   Function file for estimating the parameters of Gaussian electronic noise.

*   **`st_vpp_20250525_n.mat`**
    *   Experimental data file. Records the peak voltage locations of output pulses from a SPAD array under extremely weak light. Used for analyzing avalanche pulse and crosstalk distributions.

## Data Availability

**Note**: The complete experimental datasets are too large to host directly on GitHub (a single `.mat` file is approximately 190 MB).

**`st_vpp_20250525_n.mat`** contains only the extracted peak voltage locations from one experiment and is **not sufficient to run the main validation scripts independently**.

To request access to the **full experimental dataset** required to reproduce the findings in our paper, please contact: [xuwli@mail.ustc.edu.cn](mailto:xuwli@mail.ustc.edu.cn).
