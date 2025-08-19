This repository contains a single Python script that:
- Loads DSA projection data (`.raw`) and region masks (`.tif`)
- Computes inlet and aneurysm time–density curves (TDCs)
- Derives key angiographic parameters (BAT, PH, TTP, AUC, MTT, max dF/dt)
- Performs SVD-based deconvolution with Tikhonov regularization to estimate the impulse response function (IRF), and computes RBF, RBV, and MTT
- Compares reconstructed TDC (inlet ⊗ IRF) with measured aneurysm TDC and reports RMS error
- Visualizes masks overlay, TDCs, the convolved curve, and the IRF

> **Note**: Paths and filenames are hard-coded in the script. Update them to match your environment (see **Configure paths** below). The script itself should not be modified per your workflow requirement.

---
