# Edge of Chaos Shapes Synchronization Transitions

This repository contains the MATLAB source code for the research paper:

**Edge of chaos shapes synchronization transitions and critical slowing down in coupled auditory neurons**

**Authors:** Yirui Tang, Ziyi Cheng, Xiuying Zhou, Long Jiang, Dong Yu*, Qiming Pei* *Yangtze University*

---

## 1. File Structure

Please ensure the following three MATLAB files are located in the same directory. The program relies on the interaction between the main script and these auxiliary functions.

* **`main.m`**: The main script. It initializes model parameters, iterates through simulation conditions, and controls the main execution loop.
* **`time_trial.m`**: A function module used to compute the **synchronization metric** (e.g., synchronization error) of the coupled neuron system.
* **`Lyapunov_trial.m`**: A function module used to compute the **Largest Lyapunov Exponent (LLE)**, quantifying the chaotic dynamics of the system.

> **Note:** Do not rename these files or move them into separate subfolders, as `main.m` calls the functions directly from the current working directory.

## 2. Usage Instructions

### 2.1 Prerequisites
* MATLAB (Recommended version: R2018b or later).
* Ensure the current working directory in MATLAB is set to the folder containing these files.

### 2.2 Configuration (Optional)
In `main.m`, system parameters are defined within a structure (e.g., `params`). You may modify these values to explore different dynamical behaviors:
* **A**: Amplitude of the external acoustic stimulus.
* **T**: Period of the external acoustic stimulus.
* **Coupling parameters**: Capacitance ($k$) or Josephson current intensity ($\alpha$).

### 2.3 Running the Simulation
1.  Open `main.m` in the MATLAB editor.
2.  Run the script by clicking the **Run** button or typing `main` in the Command Window.

**What the program does:**
* Iterates over the defined range of **A-T** (Amplitude-Period) parameters.
* Calls `time_trial.m` to calculate the mean synchronization error ($\theta$).
* Calls `Lyapunov_trial.m` to calculate the maximum Lyapunov exponent ($L$).
* Generates visualization data for dynamical regimes.

### 2.4 Output
Upon completion, the program will generate and display two **3D Heatmaps** in the $A-T$ parameter space:
1.  **Mean Synchronization Error Heatmap**: Visualizes the synchronization transition regions.
2.  **Maximum Lyapunov Exponent Heatmap**: Visualizes the order-chaos transition boundaries.

## 3. Notes

* **Path Settings:** Ensure MATLAB’s "Current Folder" matches the repository location; otherwise, the script will fail to locate the function files.
* **Numerical Precision:** The simulation step size and integration method are tuned for the system's dynamics. Increasing the step size excessively to speed up computation may lead to numerical instability or inaccurate results.
* **Performance:** Sweeping a high-resolution parameter space ($A-T$) can be computationally intensive. For quick tests, consider reducing the resolution or the range of parameters in `main.m`.
