# chaos-shapes-synchronization-transitions
Program User Guide
1. File Structure
Please ensure that the following three MATLAB files are located in the same folder:
1.main.m — Main script for setting model parameters and running simulations.
2.time_trial.m— Function to compute the synchronization metric of the coupled neuron system.
3.Lyapunov_trial.m— Function to compute the maximum Lyapunov exponent of the system.
Note: All three files must be present in the same folder, and the MATLAB current working directory should be set to this folder.
2. Usage Instructions
2.1 Open the Program
1. Launch MATLAB.
2. Set the current working directory to the folder containing the program files.
3. Verify that main.m, time_trial.m, and Lyapunov_trial.m are accessible in the path.
2.2 Set Parameters
In main.m, parameters are configured using a structure (e.g., params). Typical parameters include:External drive amplitude and frequency (A, T) Coupling capacitance or other system parameters
2.3 Run the Program
Run the main script from the MATLAB command window:
The program will perform the following:
Iterate over all combinations of A-T parameters.
Call time_trial.m to calculate the average synchronization error.
Call Lyapunov_trial.m to calculate the maximum Lyapunov exponent for each parameter set.
Generate two 3D heatmaps in the A-T parameter space:1. mean synchronization error 2.Maximum Lyapunov exponent
2.4 Output
After execution, two 3D heatmaps will be displayed:
1.mean Synchronization Error Heatmap
2. Maximum Lyapunov Exponent Heatmap
3. Notes
Ensure MATLAB’s current directory is correctly set; otherwise, the main script will not be able to call the auxiliary functions.
Parameter settings must be consistent with the system’s dynamics. Using excessively large step sizes may reduce accuracy.
For long computation times, consider reducing the A-T parameter ranges or using larger step sizes.
