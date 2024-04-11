* If running process for the first time, run all steps to create the conda environment. If repeating process on the same machine, simply run step 1, 3, and 6

** Note: The process  developed based on conda version 24.3.0. The command for other versions might be different

____________________________________________________________________________________________________
1. Open the Command window by typing 'cmd' in the folder address bar in the file explorer, where
calculator files are saved. Press Enter to open the command window. Confirm that the command prompt window shows the correct folder directory.
 
____________________________________________________________________________________________________
2. Create the conda environment for OMC, run:
> conda create -n env_omc

The prompt window will show the directory where the environment will be saved. Type “y” and press Enter to proceed.
 
____________________________________________________________________________________________________
3. Activate conda environment, run:
> conda activate env_omc

When activated the name of the environment will appear before the command line.
 
____________________________________________________________________________________________________
4. Install required packages to conda environment, run: 
> conda install --file=env_packages.txt

The command window shows a list of packages to be installed. Type “y” and press Enter to proceed.
 
____________________________________________________________________________________________________
5. Install "openmatrix" package to conda environment from conda-forge channel, run: 
> conda install conda-forge::openmatrix

The command window shows a list of packages to be installed. Type “y” and press Enter to proceed.
____________________________________________________________________________________________________
6. To run the calculators in batch format simply call the "omc_batch-abm3.bat" file.
The command will show the input config file location for each calculator, runs all calculators, and saves results to the output folder.

