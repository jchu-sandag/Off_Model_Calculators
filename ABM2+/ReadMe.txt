1. To create the conda environment for OMC, run
conda env create -f environment_omc.yml

2. To activate created conda enviroment, run
conda activate omc

3. To run the calculators in batch format, direct to the corresponding folder for the omc_batch_abm2+.bat file and call
omc_batch_abm2+.bat
If user needs to run individual OMC, the command is:
python src/XX_calculator.py data/config.yml

**Note-
Temporarily, to allow the transit_ev omc to run, user will need to remove the extra white space after header "OP_Headway" in the input/trrt.csv file
