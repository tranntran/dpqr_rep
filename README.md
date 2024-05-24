# Differentially Private Quantile Regression Replication Code

## Replicate results from scratch

- Go to folder `src` for source codes for the simulation studies.
- Before running the code, make sure to create an output folder with the same subdirectories as this repo or modify the output directory in the code.
- Codes are numbered based on the corresponding simulation studies (i.e., files starting with 03 are related to simulation study 3).
- Codes marked with the same step (i.e., 01a) can be run simultaneously.
- Step `b` should be run only after all steps a have finished running.
- Sample Slurm script with specifications for the PSU HPC cluster can be found in the main folder.

## Replicate from available outputs

- As the code can take some time to run, we also include output from our simulations in the form of .Rdata to help with the replication process.
- Output is separated based on simulation study in the `output` folder.
- Please refer to files with step `b` in the `src` folder to see how to process output files.
- Note that the output files contain both accuracy and runtime results.

Feel free to reach out if you have any issues with the code.
