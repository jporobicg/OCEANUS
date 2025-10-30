#!/bin/bash
#SBATCH --job-name=oceanus_vars
#SBATCH --account=OD-234462
#SBATCH --output=oceanus_vars_%A_%a.out
#SBATCH --error=oceanus_vars_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-7


# load python
module load python

# Ensure we run from the project directory
cd /home/por07g/Documents/Code_Tools/OCEANUS

# Unbuffered Python stdout/stderr so prints appear immediately
export PYTHONUNBUFFERED=1
# Define years array
YEARS=(2017 2018 2019 2020 2021 2022 2023 2024)

# Select year from array index
YEAR=${YEARS[$SLURM_ARRAY_TASK_ID]}

echo "Starting variables processing for year $YEAR"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"

# Run variables pipeline (unbuffered)
python -u Get_variables_main.py "$YEAR"

echo "Completed variables processing for year $YEAR"
echo "End time: $(date)"


