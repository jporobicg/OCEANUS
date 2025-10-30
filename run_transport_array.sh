#!/bin/bash
#SBATCH --job-name=oceanus_transport
#SBATCH --account=OD-234462
#SBATCH --output=oceanus_transport_%A_%a.out
#SBATCH --error=oceanus_transport_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --array=0-7

# Load Python from modules (cluster environment)
module load python

# Years to process (2017–2024)
YEARS=(2017 2018 2019 2020 2021 2022 2023 2024)

# Select year by SLURM array index
YEAR=${YEARS[$SLURM_ARRAY_TASK_ID]}

echo "Starting transport (physics) processing for year $YEAR"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"

# Run physics/transport pipeline
python Get_physic_main.py "$YEAR"

echo "Completed transport processing for year $YEAR"
echo "End time: $(date)"


