#!/bin/bash
#SBATCH --job-name=oceanus_transport
#SBATCH --account=OD-234462
#SBATCH --output=oceanus_transport_%A_%a.out
#SBATCH --error=oceanus_transport_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --array=0-15

# load python
module load python
# Define years array
YEARS=(2017 2018 2019 2020 2021 2022 2023 2024)

# Map array index to year (two tasks per year)
IDX=${SLURM_ARRAY_TASK_ID}
YEAR_INDEX=$(( IDX / 2 ))
YEAR=${YEARS[$YEAR_INDEX]}

echo "Starting OCEANUS processing"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Year: $YEAR"
echo "Mode: $( [ $((IDX%2)) -eq 0 ] && echo physics || echo variables )"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"

if (( IDX % 2 == 0 )); then
  # Physics/transport run
  python Get_physic_main.py $YEAR
else
  # Variables run
  python Get_variables_main.py $YEAR
fi

echo "Processing completed for year $YEAR"
echo "End time: $(date)"
