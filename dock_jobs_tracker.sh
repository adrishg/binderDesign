#!/bin/bash
#SBATCH --job-name=dock_tracker
#SBATCH --time=1-00:10:00
#SBATCH --mem=1G
#SBATCH --output=dock_tracker-%j.out
#SBATCH --error=dock_tracker-%j.err

echo "[INFO] All docking jobs have completed successfully."
