#!/bin/bash
#SBATCH -p gpu-vladimir     # Partition name
#SBATCH --gres=gpu:1        # Request 2 GPUs
#SBATCH -t 10-2:00:00       # 1 day, adjust as needed
#SBATCH --job-name=backboneDesign # Job name
#SBATCH --mail-user=ahgonzalez@ucdavis.edu # Email after run
#SBATCH --mail-type=END

# Load environment
source /share/yarovlab/ahgz/.bashrc
conda activate SE3nv

# Command-line arguments
output_prefix=$1
input_pdb=$2
contigmap=$3
hotspot_res=$4
num_designs=${5:-10}  # Default to 10 if not provided

# Run RFdiffusion script with parameters
/share/yarovlab/ahgz/apps/RFdiffusion/scripts/run_inference.py \
    inference.output_prefix="$output_prefix" \
    inference.input_pdb="$input_pdb" \
    contigmap.contigs="$contigmap" \
    ppi.hotspot_res="$hotspot_res" \
    inference.num_designs="$num_designs" \
    denoiser.noise_scale_ca=0 \
    denoiser.noise_scale_frame=0 \
   'potentials.guiding_potentials=["type:binder_ROG,weight:100,min_dist:10"]' \
    potentials.guide_scale=2 potentials.guide_decay="quadratic"