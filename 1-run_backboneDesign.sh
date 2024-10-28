#!/bin/bash
#SBATCH -p gpu-vladimir     # Partition name
#SBATCH --gres=gpu:4        # Request 1 GPUs
#SBATCH -t 1-12:00:00       # 1 day just in case
#SBATCH --job-name=RFDiffusion #Job name
#SBATCH --mem=240G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu #mail me after run
#SBATCH --mail-type=END #mail at end of run
source /share/yarovlab/ahgz/.bashrc
conda activate SE3nv
/share/yarovlab/ahgz/apps/RFdiffusion/scripts/run_inference.py inference.output_prefix=/share/yarovlab/ahgz/Cav1.2/Binders/CaS_site/Test-1/output/ inference.input_pdb=/share/yarovlab/ahgz/Cav1.2/Binders/CaS_site/Test-1/inputs/8we7_truncated.pdb 'contigmap.contigs=[A508-705/0 A1035-1197/0 A1411-1484/0 A1493-1529/0 70-100]' 'ppi.hotspot_res=[A1110, A1112, A1115, A1116, A1126, A1123]' inference.num_designs=10 denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0