#!/bin/bash
#SBATCH -p gpu-vladimir           # Partition name
#SBATCH -t 10-12:00:00            # 10 days 12 hours just in case
#SBATCH --job-name=binder # Job name
#SBATCH --mem=160G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu # Mail me after run
#SBATCH --mail-type=END           # Mail at end of run

# Receives project_name, pdb_target, hotspots, ...

project_name="testPipeline_Cav22_AID"
project_path="/share/yarovlab/ahgz/Binders-Review/Cav22_AID_site/Test-Pipeline/"
pdb_target="/share/yarovlab/ahgz/Binders/Test_Pipeline/inputs/7miy_truncated.pdb"
hostpots=""
#Receives
#project_name, pdb_target, hotspots, ...?

source /share/yarovlab/ahgz/.bashrc

conda activate


#Returns list of design sequences that work and scatter plot, should probably return csv file with the sequences and values

# Set up project directory structure

#project_path="./$project_name"
mkdir -p "$project_path"/{1-BackboneDesign,1.5-FilteringBackbones,2-SequenceDesign,3-FoldabilityTest,results}

# Step 1: Backbone DesignBackbone Design
# --wait falg is neccessary in order to let the job end until starts the next step
#Generate 100 backbones

#Usage:
#sbatch --wait 1-run_backboneDesign.sh -o "$project_name" [-i INPUT_PDB] [-g CONTIGMAP] -n 100
#Example Cav2.2_AID
sbatch --wait 1-run_backboneDesign.sh \
    "$project_path/1-BackboneDesign/output/" \
    "/share/yarovlab/ahgz/Binders-Review/Cav22_AID_site/Test-1/inputs/7miy_truncated.pdb" \
    "[A332-406/0 A464-786/0 60-100]" \
    "[A381,A384,A385,A388,A389,A391,A392]" \
    100

# Step 1.5: Filtering Backbones
# In this case waits automatically
# Filter 1: in this case to make the C and N terminus to be in the same side and far from the interacting region between A and B, padding and tolerance might need fine tuning

python /share/yarovlab/ahgz/scripts/binderDesign/1_5_filter_NCterminus.py -d "$project_path/1-BackboneDesign/output/" -p 10.0 -t 60.0 -c 5.0 -o "$project_path/1.5-FilteringBackbones/output/"

# More filters TBA

# Step 2: Sequence Design with Soluble MPNN
#Generate sequences with Soluble MPNN chain A is always binder in RFDiffusion

sbatch --wait /share/yarovlab/ahgz/scripts/binderDesign/2-sequenceDesign_solubleMPNN_de_novo.sh -f "$project_path/1.5-FilteringBackbones/output/" -o "$project_path/2-SequenceDesign" -c 'A'

# Step 3: Foldability Test of binder only
# Part 1: actually running AF2 for all sequences, this is the slowest step in my experience
sbatch --wait /share/yarovlab/ahgz/scripts/binderDesign/3-run_FoldabilityTest.sh -s "$project_path/2-SequenceDesign" -o "$project_path/3-FoldabilityTest/output/"

# Part 2: Foldability Test Filtering and Plotting
# Run final Python script to generate plots and summary
#Plot the results by pLDDT and RMSD, returns summary with list of design squences that have pLDDT > 95 and RMSD < 1.5 A and scatterplot png

#python /share/yarovlab/ahgz/scripts/binderDesign/3-FoldabilityTest_plot_filter.py --folder_of_folders "$project_path/3-FoldabilityTest/output/" --reference_folder "$project_path/1-BackboneDesign/output" --output_csv "output.csv" --filtered_output_csv "filtered_output.csv" --summary_file "summary_foldabilityTest.txt" --plot_file "output_directory/pldds_vs_rmsd_plot.png" --plddt_threshold 95 --rmsd_threshold 2

#Step 4: AF multimer and Docking(?) step