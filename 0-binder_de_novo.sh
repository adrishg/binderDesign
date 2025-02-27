#!/bin/bash
#SBATCH -p gpu-vladimir           # Partition name
#SBATCH -t 20-12:00:00            # 20 days 12 hours just in case
#SBATCH --job-name=binder # Job name
#SBATCH --mem=32G ####testing lowering cpu mem here (needed for filtering steps that run here and not as separate job)
#SBATCH --mail-user=ahgonzalez@ucdavis.edu # Mail me after run
#SBATCH --mail-type=END           # Mail at end of run

# Receives project_name, pdb_target, hotspots, ...

##Test username github jeje

#Example Cav2.2_AID
project_name="testPipeline_Cav22_AID"
project_path="/share/yarovlab/ahgz/Binders-Review/Cav22_AID_site/Test-Pipeline/"
pdb_target="/share/yarovlab/ahgz/Binders-Review/Cav22_AID_site/Test-1/inputs/7miy_truncated.pdb"
contig_map="[A332-406/0 A464-786/0 60-100]"
hotspots="[A381,A384,A385,A388,A389,A391,A392]"
number_models=100

#Default ROG for potentials RFDiff
rog_value=5

source /share/yarovlab/ahgz/.bashrc
conda activate base

#The latest gcc version in barbera:
module load gcc/13.2.0

# Returns list of design sequences that work and scatter plot, should probably return csv file with the sequences and values

# Set up project directory structure

#project_path="./$project_name"
mkdir -p "$project_path"/{1-BackboneDesign,1.5-FilteringBackbones,2-SequenceDesign,3-FoldabilityTest,results}

# Step 1: Backbone DesignBackbone Design
# --wait falg is neccessary in order to let the job end until starts the next step
#Generate 100 backbones

#Usage:
#sbatch --wait 1-run_backboneDesign.sh -o "$project_name" [-i INPUT_PDB] [-g CONTIGMAP] -n 100
sbatch --wait /share/yarovlab/ahgz/scripts/binderDesign/1-run_backboneDesign.sh \
    "$project_path/1-BackboneDesign/output/" \
    "$pdb_target" \
    "$contig_map" \
    "$hotspots" \
    $number_models \
    $rog_value \

# Step 1.5: Filtering Backbones
# In this case waits automatically

#Adding activating base environment here since step 1 activates RFDiffusion environment that lacks pandas so error
conda activate base

# Filter 1: Radious of gyration (ROG) and also calculates "sphericality" but only filtering by ROG
python /share/yarovlab/ahgz/scripts/binderDesign/1_6_filter_Compactness.py \
    -d "$project_path/1-BackboneDesign/output/" \
    -c A \
    -o "$project_path/1.5-FilteringBackbones/output/" \
    --rg_cutoff 15.0

# Filter 2: in this case to make the C and N terminus to be in the same side and far from the interacting region between A and B, padding and tolerance might need fine tuning
python /share/yarovlab/ahgz/scripts/binderDesign/1_5_filter_NCterminus_de_novo.py \
    -d "$project_path/1.5-FilteringBackbones/output/filtered_compactness/" \
    -p 10.0 \
    -t 60.0 \
    -c 5.0 \
    -o "$project_path/1.5-FilteringBackbones/output/"

# More filters TBA

# Step 2: Sequence Design with Soluble MPNN
#Generate sequences with Soluble MPNN chain A is always binder in RFDiffusion

##### this part was with proteinMPNN swichting to ligandMPNN...
#sbatch --wait /share/yarovlab/ahgz/scripts/binderDesign/2_sequenceDesign_solubleMPNN_de_novo.sh \
#    -f "$project_path/1.5-FilteringBackbones/output/output_filtered/" \
#    -o "$project_path/2-SequenceDesign" \
#    -c 'A'
############## end of protein MPNN version

# Step 2 v2: Sequence Design with Ligand MPNN
sbatch --wait /share/yarovlab/ahgz/scripts/binderDesign/2_sequenceDesign_ligandMPNN_de_novo.sh \
#    -f "$project_path/1.5-FilteringBackbones/output/filtered_compactness/" \
#    -o "$project_path/2-SequenceDesign/" \
#    -c 'A'

# Step 3: Foldability Test of binder only
# Part 1: actually running AF2 for all sequences, this is the slowest step in my experience
sbatch --wait /share/yarovlab/ahgz/scripts/binderDesign/3_foldabilityTest_runAF2.sh \
    -s "$project_path/2-SequenceDesign/seqs/" \
    -a "$project_name" \
    -o "$project_path/3-FoldabilityTest/output/"

# Part 2: Foldability Test Filtering and Plotting
# Run final Python script to generate plots and summary
#Plot the results by pLDDT and RMSD, returns summary with list of design squences that have pLDDT > 95 and RMSD < 1.5 A and scatterplot png
conda activate base


#Step 3.5: Filtering by solubility here before going to AFmultimer

#python /share/yarovlab/ahgz/scripts/binderDesign/3-FoldabilityTest_plot_filter.py --folder_of_folders "$project_path/3-FoldabilityTest/output/" --reference_folder "$project_path/1-BackboneDesign/output" --output_csv "output.csv" --filtered_output_csv "filtered_output.csv" --summary_file "summary_foldabilityTest.txt" --plot_file "output_directory/pldds_vs_rmsd_plot.png" --plddt_threshold 95 --rmsd_threshold 2

#Step 4: AF multimer and Docking(?) step
#Part 1: Actually executing multimer job
#should we add -r for reference to only run multimer in those who passed the monomer foldability test?
sbatch --wait /share/yarovlab/ahgz/scripts/binderDesign/4_multimerTest.sh \
    -s "$project_path/2-SequenceDesign/seqs/" \
    -o "$project_path/3-FoldabilityTest/output/"

#Part 2: Different plot? since we cared about pae here?
conda activate base


