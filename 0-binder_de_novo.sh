#!/bin/bash
#SBATCH -p gpu-vladimir           # Partition name
#SBATCH -t 20-12:00:00            # 20 days 12 hours just in case
#SBATCH --job-name=binder # Job name
#SBATCH --mem=25G ####testing lowering cpu mem here (needed for filtering steps that run here and not as separate job)
#SBATCH --mail-user=ahgonzalez@ucdavis.edu # Mail me after run
#SBATCH --mail-type=END           # Mail at end of run


#############################################################################################
######## This version assumes no visual inspection for backbones, only RoG applied ########## 
#############################################################################################

# Receives project_name, pdb_target, hotspots, ...

#Example Cav2.2_AID
# Data needed for backbone generation
project_name="testPipeline_Cav22_AID"
project_path="/share/yarovlab/ahgz/Binders-Review/Cav22_AID_site/Test-Pipeline/"
pdb_target="/share/yarovlab/ahgz/Binders-Review/Cav22_AID_site/Test-1/inputs/7miy_truncated.pdb"
contig_map="[A332-406/0 A464-786/0 60-100]"
hotspots="[A381,A384,A385,A388,A389,A391,A392]"
number_models=10000

#Default ROG for potentials RFDiff
rog_value=5

# Data needed for 4-Multimer Test, could be construct with gly linkers if needed  
target_seq="KPFEIIILLTIFANCVALAIYIPFPEDDSNATNSNLERVEYLFLIIFTVEAFLKVIAYGLLFHPNAYLRNGWNLLDFIIVVVGLFSAILEQATKADGANALGGKGAGFDVKALRAFRVLRPLRLVSGVPSLQVVLNSIIKAMVPLLHIALLVLFVIIIYAIIGLELFMGKMHKTCYNQEGIADVPAEDDPSPCALETGHGRQCQNGTVCKPGWDGPKHGITNFDNFAFAMLTVFQCITMEGWTDVLYWVNDAVGRDWPWIYFVTLIIIGSFFVLNLVLGVLSGEFSKEREKAKARGDFQKLREKQQLEEDLKGYLDWITQAEDIDPENEDEGMDEEKPRNMSMPTSETESVNTENVAGGDIEGENCGARLAHRISKSKFSRYWRRWNRFCRRKCRAAVKSNVFYWLVIFLVFLNTLTIASEHYNQPNWLTEVQDTANKALLALFTAEMLLKMYSLGLQAYFVSLFNRFDCFVVCGGILETILVETKIMSPLGISVLRCVRLLRIFKITRYWNSLSNLVASLLNSVRSIASLLLLLFLFIIIFSLLGMQLFGGKFNFDEMQTRRSTFDNFPQSLLTVFQILTGEDWNSVMYDGIMAYGGPSFPGMLVCIYFIILFICGNYILLNVFLAIAVDNLADAESLTSAQKEEEEEKERKKLARTASPEKKQELVEKPAVGESKEEKIELKSITADGESPPATKINMDDLQPNENEDKSPYPNPETTGEEDEEEPEMPVGPRPRPLSELHLKEKAVPMPEASAFFIFSSNNRFRLQCHRIVNDTIFTNLILFFILLSSISLAAEDPVQHTSFRNHILFYFDIVFTTIFTIEIALKILGNADYVFTSIFTLEIILKMTAYGAFLHKGSFCRNYFNILDLLVVSVSLISFGIQSSAINVVKILRVLRVLRPLRAINRAKGLKHVVQCVFVAIRTIGNIVIVTTLLQFMFACIGVQLFKGKLYTCSDSSKQTEAECKGNYITYKDGEVDHPIIQPRSWENSKFDFDNVLAAMMALFTVSTFEGWPELLYRSIDSHTEDKGPIYNYRVEISIFFIIYIIIIAFFMMNIFVGFVIVTFQEQGEQEYKNCELDKNQRQCVEYALKARPLRRYIPKNQHQYKVWYVVNSTYFEYLMFVLILLNTICLAMQHYGQSCLFKIAMNILNMLFTGLFTVEMILKLIAFKPKGYFSDPWNVFDFLIVIGSIIDVILSETNHYFCDAWNTFDALIVVGSIVDIAITEVNPAEHTQCSPSMNAEENSRISITFFRLFRVMRLVKLLSRGEGIRTLLWTFIKSFQALPYVALLIVMLFFIYAVIGMQVFGKIALNDTTEINRNNNFQTFPQAVLLLFRCATGEAWQDIMLACMPGKKCAPESEPSNSTEGETPCGSSFAVFYFISFYMLCAFLIINLFVAVIMDN"



source /share/yarovlab/ahgz/.bashrc
conda activate base

#The latest gcc version in barbera:
module load gcc/13.2.0

# Returns list of design sequences that work and scatter plot, should probably return csv file with the sequences and values

# Set up project directory structure

#project_path="./$project_name"
mkdir -p "$project_path"/{1-BackboneDesign,1.5-FilteringBackbones,2-SequenceDesign,3-FoldabilityTest,4-MultimerTest,5-Docking,results}

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
python3 /share/yarovlab/ahgz/scripts/binderDesign/1_1_filter_RoG.py \
    --input-dir "$project_path/1-BackboneDesign/T30/" \
    --chain A \
    --output-dir "$project_path/1.5-FilteringBackbones/" \
    --rg-cutoff 11.9

# Filter 2: in this case to make the C and N terminus to be in the same side and far from the interacting region between A and B, padding and tolerance might need fine tuning


# More filters TBA

###### input for s

# Step 2: Sequence Design with Ligand MPNN
#Generate sequences with Soluble MPNN chain A is always binder in RFDiffusion
sbatch --wait /share/yarovlab/ahgz/scripts/binderDesign/2_sequenceDesign_ligandMPNN_de_novo.sh \
    -f "$project_path/1.5-FilteringBackbones/visually_inspected/" \
    -o "$project_path/2-SequenceDesign/" \
    -c 'A'


# Step 3: Foldability Test of binder only
# Part 1: actually running AF2 for all sequences, this is the slowest step in my experience
sbatch --wait /share/yarovlab/ahgz/scripts/binderDesign/3_foldabilityTest_runAF2_monomer.sh \
    -s "$project_path/2-SequenceDesign/seqs/" \
    -a "$project_name" \
    -o "$project_path/3-FoldabilityTest/af2_output/"

conda activate base
mkdir -p "$project_path"/3-FoldabilityTest/output_results/

# Part 2: Filtering and making the plot, csv file and making copy of files that passes
python3 /share/yarovlab/ahgz/scripts/binderDesign/3_foldabilityTest_filter_plot.py \
    --af-models "$project_path/3-FoldabilityTest/af2_output/" \
    --rfdiff-backbones "$project_path/1.5-FilteringBackbones/visually_inspected/" \
    --output-dir "$project_path/3-FoldabilityTest/output_results/" \
    --plddt_threshold 90\
    --rmsd_threshold 2

echo "Foldability test monomer completed for project: $project_name"

# Step 4: Multimer test
# Part 1: Submission to AF2 multimer, binder first so it will be A and target will be B
sbatch --wait /share/yarovlab/ahgz/scripts/binderDesign/4_multimerTest_runAF2.sh \
  --fasta_file "$project_path/3-FoldabilityTest/output_results/filtered_passed_seqs.fasta" \
  --output_path "$project_path/4-MultimerTest/af2_output/" \
  --target_sequence "$target_seq" \
  --template_pdb "$project_path/templates/"

#Part 2: Filtering multimer foldability test and making plot, there is also --robust version of this script
python3 /share/yarovlab/ahgz/scripts/binderDesign/4_multimerTest_filter_plot.py \
    --af-models "$project_path/4-MultimerTest/af2_output/" \
    --rfdiff-backbones "$project_path/1.5-FilteringBackbones/visually_inspected/" \
    --output-dir "$project_path/4-MultimerTest/output_results/" \
    --plddt-threshold 85\
    --rmsd-threshold 3 \

echo "Foldability test multimer filtering completed for project: $project_name"

#Step 5: Docking binder-channel for each complex that passed the multimer test
#Part 1: Submit and run all docking 1,000 poses each
sbatch --wait /share/yarovlab/ahgz/scripts/binderDesign/5-docking_submit_all.sh \
    --input-dir "$project_path/4-MultimerTest/output_results/models/"\
    --output-dir "$project_path/5-Docking/"

#Part 2: Analyze and make plots of docking results
python3 /share/yarovlab/ahgz/scripts/binderDesign/5-docking_analysis_plot.py \
   --input-dir "$project_path/5-Docking/docking_outputs/"\
   --project-name "$project_name" \
   --results-dir "$project_path/5-Docking/output_results"

#Part 3 substract lowest dG_cross and rmsd poses from silent_files 

#### testing it now!
