#!/bin/bash
#SBATCH -p gpu-vladimir           # Partition name
#SBATCH -t 20-12:00:00            # 20 days 12 hours just in case
#SBATCH --job-name=binder # Job name
#SBATCH --mem=32G ####testing lowering cpu mem here (needed for filtering steps that run here and not as separate job)
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

# Data needed for 
target_seq="KPFEIIILLTIFANCVALAIYIPFPEDDSNATNSNLERVEYLFLIIFTVEAFLKVIAYGLLFHPNAYLRNGWNLLDFIIVVVGLFSAILEQATKADGANALGGKGAGFDVKALRAFRVLRPLRLVSGVPSLQVVLNSIIKAMVPLLHIALLVLFVIIIYAIIGLELFMGKMHKTCYNQEGIADVPAEDDPSPCALETGHGRQCQNGTVCKPGWDGPKHGITNFDNFAFAMLTVFQCITMEGWTDVLYWVNDAVGRDWPWIYFVTLIIIGSFFVLNLVLGVLSGEFSKEREKAKARGDFQKLREKQQLEEDLKGYLDWITQAEDIDPENEDEGMDEEKPRNMSMPTSETESVNTENVAGGDIEGENCGARLAHRISKSKFSRYWRRWNRFCRRKCRAAVKSNVFYWLVIFLVFLNTLTIASEHYNQPNWLTEVQDTANKALLALFTAEMLLKMYSLGLQAYFVSLFNRFDCFVVCGGILETILVETKIMSPLGISVLRCVRLLRIFKITRYWNSLSNLVASLLNSVRSIASLLLLLFLFIIIFSLLGMQLFGGKFNFDEMQTRRSTFDNFPQSLLTVFQILTGEDWNSVMYDGIMAYGGPSFPGMLVCIYFIILFICGNYILLNVFLAIAVDNLADAESLTSAQKEEEEEKERKKLARTASPEKKQELVEKPAVGESKEEKIELKSITADGESPPATKINMDDLQPNENEDKSPYPNPETTGEEDEEEPEMPVGPRPRPLSELHLKEKAVPMPEASAFFIFSSNNRFRLQCHRIVNDTIFTNLILFFILLSSISLAAEDPVQHTSFRNHILFYFDIVFTTIFTIEIALKILGNADYVFTSIFTLEIILKMTAYGAFLHKGSFCRNYFNILDLLVVSVSLISFGIQSSAINVVKILRVLRVLRPLRAINRAKGLKHVVQCVFVAIRTIGNIVIVTTLLQFMFACIGVQLFKGKLYTCSDSSKQTEAECKGNYITYKDGEVDHPIIQPRSWENSKFDFDNVLAAMMALFTVSTFEGWPELLYRSIDSHTEDKGPIYNYRVEISIFFIIYIIIIAFFMMNIFVGFVIVTFQEQGEQEYKNCELDKNQRQCVEYALKARPLRRYIPKNQHQYKVWYVVNSTYFEYLMFVLILLNTICLAMQHYGQSCLFKIAMNILNMLFTGLFTVEMILKLIAFKPKGYFSDPWNVFDFLIVIGSIIDVILSETNHYFCDAWNTFDALIVVGSIVDIAITEVNPAEHTQCSPSMNAEENSRISITFFRLFRVMRLVKLLSRGEGIRTLLWTFIKSFQALPYVALLIVMLFFIYAVIGMQVFGKIALNDTTEINRNNNFQTFPQAVLLLFRCATGEAWQDIMLACMPGKKCAPESEPSNSTEGETPCGSSFAVFYFISFYMLCAFLIINLFVAVIMDN"



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
    -d "$project_path/1.5-FilteringBackbones/output/filtered_RoG/" \
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
    -f "$project_path/1.5-FilteringBackbones/output/filtered_RoG/" \
    -o "$project_path/2-SequenceDesign/" \
    -c 'A'

echo "Sequence design completed for project: $project_name"

# Step 3: Foldability Test of binder only
# Part 1: actually running AF2 for all sequences, this is the slowest step in my experience
sbatch --wait /share/yarovlab/ahgz/scripts/binderDesign/3_foldabilityTest_runAF2_monomer.sh \
    -s "$project_path/2-SequenceDesign/seqs/" \
    -a "$project_name" \
    -o "$project_path/3-FoldabilityTest/output/"

echo "Foldability test monomer completed for project: $project_name"

conda activate base
mkdir -p "$project_path"/3-FoldabilityTest/output_results/

#Filtering monomer foldability test
python /share/yarovlab/ahgz/scripts/binderDesign/3_foldabilityTest_filter_plot_robust.py \
    --af-models "$project_path/3-FoldabilityTest/output/" \
    --rfdiff-backbones "$project_path/1.5-FilteringBackbones/output/visually_inspected/" \
    --output-dir "$project_path/3-FoldabilityTest/output_results/" \
    --plddt_threshold 95\
    --rmsd_threshold 2

echo "Foldability test filtering completed for project: $project_name"

mkdir -p "$project_path"/4-MultimerTest/
mkdir -p "$project_path"/4-MultimerTest/output

sbatch --wait /share/yarovlab/ahgz/scripts/binderDesign/4_foldabilityTest_runAF2_multimer.sh \
  --fasta_file "$project_path/3-FoldabilityTest/output_results/filtered_passed_seqs.fasta" \
  --output_path "$project_path/4-MultimerTest/output/" \
  --target_sequence $target_seq \
  --template_pdb "$project_path/templates/"

echo "Foldability test multimer completed for project: $project_name"

output_dir="${project_path}/4-MultimerTest/output_results/"
# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Run foldability test
python /share/yarovlab/ahgz/scripts/binderDesign/4_multimerTest_filter_plot_robust.py \
    --af-models "${project_path}/4-MultimerTest/output/" \
    --rfdiff-backbones "${project_path}/1.5-FilteringBackbones/output/visually_inspected" \
    --output-dir "$output_dir" \
    --plddt-threshold 90.0 \
    --rmsd-threshold 2.0

echo "Foldability test multimer filtering completed for project: $project_name"
