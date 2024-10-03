Main pipeline epitopes | draft
#bin bash here


epitope_PDB= epitopeName.pdb


epitope_name = 'aliasExample'
epitope_sequence = ./share/yarovlab/ahgz/scripts/nanobodies/pdb2fasta.sh $epitope_PDB

mkdir $epitope_Name/
mkdir $epitope_name/1-BackboneDesign
mkdir $epitope_name/2-SequenceDesign
mkdir $epitope_name/3-FoldabilityTest


# 1-RFDifussion step: create backbones
#Generate N backbones
sbatch 1-run_backboneDesign.sh [-o OUTPUT_PREFIX] [-i INPUT_PDB] [-g CONTIGMAP] [-n NUM_DESIGNS]


# 1.5- Filter backbones

python analyze_pdb.py -d /path/to/pdb_directory -e LLLQKKYKDVESSLK -p 10.0 -t 60.0


# 2-Protein MPNN: Generate sequences 

sbatch 2-sequenceDesign_PoteinMPNN_soluble.sh #[-f folder_with_pdbs] [-o output_dir] [-c chains_to_design] [-e epitope_sequence]


# 3-FoldabilityTest
sbatch 3-run_FoldabilityTest.sh -s "/path/to/seq_folder" -o "/path/to/output_dir_base" -e "/path/to/conda_env" -m "your_email@example.com"


python 3-FoldabilityTest_plot_filter.py --folder_of_folders "/path/to/folder_of_folders" --reference_folder "/path/to/reference_folder" --output_csv "output.csv" --filtered_output_csv "filtered_output.csv" --summary_file "summary_foldabilityTest.txt" --plot_file "output_directory/pldds_vs_rmsd_plot.png" --plddt_threshold 95 --rmsd_threshold 2
