Main pipeline epitopes | draft
#bin bash here

#Receives
#project_name, pdb_target, hotspots, ...?


#Returns list of design sequences that work and scatter plot, should probably return csv file with the sequences and values

mkdir $project_name/
$project_path= ./$project_name/
mkdir $project_path/1-BackboneDesign
mkdir $project_path/1.5-FilteringBackbones
mkdir $project_path/2-SequenceDesign
mkdir $project_path/3-FoldabilityTest
mkdir $project_path/results


# 1-RFDifussion step: create backbones
#Generate N backbones
sbatch 1-run_backboneDesign.sh -o $project_name [-i INPUT_PDB] [-g CONTIGMAP] -n 100


# 1.5- Filter backbones, 

# Filter 1: in this case to make the C and N terminus to be in the same side and far from the interacting region between A and B, padding and tolerance might need fine tuning 
# More filters TBA
python 1_5_filter_NCterminus.py -d /path/to/pdb_directory -p 10.0 -t 60.0


# 2-Protein MPNN: Generate sequences with Soluble MPNN chain A is always binder in RFDiffusion
sbatch 2-sequenceDesign_solubleMPNN_de_novo.sh -f $project_path/1.5-FilteringBackbones/output/ -o $project_path/2-SequenceDesign -c 'A'


# 3-FoldabilityTest
#Step 1: actually running AF2 for all sequences, this is the slowest step in my experience
sbatch 3-run_FoldabilityTest.sh -s $project_path/2-SequenceDesign -o $project_path/3-FoldabilityTest/ouput/

#Plot the results by pLDDT and RMSD, returns summary with list of design squences that have pLDDT > 95 and RMSD < 1.5 A and scatterplot png
python 3-FoldabilityTest_plot_filter.py --folder_of_folders $project_path/3-FoldabilityTest/ouput/ --reference_folder $project_path/1-BackboneDesign/ouput --output_csv "output.csv" --filtered_output_csv "filtered_output.csv" --summary_file "summary_foldabilityTest.txt" --plot_file "output_directory/pldds_vs_rmsd_plot.png" --plddt_threshold 95 --rmsd_threshold 2
