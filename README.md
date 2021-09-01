# ligand_builder
automation of setting up systems for the pmxworkflow

Documentation on how to build ligand systems for pmxworkflow using
ligand_builder_v_0_3.py (version 0.3)
08_31_2021

--------------------------------------------------------------------------------
Getting started
--------------------------------------------------------------------------------

Deposit the names of ligands and smiles into an excel spreadsheet. Where smiles codes are of the form ‘smile code + chosen_ligand_name’. 
	- Your chosen ligand names should be vertical is column A of *.xlsx. 
	- Smiles should be vertical in column B of *.xlsx.

Example: 

Column A			Column B
Ligand1 name		smile code 1 Ligand1 name
Ligand2 name		smile code 2 Ligand2 name
…			…
Ligandn name		smile code n Ligandn name


For the current setup, to run the pmxworkflow, you will need to first prepare a few input files. First, the macromolecule of interest. Also, four '*.yml' files:
	ligands.yml
	targets.yml
	target.yml
  edges.yml

In development: *.yml generations and deposition to appropriate folders. 

For now, if you are running this code in the folder of your system ‘tyk2’ (system_dir). Then, the '*.yml' files should be located at:

	system_dir
	|   targets.yml
	\-- date_dir
	|   \-- 00_data
	|   |	|   target.yml
	|   |	|   ligands.yml
	|   |	|   edges.yml
	
	
date_dir is a directory name of your choosing, but the date is a logical choice for bookkeeping purposes. Example: 2021_01_01

In version 0.3, system_dir and date_dir and input variables for func in ‘ligand_builder_v*.py’ so make note of them. 

--------------------------------------------------------------------------------


‘ligand_builder_v_0_3.py’ generates '*.sdf' and '*.mol2' files with partial charges and deposit them into appropriate folders for the pmxworkflow for you. 

Import as:

from ligand_builder_v0_3 import *


Functions:

collect_excel_data(input_excel_file: str) -> str
	Reads in 'my_data.xlsx'.
	example: collect_excel_data('my_data.xlsx')

	where:
	'my_data.xlsx' = the location and name of the
			excel file with your smiles and
			names.

smi_to_sdf(toolkit_to_use: str, input_excel_file: str) -> str
	Outputs a .sdf for each ligand using openeye or openbabel

	where:
	‘toolkit_name’ = ‘openbabel’ or ‘openeye’
	'my_data.xlsx' = the location and name of the
			excel file with your smiles and
			names.


sdf_to_mol2()
	Calculates the partial charges of the structural files
    	and prints *.mol2 for each ligand.


struc_dir_setup(system_dir, date_dir, pmx_ligands_dir)
	
	Organizes the folders for *.sdf and *.mol2. Consistent with
	pmxworkflow. Iput variable system_dir is the directory for the
	simulated system. example: tyk2

	where:
	system_dir = the system you are studying, example tyk2
	date_dir = date of the project
	pmx_ligands_dir = the folder that will store ligand folders,
			 default 02_ligands

	For pmxworkflow: 
	system_dir/date_dir/pmx_ligands_dir = pwf.ligPath
	
	struc_dir_setup will also output the path variable you should 
	set for the pwf variable of pmxworkflow (may be in run.py)
    
	

generate_yml_files()
	under development


--------------------------------------------------------------------------------
Ligand_builder completes the file generation required for the pmxworkflow step:
	workflow1_parameterize_ligands
  ^ which generates input files for gromacs. 

After setting up the ligands, create docked structures.
