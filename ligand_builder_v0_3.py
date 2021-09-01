#!/usr/bin/env python
# coding: utf-8
#
# This module is created to collect different methods for ligand
# generation. *Under development*
#
#
# The data structure of the excel file should be:
#   Data in Sheet1
#   Chosen names of ligands in column A
#   Smiles stored in column B, append name (for now)
#   n lables in column A should equal n smiles in column B
#
# ----------------------------------------------------------------------
# ligand_builder v0.2, beta pre-version 1, Mary Pitman, Mobley Lab UCI
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation.
# ----------------------------------------------------------------------

__doc__="""
        Modules designed to generate ligands topologies and setup files
        for pmxworkflow via OpenBabel or Openeye. Reads data from excel
        input file.
        """
__all__ = ['collect_excel_data', 'smi_to_sdf', 'sdf_to_mol2',
           'struc_dir_setup', 'generate_yml_files']
__version__ = '0.2'
__author__ = 'M. Pitman'

import os
import sys
import subprocess as sub

import openpyxl
import yaml
from datetime import datetime
from openforcefield.topology import Molecule, Topology

# Edit for conditional import.
from openeye import oechem
from openeye.oechem import *
# ----------------------------------------------------------------------

def collect_excel_data(input_excel_file: str) -> str:
    """ Reads in 'my_data.xlsx' """
    try:
        global MY_DATA, SHEET
        MY_DATA = openpyxl.load_workbook(input_excel_file)
        SHEET = MY_DATA['Sheet1']
    except:
        print("The argument for read_excel_file is 'my_data.xlsx' .")
        sys.exit()
        
    # Test that the shape of the lists are equal.
    if len(SHEET['A']) != len(SHEET['B']):
        print("Unequal number of ligand names and smiles *.xlsx. \
              Assign a name to each smiles")
        sys.exit()
        
    # Retrieve the SMILES from the my_data.xlsx.
    global SMILES_COLLECTOR
    SMILES_COLLECTOR = ""
    for cell in SHEET['B']:
        SMILES_COLLECTOR += str('-:"') + str(cell.value) + str('" ')
          
    # Retrieve the names of each ligands from the my_data.xlsx.
    global NAME_COLLECTOR
    NAME_COLLECTOR = ""
    for cell in SHEET['A']:
        NAME_COLLECTOR += str(cell.value) + str(' ')
   
# Not totally sure the technical output type here. Need to double check
def smi_to_sdf(toolkit_to_use: str, input_excel_file: str) -> str:
    """ Outputs a .sdf for each ligand using openeye or openbabel
    example: smi_to_sdf('openbabel', 'my_data.xlsx')
    """
    collect_excel_data(input_excel_file)
    
    # Generate the .sdf container using obabel.
    if toolkit_to_use is str("openbabel"):
        import openbabel
        for cell in SHEET['B']:
            smiles_string = str('-:"') + str(cell.value) + str('" ')
            cell_shift = cell.row - 1
            lig_name = SHEET['A'][cell_shift].value
            conversion = str("obabel {} -osdf >  {}.sdf --gen3d"
                            .format(smiles_string, lig_name))
            sub.call(conversion, shell=True)
     

    # Generate the .sdf container using openeye.
    if toolkit_to_use is str("openeye"):
        for cell in SHEET['B']:
            # Shift to cell A in *.xlsx to read ligand name.
            cell_shift = cell.row - 1
            global LIG_NAME
            LIG_NAME = SHEET['A'][cell_shift].value
            smiles = str(cell.value)
            
            ifs = oemolistream()
            ifs.SetFormat(oechem.OEFormat_SMI)
            ifs.openstring(smiles)
            ofs = oemolostream('{}.sdf'.format(LIG_NAME))
            mol = OEMol()

            if (OEParseSmiles(mol, smiles) == 1):
                # Add an isomeric canonical smiles code to the bottom
                # of the .sdf file.
                OESetSDData(mol, "$SMI", OECreateIsoSmiString(mol))
                OEWriteMolecule(ofs, mol)
            else:
                sys.stderr.write("SMILES string for {} was invalid\n"
                                .format(LIG_NAME))
                sys.exit()
            ifs.close()
            ofs.close()
            
def sdf_to_mol2():
    """ Calculates the partial charges of the structural files
    and prints *.mol2 for each ligand. """
    # this can maybe be rewitting to just read the sdf in the dir so the
    for cell in SHEET['A']:
        molecule = Molecule('{}.sdf'.format(cell.value))
        #get conformers is this needed?
        #molecule.generate_conformers()
        molecule.compute_partial_charges_am1bcc()
        molecule.to_file('{}.mol2'.format(cell.value), file_format='mol2')
            
def struc_dir_setup(system_dir, date_dir, pmx_ligands_dir):
    """ Organizes the folders for *.sdf and *.mol2. Consistent with
    pmxworkflow. Iput variable system_dir is the directory for the
    simulated system. example: tyk2
    date_dir = date of the project, default 02_ligands
    pmx_ligands_dir = the folder that will store ligand info
    system_dir/date_dir/pmx_ligands_dir = pwf.ligPath
    """

    working_dir = os.getcwd()
    # Make folder structure for pmxworkflow.
    os.makedirs(f'{system_dir}/{date_dir}/{pmx_ligands_dir}', exist_ok=True)
    # If prior folders were created here, back them up.
    #A problem with this is that if you newly generated only one, it would move the other files to backup. 
    
    if len(os.listdir(f'{system_dir}/{date_dir}/{pmx_ligands_dir}')) != 0:
        os.chdir(f'{system_dir}/{date_dir}/{pmx_ligands_dir}')
        # Backup directories
        time_ran = datetime.now()
        current_time = time_ran.strftime("%H_%M_%S")
        print("Current Time =", current_time)

        # Rare but possible issue backups generates in the
        # same folder, at the precise same time. Throws error in
        # this case that "folder name" is not empty.
        # A problem with this
        for existing_dir in os.listdir():
            bu_name = "backup_{}_{}".format(current_time, existing_dir)
            os.rename(existing_dir, bu_name)
        os.chdir(f'{working_dir}')
        print(os.getcwd())
        print("Prior folders in {}/{}/{} backed up to backup_{}_*"
              .format(system_dir, date_dir, pmx_ligands_dir, current_time))
    
    # Test if no structural files found in the working directory.
    # Rewrite to pass in variable sdf or mol2 instead of duplicating.
    made_sdf_files = [f for f in os.listdir(f'{working_dir}') if f.endswith('.sdf')]
    if len(made_sdf_files) == 0:
        print("Warning: no *.sdf files found for ligands in current working directory, {}.".format(working_dir))
    made_mol2_files = [f for f in os.listdir(f'{working_dir}') if f.endswith('.mol2')]
    if len(made_mol2_files) == 0:
        print("Warning: no *.mol2 files found for ligands in current working directory, {}.".format(working_dir))
        
    # Make and move structural file directories.
    for cell in SHEET['A']:
        os.makedirs(f'{system_dir}/{date_dir}/{pmx_ligands_dir}/{cell.value}/crd', exist_ok=True)
        # If the .sdf and .mol2 files exist, move them to created directories
        if os.path.exists(f'{working_dir}/{cell.value}.sdf') == True:
            os.rename(f'{working_dir}/{cell.value}.sdf',
                      f'{system_dir}/{date_dir}/{pmx_ligands_dir}/{cell.value}/crd/{cell.value}.sdf')
            print("File {}.sdf deposited at {}/{}/{}/{}/crd"
                  .format(cell.value, system_dir, date_dir, pmx_ligands_dir, cell.value))
        if os.path.exists(f'{working_dir}/{cell.value}.mol2') == True:
            os.rename(f'{working_dir}/{cell.value}.mol2',
                      f'{system_dir}/{date_dir}/{pmx_ligands_dir}/{cell.value}/crd/{cell.value}.mol2')
            print("File {}.mol2 deposited at {}/{}/{}/{}/crd"
                  .format(cell.value, system_dir, date_dir, pmx_ligands_dir, cell.value))
    print("When you run the pmxworkflow set the pwf path input argument to {}".format(system_dir))
        
        
        
        
def generate_yml_files():

    """ Creates *.yml files for pmxworkflow RBFE calculations
    and deposit them to the set locations for pmxworkflow. """
    
"""
    ligands = {};
    for cellA, cellB in zip(SHEET['A'], SHEET['B']):
   #     ligands['---']
        ligands['name'] = '     ' + str(cellA.value),
        ligands['smiles'] = '   ' + str(cellB.value),
        ligands['docked'] = '   03_docked/{}/{}.sdf'.format(cellA.value, cellA.value)

        ligands.append(ligands)

    with open('ligands.yml', 'w') as outfile:
        yaml.dump(ligands, outfile)
       

Earlier version
  
    ligands = {}
    for cellA, cellB in zip(SHEET['A'], SHEET['B']):
        lig_entry = {
            '---'
            'name': '     ' + str(cellA.value),
            'smiles': '   ' + str(cellB.value),
            'docked': '   03_docked/{}/{}.sdf'.format(cellA.value, cellA.value),
            'measurement': 'N/A ',
            '    ki': 'N/A ',
            '    doi': 'N/A ',
            '    comment': 'N/A ',
            }

        ligands.append(lig_entry)

#    with open(r'.yml', 'w') as file:
    documents = yaml.dump(ligands, file)
  """
  
 # ^the above hasn't yet been successful. I will need to figure out how to best generate the yml files in a scriptable way. It is a little lower priority
    
#IF you want to just call this script as one thing and have it blindly do its job
# insert the if main function and just string everything together
    



        













