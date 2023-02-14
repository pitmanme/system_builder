#!/usr/bin/env python
# coding: utf-8
#
# This module is created to generate ligands primarily for the
# pmxworkflow but can be applied more generally.
#
# Takes an excel file with smiles and ligand names as input.
# The data structure of the input excel file should be:
#   Data in Sheet1
#   Chosen names of ligands in column A
#   Smiles stored in column B, append name (for now)
#   n lables in column A should equal n smiles in column B
#
# ----------------------------------------------------------------------
# ligand_builder, version 1.0, Mary Pitman, Mobley Lab UCI
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appears in supporting documentation.
# ----------------------------------------------------------------------

__doc__="""
        Modules designed to generate ligands topologies and setup files
        for pmxworkflow via Openeye. Reads data from excel
        input file.
        """
__all__ = ['collect_excel_data', 'smi_to_sdf', 'RMSD_comparison',
           'simple_sdf_to_mol2','rigid_overlay','struc_dir_setup',
           'sdf_to_mol2_multiconf','elf_pc_gen', 'ligand_reshaping'
           ]
__version__ = '1.0'
__author__ = 'M. Pitman'

import os
import sys
import subprocess as sub

import openpyxl
from datetime import datetime
from openff.toolkit.topology import Molecule, Topology
from openeye.oechem import *
from openeye.oeomega import *

import ligpose
from build_receptor import *
from oescripts import *

# ----------------------------------------------------------------------
# CHANGES NEEDED: CLEAN UP INITIAL TWO FUNCS THAT ARE LARGELY DEPRECATED
# NOW


def collect_excel_data(input_excel_file):
    """ Reads in 'my_data.xlsx' """
    try:
        global MY_DATA, SHEET
        MY_DATA = openpyxl.load_workbook(input_excel_file)
        SHEET = MY_DATA['Sheet1']
    except:
        print("The argument for collect_excel_file is 'my_data.xlsx'.")
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
    print("Smiles and ligand names imported.")
    
    
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
            conversion = str("obabel {} -osdf >  pre_{}.sdf --gen3d"
                             .format(smiles_string, lig_name))
            sub.call(conversion, shell=True)
        print("Initial structural files pre_*.sdf generated for each ligand.")
        
    # Generate the .sdf container using openeye. This structure is not aligned.
    if toolkit_to_use is str("openeye"):
        for cell in SHEET['B']:
            # Shift to cell A in *.xlsx to read ligand name.
            cell_shift = cell.row - 1
            global LIG_NAME
            LIG_NAME = SHEET['A'][cell_shift].value
            smiles = str(cell.value)
            
            #ifs = oemolistream()
            #ifs.SetFormat(OEFormat_ISM)
            #ifs.openstring(smiles)
            ofs = oemolostream('pre_{}.sdf'.format(LIG_NAME))
            mol = oechem.OEGraphMol()
            #oechem.OESmilesToMol(mol, smiles)
            '''
            for atom in mol.GetAtoms():
                if atom.IsChiral() and not atom.HasStereoSpecified(oechem.OEAtomStereo_Tetrahedral):
                    print(1)
                    v = []
                    for neigh in atom.GetAtoms():
                        v.append(neigh)
                    atom.SetStereo(v, oechem.OEAtomStereo_Tetra, oechem.OEAtomStereo_Left)
                    
                    stereovalue = atom.GetStereo(v, OEAtomStereo_Tetrahedral)

                    print("Atom:",atom.GetIdx()," ")

                    if stereovalue == OEAtomStereo_RightHanded:
                        print("Right Handed")
                    elif stereovalue == OEAtomStereo_LeftHanded:
                        print("Left Handed")
            '''
            if (OEParseSmiles(mol, smiles) == 1):
                # Add an isomeric canonical smiles code to the bottom
                # of the .sdf file.
                oechem.OEAssignAromaticFlags(mol)
                OESetSDData(mol, "$SMI", OECreateIsoSmiString(mol))
                OEWriteMolecule(ofs, mol)
            else:
                sys.stderr.write("SMILES string for {} was invalid\n"
                                .format(LIG_NAME))
                sys.exit()
            #ifs.close()
            ofs.close()
        print("Initial structural files pre_*.sdf generated for each ligand.")

# Need to add an option for when stereochemistry is not known
# Enumerate -> create multiple?
def elf_pc_gen(smiles):
    ''' Outputs a charged molecule using am1bcc method
    and ELF conformer selection:
    (NEW) Reads the smile code
    (NEW) writes out temp sdf file
    (1) expands conformations
    (2) performs ELF conformer selection
    (3) assigns partial charges to each ELF conformer
    (4) outputs averaged charged molecule.
    '''
    
        
    # Convert to Openff molecule and generate conformations.
    off_mol = Molecule.from_smiles(smiles)

    # Temporary for testing.
    n = off_mol.name
    off_mol.to_file('step1_off_{}.sdf'.format(n), file_format = "SDF")
    
    #off_mol = Molecule.from_openeye(molecule)
    Molecule.assign_fractional_bond_orders(off_mol, bond_order_model= "am1-wiberg")
    
    Molecule.generate_conformers(off_mol, n_conformers = 100)
    #Molecule.enumerate_stereoisomers(off_mol)
    Molecule.generate_conformers(off_mol, n_conformers = 100)
    print('Expanding molecule to {} conformations.'
          .format(off_mol.n_conformers))
    Molecule.apply_elf_conformer_selection(off_mol, limit = 10)
    
    # Temporary for testing.
    off_mol.to_file('step2_off_{}.sdf'.format(n), file_format = "SDF")
    
    print('Elf conformer selection chose {} conformations.'
          .format(off_mol.n_conformers))
    nconformers = off_mol.n_conformers
        
    # Assign atom names to oe_mol if needed.
    oe_mol = Molecule.to_openeye(off_mol)
    assign_names = False
    for atom in oe_mol.GetAtoms():
        if atom.GetName()=='':
            assign_names = True
    if assign_names:
        OETriposAtomNames(oe_mol)
        
    # Set up storage for partial charges.
    partial_charges = dict()
    for atom in oe_mol.GetAtoms():
        name = atom.GetName()
        partial_charges[name] = 0.0

    # Iterate over conformations and calculate partial charges.
    conf_i = 0
    print('Calculating partial charges over ELF conformations...')
    for conformation in oe_mol.GetConfs():
        conf_i += 1
        charged_molecule = OEMol(conformation)
        off_cmol = Molecule.from_openeye(charged_molecule)
        off_cmol.compute_partial_charges_am1bcc(use_conformers=off_mol.conformers
                                                [conf_i-1:conf_i])
        charged_molecule = Molecule.to_openeye(off_cmol)
           
        # Sum partial charges.
        for atom in charged_molecule.GetAtoms():
            name = atom.GetName()
            partial_charges[name] += atom.GetPartialCharge()
           
    # Calculate average partial charge and output charged_molecule.
    charged_molecule = OEMol(oe_mol.GetConf(OEHasConfIdx(0)))
    for atom in charged_molecule.GetAtoms():
        name = atom.GetName()
        atom.SetPartialCharge(partial_charges[name] / nconformers)
    
    return charged_molecule


# Add in the new way to calc partial charges
def sdf_to_mol2_multiconf():
    """ Calculates the partial charges of the structural files
    and prints *.mol2 for each ligand. """
    for cell in SHEET['A']:
        if (os.path.isfile('pre_{}.sdf'.format(cell.value)) == True):
            #infile = '/Users/mpitman/work/system_tests/Naf/stereo_test/subtest/step1_off_stereo1.sdf'
            infile = 'pre_{}.sdf'.format(cell.value)
            ifs = oemolistream(infile)
            outfile = 'pre_{}.mol2'.format(cell.value)
            ofs = oemolostream(outfile)
            mol = OEGraphMol()
            OEReadMolecule(ifs, mol)
            charged_molecule = elf_pc_gen(mol)
            OEWriteMolecule(ofs, charged_molecule)
        ifs.close()
        ofs.close()


def single_sdf_to_charged_mol(cell):
    """ Calculates the partial charges of the structural files
    and prints *.mol2 for each ligand. """
    for cell in SHEET['B']:
        # Shift to cell A in *.xlsx to read ligand name.
        cell_shift = cell.row - 1
        global LIG_NAME
        LIG_NAME = SHEET['A'][cell_shift].value
        smiles = str(cell.value)
        charged_molecule = elf_pc_gen(mol, smiles)
    return charged_molecule
 
 
def simple_sdf_to_mol2():
    """ Calculates the partial charges of the structural files
    and prints *.mol2 using basic charge generation """
    for cell in SHEET['A']:
        molecule = Molecule('{}.sdf'.format(cell.value))
        molecule.compute_partial_charges_am1bcc()
        molecule.to_file('{}.mol2'.format(cell.value),
                         file_format='mol2')
 
# Implement this later to do a test.
def RMSD_comparison(ref, fitted):
    pair = OEMatch()
    for i_ref, i_fitted in zip(ref.GetAtoms(), fitted.GetAtoms()):
        match.AddPair(i_ref, i_fitted)
        OERMSD(ref, fitted, pair)
        overlay = True
        x = OERMSD(ref, fitted, pair, overlay)
        print(x)


def ligand_reshaping_prior():
    """ Calculates the partial charges based on expanded
    conformations and then outputs the corrected ligand
    locations in *.sdf and *.mol2.
    """
    receptor = create_bound_receptor('/Users/mpitman/work/system_tests/7cmd/structure_files/uncapped.pdb')
    for cell in SHEET['A']:
        #generate the charged molecule from expanded conformations, am1bcc
        charged_molecule = single_sdf_to_charged_mol(cell)
        #and align the charged_molecule
        pose_molecule_from_charged(receptor, cell, charged_molecule,
                                   n_conformations=25, n_poses=1)
        print('Ligand {} aligned onto receptor bound ligand'
              .format(cell.value))
              
def ligand_docking(pdb_path, coords):
    """ Calculates the partial charges based on expanded
    conformations and then outputs the corrected ligand
    locations in *.sdf and *.mol2.
    """
    receptor = docking_receptor(pdb_path, coords)
    for cell in SHEET['A']:
        #generate the charged molecule from expanded conformations, am1bcc
        charged_molecule = single_sdf_to_charged_mol(cell)
        #and align the charged_molecule
        pose_molecule_from_charged(receptor, cell, charged_molecule,
                                   n_conformations=25, n_poses=1)
        print('Ligand {} aligned onto receptor bound ligand'
              .format(cell.value))

    
#def ligand_reshaping(pdb_path, mtz_path, method):
def ligand_reshaping(pdb_path, method):
    """ Calculates the partial charges or each ligand using ELF conformer
    selection and then averaging over partial charges generated with
    am1bcc. The poses the ligand using ligpose.py OEPosit() method.
    Outputs corrected ligand locations in *.sdf and *.mol2.
    """
    for cell in SHEET['B']:
        # Shift to cell A in *.xlsx to read ligand name.
        cell_shift = cell.row - 1
        global LIG_NAME
        LIG_NAME = SHEET['A'][cell_shift].value
        smiles = str(cell.value)
        #generate the charged molecule from expanded conformations, am1bcc
        charged_molecule = elf_pc_gen(smiles)
        #and align the charged_molecule
        print('Preparing ligand {} for flexible positioning'
              ' using OEPosit().'
              .format(cell.value))
        ligpose.main(pdb_path, charged_molecule, method)
        #ligpose.main(pdb_path, mtz_path, charged_molecule, method)
        print('Ligand {} aligned onto receptor bound ligand'
              .format(cell.value))
              

def rigid_overlay():
    ''' This function will do a rigid rotation,
    calculate the common core and the maximum common
    substructure of the ligand to be positioned and
    the resolved ligand used as reference.
    '''
    working_dir = os.getcwd()
    os.makedirs(f'tmp_outputs', exist_ok=True)
    
    # Read in the bound complex reference structure
    complex_in = oemolistream("complex.pdb")
    mol = OEGraphMol()
    OEReadMolecule(complex_in, mol)

    lig = OEGraphMol()
    prot = OEGraphMol()
    wat = OEGraphMol()
    other = OEGraphMol()

    # Extract from the complex pdb file the subparts
    OESplitMolComplex(lig, prot, wat, other, mol)
    ofs = oemolostream("ligand.pdb")
    OEWriteConstMolecule(ofs, lig)
    ofs.close()
    
    for cell in SHEET['A']:
        #currently using only one file for testing
        lig_target = oemolistream("{}.mol2".format(cell.value))
        os.rename(f'{cell.value}.mol2',
                  f'tmp_outputs/{cell.value}.mol2')
        # Should put in a flag to do this if it exists
        os.rename(f'{cell.value}.sdf',
                  f'tmp_outputs/{cell.value}.sdf')
        lig_move = OEGraphMol()
        OEReadMolecule(lig_target, lig_move)
   
        # I might need to generate conformations?
        sub.call("python oescripts/mcs3dalign.py ligand.pdb {}.mol2\
                 {}_aligned.mol2".format(cell.value), shell=True)

    # Find maximum common substructure smi in dataset of ligands
    # Doesn't work for disconnected common substructures.
    file = open("smiles.smi", "w")
    #generate the lig smile and write to file
    for cell in SHEET['B']:
        file.write(str(cell.value) + "\n")
    file.close()
    sub.call("python3 oescripts/getcore.py -in smiles.smi", shell=True)
    print("The common substructure found was saved as \
          common_substructure.smi")
    print("Initial *.sdf and *mol2 files were moved to tmp_outputs. ")
    
 
def struc_dir_setup(system_dir, date_dir, pmx_ligands_dir):
    """ Organizes the folders for *.sdf and *.mol2. Consistent with
    pmxworkflow. Input variables:
    system_dir is the directory for the simulated system.
        example: tyk2
    date_dir = date of the project,
        default 02_ligands
    pmx_ligands_dir = the folder that will store ligand info
    Relevant to pmxworkflow:
        system_dir/date_dir/pmx_ligands_dir = pwf.ligPath
    """

    working_dir = os.getcwd()
    # Make folder structure for pmxworkflow.
    os.makedirs(f'{system_dir}/{date_dir}/{pmx_ligands_dir}',
                exist_ok=True)
    # If prior folders were created here, back them up.
    if len(os.listdir(f'{system_dir}/{date_dir}/{pmx_ligands_dir}')) != 0:
        os.chdir(f'{system_dir}/{date_dir}/{pmx_ligands_dir}')
        # Backup directories
        time_ran = datetime.now()
        current_time = time_ran.strftime("%H_%M_%S")

        for existing_dir in os.listdir():
            bu_name = "backup_{}_{}".format(current_time, existing_dir)
            os.rename(existing_dir, bu_name)
        print("Your current working directory where files from this"
              " run were ouput is:")
        os.chdir(f'{working_dir}')
        print(os.getcwd())
        print("\nPrior folders in {}/{}/{}\nwere backed up to backup_{}_*"
              .format(system_dir, date_dir, pmx_ligands_dir, current_time))
    
    # Tests.
    made_sdf_files = [f for f in os.listdir(f'{working_dir}')
                      if f.endswith('.sdf')]
    if len(made_sdf_files) == 0:
        print("Warning: no *.sdf files found for ligands in current working \
              directory, {}.".format(working_dir))
    made_mol2_files = [f for f in os.listdir(f'{working_dir}')
                       if f.endswith('.mol2')]
    if len(made_mol2_files) == 0:
        print("Warning: no *.mol2 files found for ligands in current working \
              directory, {}.".format(working_dir))
        
    # Make and move structural file directories.
    for cell in SHEET['A']:
        os.makedirs(f'{system_dir}/{date_dir}/{pmx_ligands_dir}/{cell.value}/crd',
                    exist_ok=True)
        # If the .sdf and .mol2 files exist, move them to created directories
        if os.path.exists(f'{working_dir}/{cell.value}.sdf') == True:
            os.rename(f'{working_dir}/{cell.value}.sdf',
                      f'{system_dir}/{date_dir}/{pmx_ligands_dir}/{cell.value}/crd/{cell.value}.sdf')
            print("File {}.sdf deposited at {}/{}/{}/{}/crd"
                  .format(cell.value, system_dir, date_dir, pmx_ligands_dir,
                  cell.value)
                  )
        if os.path.exists(f'{working_dir}/{cell.value}.mol2') == True:
            os.rename(f'{working_dir}/{cell.value}.mol2',
                      f'{system_dir}/{date_dir}/{pmx_ligands_dir}/{cell.value}/crd/{cell.value}.mol2')
            print("File {}.mol2 deposited at {}/{}/{}/{}/crd"
                  .format(cell.value, system_dir, date_dir, pmx_ligands_dir,
                  cell.value)
                  )
    print("\n**   When you run the pmxworkflow set the pwf path input argument"
          " to {}  **\n".format(system_dir))
    print("Ligand generation and system building completed at {}.\n"
          "------------------------------------------------------\n\n"
          .format(current_time))
    
