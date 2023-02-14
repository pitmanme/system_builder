#!/usr/bin/env python
# coding: utf-8
#
# ----------------------------------------------------------------------
# Run this script to test only ligand charge generations.
# ----------------------------------------------------------------------
__doc__="""
        Run this script to test only ligand charge generations.
        
        arguments for main(): (pdb_path, mtz_path, lig_in, method)
        Input arguments are str. Example: 'shapefit'
        
        input_excel_file =  Location of the input excel file with ligand names [col A] and ligand smile + name in [col B].
        
        pdb_path: path to the *.pdb file containing your reference
                  complex structure.
                
        
        Outputs:
            - *.sdf and *.mol2 files for building simulations.
            - *.svg files to visualize posed ligand interactions and
              compare to the reference complex interactions.
         
        """

__author__ = 'M. Pitman'

import sys

import ligbuild
from ligpose import Logger

# ----------------------------------------------------------------------
# Define input variables
# ----------------------------------------------------------------------

# Location of the input excel file with ligand names [col A] and
# ligand smile + name in [col B].
input_excel_file = 'smi_codes.xlsx'
toolkit_to_use = 'openeye'

# Reference complex structures for ligand positioning.
pdb_path = '7cmd.pdb'

# ----------------------------------------------------------------------

# Set up logfile output
sys.stdout = Logger()
print("Log file generated at ligpose.log.")

# Read in excel file.
ligbuild.collect_excel_data(input_excel_file)

# Generate pre_*.sdf loops through mols in excel input.
ligbuild.smi_to_sdf(toolkit_to_use, input_excel_file)

#Calculates charges
ligbuild.sdf_to_mol2_multiconf()

