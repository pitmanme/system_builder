#!/usr/bin/env python
# coding: utf-8
#
# ----------------------------------------------------------------------
# system_builder.py, version 1.0, Mary Pitman, Mobley Lab UCI
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appears in supporting documentation.
# Later this should be updated for argument parsing and command line
# input option.
# ----------------------------------------------------------------------
__doc__="""
        Modules designed to posit ligands using the OEPosit space of
        Openeye. This module performs ligand posing through
        information from existing complex structures or through
        docking.
        
        arguments for main(): (pdb_path, mtz_path, lig_in, method)
        Input arguments are str. Example: 'shapefit'
        
        input_excel_file =  Location of the input excel file with ligand names [col A] and ligand smile + name in [col B].
        
        pdb_path: path to the *.pdb file containing your reference
                  complex structure.
        mtz_path: *.mtz file for the reference complex structure.
                  *.mtz files can be downloaded from the PDB.
        method:   the method you wish to use for OEPosit() to
                  perform the flexible ligand fitting.
                  
                  Methods:
                  shapefit, mcs, hybrid, fred.
                
        
        Outputs:
            - *.sdf and *.mol2 files for building simulations.
            - *.svg files to visualize posed ligand interactions and
              compare to the reference complex interactions.
            - a receptor.pdb to generate the updated bound complex.
        """
        
__version__ = '1.0'
__author__ = 'M. Pitman'

import sys

import ligbuild
from ligpose import Logger

# ----------------------------------------------------------------------
# Define input variables
# ----------------------------------------------------------------------

# Location of the input excel file with ligand names [col A] and
# ligand smile + name in [col B].
input_excel_file = '/Users/mpitman/work/system_tests/nos/6xk8/bNOS_inhibs.xlsx'
toolkit_to_use = 'openeye'

# Reference complex structures for ligand positioning.
pdb_path = '/Users/mpitman/work/system_tests/nos/6xk8/6xk8.pdb'
mtz_path = '/Users/mpitman/work/system_tests/Nos/6xk8/6xk8_phases.mtz'
method = 'shapefit'

# The directory input variables
system_dir = '/Users/mpitman/work/rfe/tests/6xk8'
date_dir = '2021-11-03_6xk8'
pmx_ligands_dir = '02_ligands'

# ----------------------------------------------------------------------

# Set up logfile output
sys.stdout = Logger()
print("Log file generated at ligpose.log.")

# Read in excel file.
ligbuild.collect_excel_data(input_excel_file)

# Generate pre_*.sdf loops through mols in excel input.
ligbuild.smi_to_sdf(toolkit_to_use, input_excel_file)

# Perform flexibe ligand positioning and calculate ligand charges.
#ligbuild.ligand_reshaping(pdb_path, mtz_path, method)
ligbuild.ligand_reshaping(pdb_path, method)

# Generate the needed folders and move files to correct locations.
# Existing redundant files will be backed up as reported on terminal
# output in logfile.log.
ligbuild.struc_dir_setup(system_dir, date_dir, pmx_ligands_dir)

