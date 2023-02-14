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
        Modules designed to dock ligands using the OEPosit space of
        Openeye. This module performs ligand posing through
        information from existing complex structures or through
        docking.
        
        input_excel_file =  Location of the input excel file with ligand names [col A] and ligand smile + name in [col B].
        
        pdb_path: path to the *.pdb file containing your reference
                  complex structure.
        
        coords: array of guess coordinates
                
        
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
input_excel_file = '/Users/mpitman/work/aimee_figs/683_dock.xlsx'
toolkit_to_use = 'openeye'

# Reference complex structures for ligand positioning.
pdb_path = '/Users/mpitman/work/aimee_figs/importinb.pdb'

coords = [29.590, 63.695, 45.154]
coords1 = [37.449, 11.236, 22.176]
coords2 = [56.307, -17.271, 9.723]
coords3 = [63.406, -15.448, 4.408]
coords4 = [29.554, 2.195, 30.429]

# ----------------------------------------------------------------------

# Set up logfile output
#sys.stdout = Logger()
#print("Log file generated at ligpose.log.")

# Read in excel file.
ligbuild.collect_excel_data(input_excel_file)

# Generate pre_*.sdf loops through mols in excel input.
ligbuild.smi_to_sdf(toolkit_to_use, input_excel_file)


# Docks ligand based on guess
ligbuild.ligand_docking(pdb_path, coords)



