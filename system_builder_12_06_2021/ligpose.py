#!/usr/bin/env python
# coding: utf-8
#
# The main change that is needed here is that it needs to:
#    read in the charged molecule from ligand_builder
# ----------------------------------------------------------------------
# ligpose.py, version 1.0, Mary Pitman, Mobley Lab UCI
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
        Modules designed to posit ligands using the OEPosit space of
        Openeye. This module performs ligand posing through
        information from existing complex structures or through
        docking.
        
        arguments for main(): (pdb_path, mtz_path, lig_in, method)
        Input arguments are str. Example: 'shapefit'
        
        pdb_path: path to the *.pdb file containing your reference
                  complex structure.
        mtz_path: *.mtz file for the reference complex structure.
                  *.mtz files can be downloaded from the PDB.
        lig_in:   The *.mol2 file generated from ligbuild.py or
                  through other means. Can take *.pdb and *.sdf.
                  Any 3D structure.
        method:   the method you wish to use for OEPosit().
                  shapefit, mcs, hybrid, fred."
                
        
        Outputs:
            - *.sdf and *.mol2 files for building simulations.
            - *.svg files to visualize posed ligand interactions and
              compare to the reference complex interactions.
            - a receptor.pdb to generate the updated bound complex.
        """
__version__ = '1.0'
__author__ = 'M. Pitman'

import os
import sys
import re
import subprocess as sub

from openeye import oechem
from openeye import oespruce
from openeye import oedocking
from openeye import oegrid

from oescripts import du2liginters
# ----------------------------------------------------------------------

def extract_ligand(du, ofs):
    ''' Extracts and checks for a ligand in a complex.'''
    ligand = oechem.OEGraphMol()
    if not du.HasLigand():
        oechem.OEThrow.Fatal("Error: There is no ligand in the"
                             " OEDesignUnit."
                             )
    oechem.OEWriteMolecule(ofs, ligand)
    ofs.close()
    
    return ligand
    

def build_du(pdb_path):
#def build_du(pdb_path, mtz_path):
    ''' Builds the initial design unit. '''
    # Read the input pdb file
    iname = pdb_path
    ifs = oechem.oemolistream()
    
    if not ifs.open(iname):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % iname)

    if not oechem.OEIs3DFormat(ifs.GetFormat()):
        oechem.OEThrow.Fatal("%s is not in a 3D format." % iname)
    sopt = oechem.OESplitMolComplexOptions()
    '''
    edfile = mtz_path
    ed = oegrid.OESkewGrid()
    if not oegrid.OEReadMTZ(edfile, ed, oegrid.OEMTZMapType_Fwt):
        oechem.OEThrow.Fatal(
            "Unable to read electron density file %s" % edfile
        )
    include_ed = True
    '''
    ifs.SetFlavor(
        oechem.OEFormat_PDB,
        oechem.OEIFlavor_PDB_Default
        | oechem.OEIFlavor_PDB_DATA
        | oechem.OEIFlavor_PDB_ALTLOC,
    )
    
    complexmol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(ifs, complexmol):
        oechem.OEThrow.Fatal("Unable to read %s." % iname)

    ifs.close()
    
    metadata = oespruce.OEStructureMetadata()
    opts = oespruce.OEMakeDesignUnitOptions()
    
    #dus = oespruce.OEMakeDesignUnits(complexmol, ed, metadata, opts)
    dus = oespruce.OEMakeDesignUnits(complexmol, metadata, opts)
    du_ofile = os.path.basename(iname)[:-4] + "_DU.oedu"
    print("\n--------------------------------------------------\n"
          "Saved a design unit binary file to {}"
          .format(du_ofile))
    for i, du in enumerate(dus):
        oechem.OEWriteDesignUnit(du_ofile, du)
    oespruce.OEProtonateDesignUnit(du)
    
    return du, du_ofile
 

def make_receptor(oedu_path):
    '''
    Defines how to make the receptor that will be
    added to the design unit.
    '''
    recOpts = oedocking.OEMakeReceptorOptions()
    iname = oedu_path
    ifs = oechem.oeifstream()
    ifs.open(iname)

    rec_ofile = 'receptor_out.oedu'
    ofs = oechem.oeofstream(rec_ofile)

    du = oechem.OEDesignUnit()
    while oechem.OEReadDesignUnit(ifs, du):
        if oedocking.OEMakeReceptor(du, recOpts):
            oechem.OEWriteDesignUnit(ofs, du)
        else:
            oechem.OEThrow.Warning("Failed to make receptor.")
    
    ifs.close()
    ofs.close()
    
    return rec_ofile
    

def flexible_overlay(pdb_path, lig_in, method):
#def flexible_overlay(pdb_path, mtz_path, lig_in, method):
    '''
    Positions the ligand, sets options for OEPosit(),
    and outputs ligand and receptor structure files.
    '''
    #du, du_ofile = build_du(pdb_path, mtz_path)
    du, du_ofile = build_du(pdb_path)
    rec_ofile = make_receptor(du_ofile)
    
    # Read output from make_receptor()
    ifs = oechem.oeifstream()
    if not ifs.open(rec_ofile):
        oechem.OEThrow.Fatal("Unable to open %s for reading"
                             % rec_ofile)
                             
    du = oechem.OEDesignUnit()
    oechem.OEReadDesignUnit(ifs, du)
    if du.HasReceptor():
        print("Receptor built for the design unit, du.")
    if not du.HasReceptor():
        print("The design unit does not have a receptor built")
        
    # Set posing options.
    if method == 'shapefit':
        pose_method = oedocking.OEPositMethod_SHAPEFIT
    elif method == 'mcs':
        pose_method = oedocking.OEPositMethod_MCS
    elif method == 'hybrid':
        pose_method = oedocking.OEPositMethod_HYBRID
    elif method == 'fred':
        pose_method = oedocking.OEPositMethod_FRED
    else:
        print("Method argument required for flexible overlay."
              " args: pdb_path, mtz_path, lig_file, method"
              " methods(str): shapefit, mcs, hybrid, fred."
              )
    positOpts = oedocking.OEPositOptions()
    positOpts.SetFullConformationSearch(True)
    positOpts.SetPositMethods(pose_method)
    positOpts.SetIgnoreNitrogenStereo(True)
    
    poser = oedocking.OEPosit(positOpts)
    poser.AddReceptor(du)
    
    print("Posing will be generated with the method {} "
          "over a full conformational search, with ambiguous "
          "nitrogen stereocenters ignored.\n "
          .format(method))
    
    ifs.close()
          
    # Read in the ligand to be positioned.
    try:
        lig_to_pose = oechem.OEMol(lig_in)
    except TypeError:
        ims = oechem.oemolistream(lig_in)
        lig_to_pose = oechem.OEMol()
        oechem.OEReadMolecule(ims, lig_to_pose)
    
        if not ims.open(lig_in):
            oechem.OEThrow.Fatal("Unable to open %s for reading"
                                 % lig_in)
        ims.close()
        
    if lig_to_pose.GetTitle() == '':
        lig_to_pose.SetTitle('posed_lig')
        print("The ligand to be posed did not have a title. "
              "Ligand title set to posed_lig.")

    # Do the posing.
    print("Posing in progress ...")
    result = oedocking.OESinglePoseResult()
    poser_output = poser.Dock(result, lig_to_pose)

    if poser_output == oedocking.OEDockingReturnCode_Success:
        # Parse the DU.
        posed_du = result.GetDesignUnit()
        posed_oemol = oechem.OEGraphMol()
        posed_du.GetLigand(posed_oemol)
        du_protein = oechem.OEGraphMol()
        posed_du.GetProtein(du_protein)
        posed_du.SetDoubleData(poser.GetName(), result.GetProbability())
        print("\nThe ligand %s was successfully positioned."
              % lig_to_pose.GetTitle())
        
        # Save *.sdf and *.mol2 files of the ligand.
        extensions = ('sdf', 'mol2')
        for i in extensions:
            ofile = 'shapefit_{}.{}'.format(lig_to_pose.GetTitle(), i)
            ostream = oechem.oemolostream()
            ostream.open(ofile)
            oechem.OEWriteMolecule(ostream, posed_oemol)
            ostream.close()
        lig_out = 'shapefit_{}.sdf'.format(lig_to_pose.GetTitle())
        
        # Save *.pdb of the receptor.
        ofile_protein = 'receptor_protein.pdb'
        ostream = oechem.oemolostream()
        ostream.open(ofile_protein)
        oechem.OEWriteMolecule(ostream, du_protein)
        ostream.close()
        
        # Prints the posed ligand interactions with the protein.
        print("--------------------------------------------------\n"
              "List of interactions between the posed ligand\nand "
              "the protein:"
              "\n--------------------------------------------------")
        du2liginters(posed_du)
        
    else:
        error_mes = oedocking.OEDockingReturnCodeGetName(poser_output)
        oechem.OEThrow.Warning("%s: %s" % (lig_to_pose.GetTitle(),
                               error_mes))

    return lig_out
  
  
def visualize_contacts(pdb_path, lig_out):
    '''
    Uses an openeye script contained in oescripts to create
    a visualization of how the ligand interacts with the structure
    in the reference complex and the generated posed complex.
    
    example of how to open:
        open -a "Google Chrome" ouput.svg
    '''
    # Create the visualization of the reference complex
    # interactions.
    resolved_inter_ofile = os.path.basename(pdb_path)[:-4] \
                           + "_interactions.svg"
    sub.call("python3 oescripts/complex2img.py -complex {} -out {}"
             .format(pdb_path, resolved_inter_ofile), shell = True)
    print("\n--------------------------------------------------\n"
          "The interactions between the reference complex "
          "structure and resolved ligand were visualized at:\n{}"
          .format(resolved_inter_ofile))
          
    # Create the visualization of the posed ligand interactions.
    ofile2 = 'receptor_protein.pdb'
    #ofile2 = 'receptor_out.oedu'
    base_ofile = lig_out.split(".")[0]
    inter_ofile = '{}_interactions.svg'.format(base_ofile)
    sub.call("python3 oescripts/complex2img.py -protein {} "
             "-ligand {} -out {}"
             .format(ofile2, lig_out, inter_ofile), shell = True)
    print("The interactions between the posed ligand and "
          "the receptor were visualized at:\n{}\n"
          "--------------------------------------------------"
          .format(inter_ofile))


def clean_up(files):
    '''
    Removes *.oedu files generated by ligpose.py.
    Input: list of files to remove.
    '''
    print("Cleaning up files from ligand posing...")
    for i in files:
        if os.path.exists(i):
            os.remove(i)
            print("Design unit binary file {} removed by"
                  " ligpose.clean_up".format(i))
        else:
            print("The file {} does not exist for clean up."
                  .format(i))


class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("ligpose.log", "a")

    def write(self, outputs):
        self.terminal.write(outputs)
        self.log.write(outputs)

    def flush(self):
        pass
 
 
#def main(pdb_path, mtz_path, lig_in, method):
def main(pdb_path, lig_in, method):
    #lig_out = flexible_overlay(pdb_path, mtz_path, lig_in, method)
    lig_out = flexible_overlay(pdb_path, lig_in, method)
    visualize_contacts(pdb_path, lig_out)
    # Remove the generated binary files.
    to_clean = [os.path.basename(pdb_path)[:-4] + "_DU.oedu",
                'receptor_out.oedu'
                ]
    clean_up(to_clean)


if __name__ == '__main__':
    main()
