# clean up so you don't have to put the xlsx in
# -------------------------------------------------------------------------------
import os
import sys

from openeye.oechem import *
from openeye.oedocking import *
from openeye.oeomega import *
from openeye import oedocking
from openeye import oechem
from openeye import oeomega

import openpyxl
from oescripts import *

#---------------------------------------------------------------------------
def collect_excel_data(input_excel_file: str) -> str:
    """ Reads in 'my_data.xlsx' """
    try:
        global MY_DATA, SHEET
        MY_DATA = openpyxl.load_workbook(input_excel_file)
        SHEET = MY_DATA['Sheet1']
    except:
        print("The argument for read_excel_file is 'my_data.xlsx' .")
        sys.exit()

# Takes in one complex to form a receptor and then fits multiple ligands
def shape_fitter(file):
    # For testing.
    collect_excel_data('/Users/mpitman/work/rfe/tests/tyk2_2/tyk2_2_lig_data.xlsx')
    
    # Input is a string
    working_dir = '.'
    
    # Load structual pdb file
    #ifile = oemolistream(os.path.join(working_dir, file))
    ifile = oemolistream(file)
    complex_mol = OEMol()
    OEReadMolecule(ifile, complex_mol)
    ifile.close()

    # Make receptors of the hosts for use in docking; takes about 4 minutes on my Mac
def create_receptor(pdb_path):
    """Create an OpenEye receptor from a PDB file.
    Parameters
    ----------
    protein_pdb_path : str
        Path to the receptor PDB file.
    box : 1x6 array of float
        The minimum and maximum values of the coordinates of the box
        representing the binding site [xmin, ymin, zmin, xmax, ymax, zmax].
    Returns
    -------
    receptor : openeye.oedocking.OEReceptor
        The OpenEye receptor object.
    """
    input_mol_stream = oemolistream(pdb_path)
    complex_mol = OEGraphMol()
    OEReadMolecule(input_mol_stream, complex_mol)

    com = OEFloatArray(3)
    OEGetCenterOfMass(complex_mol, com)
    receptor = OEGraphMol()
    oedocking.OEMakeReceptor(receptor, complex_mol, com[0], com[1], com[2])
    # This will need to be edited for when a whole path is given
    receptor_name = pdb_path.split(".")[0]
    # I might need to be outputting an o
    #ofs = oemolistream('test.oeb')
    oedocking.OEWriteReceptorFile(receptor,'test.oeb')

    return receptor
    
def docking_receptor(pdb_path, coords):
    """Create an OpenEye receptor from a PDB file.
    Parameters
    ----------
    protein_pdb_path : str
        Path to the receptor PDB file.
    box : 1x6 array of float
        The minimum and maximum values of the coordinates of the box
        representing the binding site [xmin, ymin, zmin, xmax, ymax, zmax].
    Returns
    -------
    receptor : openeye.oedocking.OEReceptor
        The OpenEye receptor object.
    """
    input_mol_stream = oemolistream(pdb_path)
    complex_mol = OEGraphMol()
    OEReadMolecule(input_mol_stream, complex_mol)
    
    scaler = 10
    box = [coords[0] - scaler, coords[1] - scaler, coords[2]- scaler, coords[0] + scaler, coords[1] + scaler, coords[2] + scaler]
    box = oedocking.OEBox(*box)
    receptor = OEGraphMol()
    #oedocking.OEMakeReceptor(receptor, complex_mol, coords[0], coords[1], coords[2])
    oedocking.OEMakeReceptor(receptor, complex_mol, box)
    # This will need to be edited for when a whole path is given
    receptor_name = pdb_path.split(".")[0]
    print('Receptor generated from {} for docking of ligands'.format(pdb_path))
    
    return receptor
    

    
def create_bound_receptor(pdb_path):
    """Create an OpenEye receptor from a PDB file and ligand mol2 file.
    Parameters
    ----------
    protein_pdb_path : str
        Path to the receptor PDB file.
    ligand_file_path : str
        path to the ligand file. Can be any file supported by openeye
    Returns
    -------
    receptor : openeye.oedocking.OEReceptor
        The OpenEye receptor object.
    """
    # Load in the pdb
    input_mol_stream = oemolistream(pdb_path)
    complex_mol = OEGraphMol()
    OEReadMolecule(input_mol_stream, complex_mol)
    
    # Find ligand
    lig = oechem.OEGraphMol()
    prot = oechem.OEGraphMol()
    wat = oechem.OEGraphMol()
    other = oechem.OEGraphMol()
    # Extract from the complex pdb file the subparts
    oechem.OESplitMolComplex(lig, prot, wat, other, complex_mol)

    receptor = oechem.OEGraphMol()
    oedocking.OEMakeReceptor(receptor, complex_mol, lig)
    print('Receptor generated from {} for alignment of ligands'.format(pdb_path))
    
    return receptor
        
  
def generate_conformers(molecule, max_confs=800, strictStereo=True, ewindow=15.0, rms_threshold=1.0, strictTypes = True):
    """Generate conformations for the supplied molecule
  
    -----
    This is adapted from openmoltools! To avoid version incompatabiltiies.
    Question: is this still required?
    """

    molcopy = oechem.OEMol(molecule)
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(max_confs)
    omega.SetIncludeInput(True)
    omega.SetCanonOrder(False)

    omega.SetSampleHydrogens(True)
    omega.SetEnergyWindow(ewindow)
    omega.SetRMSThreshold(rms_threshold)

    omega.SetStrictStereo(strictStereo)
    omega.SetStrictAtomTypes(strictTypes)

    omega.SetIncludeInput(False)  # don't include input
    if max_confs is not None:
        omega.SetMaxConfs(max_confs)

    status = omega(molcopy)  # generate conformation
    if not status:
        raise(RuntimeError("omega returned error code %d" % status))

    return molcopy

    
# Maybe have this read in the charged molecules returned from other func
# Pass in charged_molecule from ligand_builder_v0_5 assign_partial_charges() func
def pose_molecule_from_charged(receptor, cell, charged_molecule, n_conformations=25, n_poses=1):
    
    """Run the multi-conformer docker.
    Parameters
    ----------
    receptor : openeye.oedocking.OEReceptor
        The openeye receptor.
    molecule_smiles : str
        The SMILES string of the molecule.
    n_conformations : int, optional
        The number of omega conformations to pass to the multi-conformer
        docker (default is 10).
    n_poses : int, optional
        Number of binding poses to return.
    Returns
    -------
    docked_oemol : openeye.oechem.OEMol
        The docked multi-conformer OpenEye molecule.
    """
    #generate multi conformer mol2 or molecule with charges
    # SEE IF THERE IS SOMETHING IN OPENEYE OR OFF TO DO THIS INSTEAD
    molecule_oemol = generate_conformers(charged_molecule, max_confs=n_conformations)

    # Write out the predocked version for testing
    outfile = 'predocking_{}.mol2'.format(cell.value)
    ostream = oemolostream()
    ostream.open(outfile)
    OEWriteMolecule(ostream, molecule_oemol)
    ostream.close()

    poser = oedocking.OEPosit()
    poser.Initialize(receptor)
    
    # Initiate and do the initial posing
    posed_oemol = oechem.OEMol()
    poser.DockMultiConformerMolecule(posed_oemol, molecule_oemol, n_poses)
    
    # Write out the docked mol2 file
    ofile = '{}.mol2'.format(cell.value)
    ostream = oemolostream()
    ostream.open(ofile)
    OEWriteMolecule(ostream, posed_oemol)
    ostream.close()
    
    # Write out the docked sdf file
    ofile = '{}.pdb'.format(cell.value)
    ostream = oemolostream()
    ostream.open(ofile)
    OEWriteMolecule(ostream, posed_oemol)
    ostream.close()
    
    # Write out the receptor after docking
    orfile = 'receptor_postdocking_of_{}.pdb'.format(cell.value)
    ostream = oemolostream()
    ostream.open(orfile)
    OEWriteMolecule(ostream, receptor)
    ostream.close()

    return posed_oemol

def test2_recipe():
    receptor = create_bound_receptor('struc.pdb')
    molecule_smiles = "[H]c1c(c(c(c(c1[H])Cl)C(=O)N([H])c2c(c(nc(c2[H])N([H])C(=O)C([H])([H])[H])[H])[H])Cl)[H]"
    pose_molecule_from_charged(receptor, cell, charged_molecule, n_conformations=25, n_poses=1)
    

