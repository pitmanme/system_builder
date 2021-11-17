import sys
from openeye import oechem

input = '/Users/mpitman/work/system_tests/Nos/6xk8/pmxworkflow/2021-10-19/02_ligands/V54/crd/V54.sdf'
output = '/Users/mpitman/work/system_tests/Nos/6xk8/pmxworkflow/2021-10-19/02_ligands/V54/crd/isocan_V54.smi'

def gen_isocan(input, output):
    ''' short script to ouput the protinated isocan *.mol2'''
    ifs = oechem.oemolistream(input)
    ofs = oechem.oemolostream(output)

    ifs.SetFormat(oechem.OEFormat_SDF)
    ofs.SetFormat(oechem.OEFormat_ISM)

    for mol in ifs.GetOEGraphMols():
        oechem.OESetSDData(mol, "$SMI", oechem.OECreateIsoSmiString(mol))
        oechem.OEWriteMolecule(ofs, mol)
    ifs.close()
    ofs.close()

gen_isocan(input, output)
