#!/usr/bin/env python
# Excerpt from Openeye
# (C) 2017 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms.
#
# This is a modified version of the sample code for module use

#############################################################################
# Print protein-ligand interactions
#############################################################################
import sys

from openeye import oechem

#pdb_path = '/Users/mpitman/work/system_tests/7cmd/structure_files/uncapped.pdb'

def du2liginters(du):

    ligand = oechem.OEGraphMol()
    du.GetLigand(ligand)
    
    protein = oechem.OEGraphMol()
    du.GetProtein(protein)


    if ligand.NumAtoms() == 0:
        oechem.OEThrow.Fatal("Cannot separate complex!")

    # Perceive interactions

    asite = oechem.OEInteractionHintContainer(protein, ligand)
    if not oechem.OEIsValidActiveSite(asite):
        oechem.OEThrow.Fatal("Cannot initialize active site!")

    oechem.OEPerceiveInteractionHints(asite)

    print("Number of interactions:", asite.NumInteractions())
    for itype in oechem.OEGetActiveSiteInteractionHintTypes():
        numinters = asite.NumInteractions(oechem.OEHasInteractionHintType(itype))
        if numinters == 0:
            continue
        print("%d %s :" % (numinters, itype.GetName()))

        inters = [s for s in asite.GetInteractions(oechem.OEHasInteractionHintType(itype))]
        print("\n".join(sorted(GetInteractionString(s) for s in inters)))

    print("\nResidue interactions:")
    for res in oechem.OEGetResidues(asite.GetMolecule(oechem.OEProteinInteractionHintComponent())):
        PrintResidueInteractions(asite, res)

    print("\nLigand atom interactions:")
    for atom in asite.GetMolecule(oechem.OELigandInteractionHintComponent()).GetAtoms():
        PrintLigandAtomInteractions(asite, atom)


def GetResidueName(residue):
    return "%3s %4d %s" % (residue.GetName(), residue.GetResidueNumber(), residue.GetChainID())


def GetInteractionString(inter):

    fragstrs = []
    for frag in [inter.GetBgnFragment(), inter.GetEndFragment()]:
        if frag is None:
            continue
        fragtype = frag.GetComponentType()
        if fragtype == oechem.OELigandInteractionHintComponent():
            fragstrs.append("ligand:" + " ".join(sorted(str(a) for a in frag.GetAtoms())))
        if fragtype == oechem.OEProteinInteractionHintComponent():
            fragstrs.append("protein:" + " ".join(sorted(str(a) for a in frag.GetAtoms())))

    return " ".join(sorted(f for f in fragstrs))


def PrintResidueInteractions(asite, residue):

    ligatomnames = set()
    for inter in asite.GetInteractions(oechem.OEHasResidueInteractionHint(residue)):
        ligfrag = inter.GetFragment(oechem.OELigandInteractionHintComponent())
        if ligfrag is None:
            continue
        for latom in ligfrag.GetAtoms():
            ligatomnames.add(str(latom))

    if len(ligatomnames) != 0:
        print(GetResidueName(residue), ": ", " ".join(sorted(a for a in ligatomnames)))


def PrintLigandAtomInteractions(asite, atom):

    resnames = set()
    for inter in asite.GetInteractions(oechem.OEHasInteractionHint(atom)):
        profrag = inter.GetFragment(oechem.OEProteinInteractionHintComponent())
        if profrag is None:
            continue
        for patom in profrag.GetAtoms():
            residue = oechem.OEAtomGetResidue(patom)
            resnames.add(GetResidueName(residue))

    if len(resnames) != 0:
        print(atom, ":", " ".join(sorted(r for r in resnames)))


#print_ligand_interactions(pdb_path)
