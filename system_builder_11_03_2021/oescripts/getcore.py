#!/usr/bin/env python3
# (C) 2017 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.
#
#
# This version contains edits by M Pitman. 
#############################################################################
# Returns the "best" i.e. largest common fragment of a set of molecules
#############################################################################

import sys
from openeye import oechem
from openeye import oegraphsim


def main(argv=[__name__]):

    itf = oechem.OEInterface(InterfaceData)

    if not oechem.OEParseCommandLine(itf, argv):
        return 1

    iname = itf.GetString("-in")

    # check input file

    ifs = oechem.oemolistream()
    if not ifs.open(iname):
        oechem.OEThrow.Fatal("Cannot open input file!")

    # read molecules

    mollist = []
    for idx, mol in enumerate(ifs.GetOEGraphMols()):
        oechem.OESuppressHydrogens(mol)
        mollist.append(oechem.OEGraphMol(mol))

    minbonds = itf.GetInt("-minbonds")
    maxbonds = itf.GetInt("-maxbonds")

    bestcore = GetCoreFragment(mollist, minbonds, maxbonds)
    # Here are my edits
    print("Core fragment = %s" % oechem.OEMolToSmiles(bestcore))
    f = open("common_substructure.smi","w")
    f.write(oechem.OEMolToSmiles(bestcore))
    f.close()
 #   print('%s' % oechem.OEMolToSmiles(bestcore), file=f)
    return 0


def GetReferenceMolecule(mollist):

    refmol = None
    for mol in mollist:
        if refmol is None or mol.NumAtoms() < refmol.NumAtoms():
            refmol = mol

    return refmol


def GetFragments(mol, minbonds, maxbonds):

    frags = []
    fptype = oegraphsim.OEGetFPType("Tree,ver=2.0.0,size=4096,bonds=%d-%d,atype=AtmNum,btype=Order"
                                    % (minbonds, maxbonds))

    for abset in oegraphsim.OEGetFPCoverage(mol, fptype, True):
        fragatompred = oechem.OEIsAtomMember(abset.GetAtoms())

        frag = oechem.OEGraphMol()
        adjustHCount = True
        oechem.OESubsetMol(frag, mol, fragatompred, adjustHCount)
        oechem.OEFindRingAtomsAndBonds(frag)
        frags.append(oechem.OEGraphMol(frag))

    return frags


def GetCommonFragments(mollist, frags,
                       atomexpr=oechem.OEExprOpts_DefaultAtoms,
                       bondexpr=oechem.OEExprOpts_DefaultBonds):

    corefrags = []

    for frag in frags:

        ss = oechem.OESubSearch(frag, atomexpr, bondexpr)
        if not ss.IsValid():
            continue

        validcore = True
        for mol in mollist:
            validcore = ss.SingleMatch(mol)
            if not validcore:
                break

        if validcore:
            corefrags.append(frag)

    return corefrags


def GetCoreFragment(mols, minbonds, maxbonds,
                    atomexpr=oechem.OEExprOpts_DefaultAtoms,
                    bondexpr=oechem.OEExprOpts_DefaultBonds):

    print("Number of molecules = %d" % len(mols))

    refmol = GetReferenceMolecule(mols)

    frags = GetFragments(refmol, minbonds, maxbonds)
    if len(frags) == 0:
        oechem.OEThrow.Error("No fragment is enumertated with bonds %d-%d!" % (minbonds, maxbonds))

    print("Number of fragments = %d" % len(frags))

    commonfrags = GetCommonFragments(mols, frags, atomexpr, bondexpr)
    if len(commonfrags) == 0:
        oechem.OEThrow.Error("No common fragment is found!")

    print("Number of common fragments = %d" % len(commonfrags))

    core = None
    for frag in commonfrags:
        if core is None or GetFragmentScore(core) < GetFragmentScore(frag):
            core = frag

    return core


def GetFragmentScore(mol):

    score = 0.0
    score += 2.0 * oechem.OECount(mol, oechem.OEAtomIsInRing())
    score += 1.0 * oechem.OECount(mol, oechem.OENotAtom(oechem.OEAtomIsInRing()))

    return score


InterfaceData = '''
!CATEGORY "input options"

    !PARAMETER -in
      !ALIAS -i
      !TYPE string
      !REQUIRED true
      !KEYLESS 1
      !VISIBILITY simple
      !BRIEF Input filename
    !END

    !PARAMETER -minbonds
      !ALIAS -minb
      !TYPE int
      !REQUIRED false
      !VISIBILITY simple
      !LEGAL_RANGE 1 10
      !DEFAULT 5
      !BRIEF Minimum number of bonds in core
    !END

    !PARAMETER -maxbonds
      !ALIAS -maxb
      !TYPE int
      !REQUIRED false
      !VISIBILITY simple
      !LEGAL_RANGE 4 18
      !DEFAULT 18
      !BRIEF Maximun number of bonds in core
    !END

!END
'''

if __name__ == "__main__":
    sys.exit(main(sys.argv))
