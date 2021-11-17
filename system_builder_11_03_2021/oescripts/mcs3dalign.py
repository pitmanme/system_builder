#!/usr/bin/env python
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

#############################################################################
# Align two compounds based on the maximum common substructure
#############################################################################
import sys
from openeye import oechem


def MCSAlign(refmol, fitmol, ofs):
    atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_Aromaticity
    bondexpr = 0
    mcss = oechem.OEMCSSearch(oechem.OEMCSType_Exhaustive)
    mcss.Init(refmol, atomexpr, bondexpr)
    mcss.SetMCSFunc(oechem.OEMCSMaxBondsCompleteCycles())

    rmat = oechem.OEDoubleArray(9)
    trans = oechem.OEDoubleArray(3)
    unique = True
    overlay = True
    for match in mcss.Match(fitmol, unique):
        rms = oechem.OERMSD(mcss.GetPattern(), fitmol, match, overlay, rmat, trans)
        if rms < 0.0:
            oechem.OEThrow.Warning("RMS overlay failure")
            continue
        oechem.OERotate(fitmol, rmat)
        oechem.OETranslate(fitmol, trans)
        oechem.OEWriteMolecule(ofs, fitmol)


def main(argv=[__name__]):
    if len(argv) != 4:
        oechem.OEThrow.Usage("%s <refmol> <fitmol> <outfile>" % argv[0])

    reffs = oechem.oemolistream()
    if not reffs.open(argv[1]):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % argv[1])
    if not oechem.OEIs3DFormat(reffs.GetFormat()):
        oechem.OEThrow.Fatal("Invalid input format: need 3D coordinates")
    refmol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(reffs, refmol):
        oechem.OEThrow.Fatal("Unable to read molecule in %s" % argv[1])
    if not refmol.GetDimension() == 3:
        oechem.OEThrow.Fatal("%s doesn't have 3D coordinates" % refmol.GetTitle())

    fitfs = oechem.oemolistream()
    if not fitfs.open(argv[2]):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % argv[2])
    if not oechem.OEIs3DFormat(fitfs.GetFormat()):
        oechem.OEThrow.Fatal("Invalid input format: need 3D coordinates")

    ofs = oechem.oemolostream()
    if not ofs.open(argv[3]):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % argv[3])
    if not oechem.OEIs3DFormat(ofs.GetFormat()):
        oechem.OEThrow.Fatal("Invalid output format: need 3D coordinates")

    oechem.OEWriteConstMolecule(ofs, refmol)
    oechem.OESuppressHydrogens(refmol)

    for fitmol in fitfs.GetOEGraphMols():
        if not fitmol.GetDimension() == 3:
            oechem.OEThrow.Warning("%s doesn't have 3D coordinates" % fitmol.GetTitle())
            continue
        MCSAlign(refmol, fitmol, ofs)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
