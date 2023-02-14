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
# This version has minor edits: 
# usage: prompt > python3 simcalc2img.py -query query.mol -target target.mol -out similarity.png
#############################################################################
# Depicting molecule similarity
#############################################################################

import sys
from openeye import oechem
from openeye import oedepict
from openeye import oegrapheme
from openeye import oegraphsim


def main(argv=[__name__]):

    itf = oechem.OEInterface(InterfaceData)
    oedepict.OEConfigureImageWidth(itf, 1000.0)
    oedepict.OEConfigureImageHeight(itf, 500.0)
    oegraphsim.OEConfigureFingerPrint(itf, oegraphsim.OEGetFPType(oegraphsim.OEFPType_Tree))

    if not oechem.OEParseCommandLine(itf, argv):
        return 1

    qname = itf.GetString("-query")
    tname = itf.GetString("-target")
    oname = itf.GetString("-out")

    # check input/output files

    qifs = oechem.oemolistream()
    if not qifs.open(qname):
        oechem.OEThrow.Fatal("Cannot open query input file!")

    tifs = oechem.oemolistream()
    if not tifs.open(tname):
        oechem.OEThrow.Fatal("Cannot open target input file!")

    ext = oechem.OEGetFileExtension(oname)
    if not oedepict.OEIsRegisteredImageFile(ext):
        oechem.OEThrow.Fatal("Unknown image type!")

    ofs = oechem.oeofstream()
    if not ofs.open(oname):
        oechem.OEThrow.Fatal("Cannot open output file!")

    # read molecules

    qmol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(qifs, qmol):
        oechem.OEThrow.Fatal("Cannot read query molecule!")
    tmol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(tifs, tmol):
        oechem.OEThrow.Fatal("Cannot read target molecule!")

    # get fingerprint type

    fptype = oegraphsim.OESetupFingerPrint(itf)
    print("Using fingerprint type %s" % fptype.GetFPTypeString())

    # create image

    width, height = oedepict.OEGetImageWidth(itf), oedepict.OEGetImageHeight(itf)
    image = oedepict.OEImage(width, height)

    # setup depiction options

    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    oedepict.OESetup2DMolDisplayOptions(opts, itf)
    opts.SetBondWidthScaling(True)

    # depict molecules with overlaps

    DepictMoleculeOverlaps(image, qmol, tmol, fptype, opts)

    oedepict.OEWriteImage(ofs, ext, image)


def SetFingerPrintSimilarity(qmol, tmol, fptype, tag, maxvalue=0):

    qbonds = oechem.OEUIntArray(qmol.GetMaxBondIdx())
    tbonds = oechem.OEUIntArray(tmol.GetMaxBondIdx())

    for match in oegraphsim.OEGetFPOverlap(qmol, tmol, fptype):
        for bond in match.GetPatternBonds():
            qbonds[bond.GetIdx()] += 1
        for bond in match.GetTargetBonds():
            tbonds[bond.GetIdx()] += 1

    maxvalue = max(maxvalue, max(qbonds))
    maxvalue = max(maxvalue, max(tbonds))

    for bond in qmol.GetBonds():
        bond.SetData(tag, qbonds[bond.GetIdx()])
    for bond in tmol.GetBonds():
        bond.SetData(tag, tbonds[bond.GetIdx()])

    return maxvalue


class ColorBondByOverlapScore(oegrapheme.OEBondGlyphBase):
    def __init__(self, cg, tag):
        oegrapheme.OEBondGlyphBase.__init__(self)
        self.colorg = cg
        self.tag = tag

    def RenderGlyph(self, disp, bond):

        bdisp = disp.GetBondDisplay(bond)
        if bdisp is None or not bdisp.IsVisible():
            return False

        if not bond.HasData(self.tag):
            return False

        linewidth = disp.GetScale() / 3.0
        color = self.colorg.GetColorAt(bond.GetData(self.tag))
        pen = oedepict.OEPen(color, color, oedepict.OEFill_Off, linewidth)

        adispB = disp.GetAtomDisplay(bond.GetBgn())
        adispE = disp.GetAtomDisplay(bond.GetEnd())

        layer = disp.GetLayer(oedepict.OELayerPosition_Below)
        layer.DrawLine(adispB.GetCoords(), adispE.GetCoords(), pen)

        return True

    def ColorBondByOverlapScore(self):
        return ColorBondByOverlapScore(self.colorg, self.tag).__disown__()


def DepictMoleculeOverlaps(image, qmol, tmol, fptype, opts):

    tag = oechem.OEGetTag("fpoverlap")
    maxvalue = SetFingerPrintSimilarity(qmol, tmol, fptype, tag)

    colorg = oechem.OELinearColorGradient()
    colorg.AddStop(oechem.OEColorStop(0.0, oechem.OEGreenTint))
    colorg.AddStop(oechem.OEColorStop(0.5, oechem.OECyan))
    colorg.AddStop(oechem.OEColorStop(maxvalue, oechem.OEBlue))
    bondglyph = ColorBondByOverlapScore(colorg, tag)

    oedepict.OEPrepareDepiction(qmol)
    overlaps = oegraphsim.OEGetFPOverlap(qmol, tmol, fptype)
    oedepict.OEPrepareMultiAlignedDepiction(tmol, qmol, overlaps)

    grid = oedepict.OEImageGrid(image, 1, 2)
    grid.SetMargin(oedepict.OEMargin_Bottom, 10)
    opts.SetDimensions(grid.GetCellWidth(), grid.GetCellHeight(), oedepict.OEScale_AutoScale)
    opts.SetAtomColorStyle(oedepict.OEAtomColorStyle_WhiteMonochrome)

    molscale = min(oedepict.OEGetMoleculeScale(qmol, opts),
                   oedepict.OEGetMoleculeScale(tmol, opts))
    opts.SetScale(molscale)

    qdisp = oedepict.OE2DMolDisplay(qmol, opts)
    oegrapheme.OEAddGlyph(qdisp, bondglyph, oechem.IsTrueBond())
    oedepict.OERenderMolecule(grid.GetCell(1, 1), qdisp)

    tdisp = oedepict.OE2DMolDisplay(tmol, opts)
    oegrapheme.OEAddGlyph(tdisp, bondglyph, oechem.IsTrueBond())
    oedepict.OERenderMolecule(grid.GetCell(1, 2), tdisp)

    qfp = oegraphsim.OEFingerPrint()
    oegraphsim.OEMakeFP(qfp, qmol, fptype)

    tfp = oegraphsim.OEFingerPrint()
    oegraphsim.OEMakeFP(tfp, tmol, fptype)

    score = oegraphsim.OETanimoto(qfp, tfp)

    font = oedepict.OEFont(oedepict.OEFontFamily_Default, oedepict.OEFontStyle_Default, 28,
                           oedepict.OEAlignment_Center, oechem.OEBlack)
    center = oedepict.OE2DPoint(image.GetWidth() / 2.0, image.GetHeight() - 20)
    image.DrawText(center, "Tree Tanimoto score = %.3f" % score, font)


InterfaceData = '''
!CATEGORY "input/output options"

  !PARAMETER -query
    !ALIAS -q
    !TYPE string
    !REQUIRED true
    !VISIBILITY simple
    !BRIEF input query molecule filename
  !END

  !PARAMETER -target
    !ALIAS -t
    !TYPE string
    !REQUIRED true
    !VISIBILITY simple
    !BRIEF input target molecule filename
  !END

  !PARAMETER -out
    !ALIAS -o
    !TYPE string
    !REQUIRED true
    !VISIBILITY simple
    !BRIEF Output filename
  !END
!END
'''

if __name__ == "__main__":
    sys.exit(main(sys.argv))
