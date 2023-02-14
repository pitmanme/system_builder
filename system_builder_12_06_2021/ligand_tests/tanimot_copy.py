from openeye.oechem import *
from openeye.oeiupac import *
from openeye.oeomega import *
from openeye.oeshape import *
from openeye.oedepict import *
from openeye import oegraphsim
from openeye import oechem
from openeye import oedepict
from openeye import oegrapheme
from openeye import oeshape
from openeye import oegraphsim
import sys

# This should be changed to either read from files or
# to also be able to read from molecules passed in.
def main(argv=[__name__]):

    itf = oechem.OEInterface(InterfaceData)

    if not oechem.OEParseCommandLine(itf, argv):
        return 1

    fname = itf.GetString("-in")

    # Check input file.
    ifs = oechem.oemolistream()
    if not ifs.open(fname):
        oechem.OEThrow.Fatal("Cannot open input file!")
        
    tname = itf.GetString("-ref")
    
    # Check input file.
    ifs2 = oechem.oemolistream()
    if not ifs2.open(tname):
        oechem.OEThrow.Fatal("Cannot open input file!")

    # Check input file.
    ifs2 = oechem.oemolistream()
    if not ifs2.open(tname):
        oechem.OEThrow.Fatal("Cannot open reference file!")

    oname = itf.GetString("-out")
    
    # Reads in the mol you will overlay 'query' (fname),
    # the target mol, and will output the overlayed query.
    bestcore = shape_color_overlay(fname, tname, oname)
    visualize_overlap('test.png')


# check input/output files
def shape_color_overlay(fname, tname, oname):
    ''' input (str):
    fitting mol name, fname,
    target name, tname
    output file name, oname (needs to be oeb or oeb.gz)
    '''
    fifs = oechem.oemolistream()
    if not fifs.open(fname):
        oechem.OEThrow.Fatal("Cannot open query input file!")

    tifs = oechem.oemolistream()
    if not tifs.open(tname):
        oechem.OEThrow.Fatal("Cannot open target input file!")

    ofs = oechem.oemolostream()
    if not ofs.open(oname):
        oechem.OEThrow.Fatal("Cannot open output file!")
        
    # read molecules
    fmol = oechem.OEMol()
    if not oechem.OEReadMolecule(fifs, fmol):
        oechem.OEThrow.Fatal("Cannot read query molecule!")
    tmol = oechem.OEMol()
    if not oechem.OEReadMolecule(tifs, tmol):
        oechem.OEThrow.Fatal("Cannot read target molecule!")
    fifs.close()
    tifs.close()
    
    # Generate 3D confs
    omega = OEOmega()
    nconfs = 100
    omega.SetMaxConfs(nconfs)
    omega.SetStrictStereo(False) #false = random stereoisomer if unk
    omega.SetStrictAtomTypes(False)
    omega(fmol)
    #omega(tmol)
    
    # Perform shape overlay
    # Setup ROCS to provide specified number of conformers per hit
    options = OEROCSOptions()
    options.SetNumBestHits(1)
    options.SetConfsPerHit(nconfs)
    rocs = OEROCS(options)
    rocs.AddMolecule(fmol)
    
    # Loop over results and output
    tanimotos = []
    for r in rocs.Overlay(tmol):
        outmol = r.GetOverlayConf() #Use GetOverlayConf to get just the best; GetOverlayConfs for all
        OERemoveColorAtoms(outmol)
        OEAddExplicitHydrogens(outmol)
        OEWriteMolecule(ofs, outmol)
        score = r.GetTanimotoCombo()
        print("title: %s  tanimoto combo = %.2f" % (outmol.GetTitle(), r.GetTanimotoCombo()))
        tanimotos.append(score)
    ofs.close()
    return tanimotos, tmol, outmol
    
def visualize_overlap(outfile):
        
    ofs = oechem.oeofstream()
    if not ofs.open(outfile):
        oechem.OEThrow.Fatal("Cannot open output file!")
        
    # create image
    image = oedepict.OEImage(500, 500)
    
    # setup depiction options
    opts = oedepict.OE2DMolDisplayOptions(500, 500, oedepict.OEScale_AutoScale)
    oedepict.OESetup2DMolDisplayOptions(opts, itf)
    opts.SetBondWidthScaling(True)
    
InterfaceData = '''
!CATEGORY "input options"

    !PARAMETER -in
      !ALIAS -i
      !TYPE string
      !REQUIRED true
      !KEYLESS 1
      !VISIBILITY simple
      !BRIEF Queury ligand input filename
    !END

    !PARAMETER -ref
      !ALIAS -r
      !TYPE string
      !REQUIRED true
      !VISIBILITY simple
      !BRIEF Reference ligand filename
    !END

    !PARAMETER -out
      !ALIAS -o
      !TYPE string
      !REQUIRED false
      !VISIBILITY simple
      !DEFAULT out.sdf
      !BRIEF Maximun number of bonds in core
    !END

!END
'''

if __name__ == "__main__":
    sys.exit(main(sys.argv))
    
#tanimotos, tmol, outmol = shape_color_overlay('/Users/mpitman/work/system_tests/Naf/poulos_lead1.pdb', '/Users/mpitman/work/system_tests/Naf/6xk3_G_V4V.sdf', 'NOS_fitted1.sdf')
#visualize_overlap(tanimotos, tmol, outmol, 'vis.png')




