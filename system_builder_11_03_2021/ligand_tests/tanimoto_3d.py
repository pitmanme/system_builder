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
    
def visualize_overlap(tanimotos, tmol, outmol, outfile):
    ofs = oechem.oeofstream()
    if not ofs.open(outfile):
        oechem.OEThrow.Fatal("Cannot open output file!")
        
    # create image
    image = oedepict.OEImage(500, 500)
    
    # setup depiction options
    opts = oedepict.OE2DMolDisplayOptions(500, 500, oedepict.OEScale_AutoScale)
    oedepict.OESetup2DMolDisplayOptions(opts)
    opts.SetBondWidthScaling(True)
    

tanimotos, tmol, outmol = shape_color_overlay('/Users/mpitman/work/system_tests/Naf/poulos_lead1.pdb', '/Users/mpitman/work/system_tests/Naf/6xk3_G_V4V.sdf', 'NOS_fitted1.sdf')




