#This shows how to inegrate ligpose.py into my other scripts


import ligpose
# ----------------------------------------------------------------------
# Here for testing.
#pdb_path = '/Users/mpitman/work/system_tests/7cmd/7cmd.pdb'
#mtz_path = '/Users/mpitman/work/system_tests/7cmd/7cmd_phases.mtz'
#lig_in = '/Users/mpitman/work/system_tests/7cmd/lead_from_shapefit.sdf'
method = 'shapefit'

# NOS system

pdb_path = '/Users/mpitman/work/system_tests/Naf/shapefit_test/6xk3.pdb'
mtz_path = '/Users/mpitman/work/system_tests/Naf/shapefit_test/6xk3_phases.mtz'
lig_in = '/Users/mpitman/work/system_tests/Naf/poulos_lead1.pdb'

# ----------------------------------------------------------------------
ligpose.main(pdb_path, mtz_path, lig_in, method)
