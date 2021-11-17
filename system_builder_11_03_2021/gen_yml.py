import yaml


        
def generate_yml_files():

    """ Creates *.yml files for pmxworkflow RBFE calculations
    and deposit them to the set locations for pmxworkflow. """
    
"""
    ligands = {};
    for cellA, cellB in zip(SHEET['A'], SHEET['B']):
   #     ligands['---']
        ligands['name'] = '     ' + str(cellA.value),
        ligands['smiles'] = '   ' + str(cellB.value),
        ligands['docked'] = '   03_docked/{}/{}.sdf'.format(cellA.value, cellA.value)

        ligands.append(ligands)

    with open('ligands.yml', 'w') as outfile:
        yaml.dump(ligands, outfile)
       

Earlier version
  
    ligands = {}
    for cellA, cellB in zip(SHEET['A'], SHEET['B']):
        lig_entry = {
            '---'
            'name': '     ' + str(cellA.value),
            'smiles': '   ' + str(cellB.value),
            'docked': '   03_docked/{}/{}.sdf'.format(cellA.value, cellA.value),
            'measurement': 'N/A ',
            '    ki': 'N/A ',
            '    doi': 'N/A ',
            '    comment': 'N/A ',
            }

        ligands.append(lig_entry)

#    with open(r'.yml', 'w') as file:
    documents = yaml.dump(ligands, file)
  """
  
 # ^the above hasn't yet been successful. I will need to figure out how to best generate the yml files in a scriptable way. It is a little lower priority
    
#IF you want to just call this script as one thing and have it blindly do its job
# insert the if main function and just string everything together
