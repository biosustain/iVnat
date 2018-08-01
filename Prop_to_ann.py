from pandas import read_csv, DataFrame
import pandas as pd
import os
import re #Both patterns and strings to be searched can be Unicode strings as well as 8-bit strings.
import math
import cobra
import cobra.test
from __future__ import print_function
from os.path import join
from cobra.io import write_sbml_model

from IPython.core.interactiveshell import InteractiveShell

InteractiveShell.ast_node_interactivity = "all"

data_dir = "/Users/lizrad/Dev/iVnat"
print("files found: ")
print(", ".join(i for i in os.listdir(data_dir) if i.startswith("iVnat")))
model=cobra.io.read_sbml_model(join(data_dir, "iVnat.xml"))
#model=cobra.io.read_legacy_sbml("C:\\Users\Asus\Documents\Vibrio_project_literature\GCF_001456255.1_rast_metabolic_model.SBML\GCF_001456255.1_rast_metabolic_model.xml")
model



for met in model.metabolites:
    print(met.annotation, met.id)

### Read chem_prop

chem_prop = read_csv("/Users/lizrad/Documents/Vibrio_folder/chem_prop.tsv" , sep='\t' , low_memory=False)
chem_prop


### Reorder chem_prop to dict

chem_prop_reordered = chem_prop.set_index("MNX_ID")
chem_prop_reordered
df = chem_prop_reordered.fillna('MISSING')
chem_prop_dict = df.to_dict('index')
chem_prop_dict

for met in model.metabolites: 
    for prop in chem_prop_dict.keys():
        if prop==met.annotation["metanetx.chemical"]:
            met.annotation["inchi"]=chem_prop_dict[prop]["InChI"]
            met.annotation["inchikey"]=chem_prop_dict[prop]["InChIKey"]
            met.annotation["smiles"]=chem_prop_dict[prop]["SMILES"]
            #met.annotation["Description"]=chem_prop_dict[prop]["Description"]
            met.formula=chem_prop_dict[prop]["Formula"]
            model.repair()
            print(met.annotation)
            print(met.id,met.formula)
            

diff_charge=[]
same_charge=[]
for met in model.metabolites: 
    for prop in chem_prop_dict.keys():
        if prop==met.annotation["metanetx.chemical"]:
            if met.charge==chem_prop_dict[prop]["Charge"]:
                same_charge.append(met.id)
                print("same", met.id,met.charge,chem_prop_dict[prop]["Charge"] )
            else:
                diff_charge.append(met.id)
                print("different", met.id,met.charge,chem_prop_dict[prop]["Charge"] )



len(diff_charge)
len(same_charge)

Conclusion: Do not change the seed charges. Stick to the resource in which the model was generated.

print(len(model.metabolites))

### Select the chebi sources and add to annotation

#i=0
for met in model.metabolites: 
    for prop in chem_prop_dict.keys():
        if prop==met.annotation["metanetx.chemical"]:
            name, number=chem_prop_dict[prop]["Source"].split(":")
            if name=="chebi":
                met.annotation["chebi"]="CHEBI:"+number
                #i=i+1
                print(met.annotation)
                model.repair()
                
#i
                
            

### Number of metabolites without Bigg ID

not_in_bigg=[]
for met in model.metabolites:
    if "bigg.metabolite" not in met.annotation:
        not_in_bigg.append(met.id)
len(not_in_bigg)        

for met in model.metabolites:
    print(met.annotation, met.id)

model.metabolites.NAD_c0.annotation

write_sbml_model(model, "/Users/lizrad/Dev/iVnat/iVnat.xml")