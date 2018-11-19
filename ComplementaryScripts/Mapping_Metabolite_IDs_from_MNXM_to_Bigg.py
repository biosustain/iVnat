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

### Read chem_xref

chem_xref = read_csv("/Users/lizrad/Documents/Vibrio_folder/chem_xref.tsv" , sep='\t' , low_memory=False)
chem_xref

### Bigg IDs to dict key

df_bigg = chem_xref.loc[2:10216,'XREF':'MNX_ID']
pivoted_bigg = df_bigg.pivot_table(index='XREF',columns='XREF_ID',values='MNX_ID',aggfunc = lambda x: x)
pivoted_bigg
#pivoted_bigg=pivoted_bigg.drop(['metanetx'])
chem_xref_MNX_bigg_dict = pivoted_bigg.to_dict('index')
chem_xref_MNX_bigg_dict2 = chem_xref_MNX_bigg_dict['bigg']
chem_xref_MNX_bigg_dict2

### Bigg IDs to annotation field

for met in model.metabolites:
    data=''
    if 'metanetx.chemical' in met.annotation:
        data=met.annotation['metanetx.chemical'] 
        
        for xref in chem_xref_MNX_bigg_dict2.keys():
            if chem_xref_MNX_bigg_dict2[xref]==data:
                met.annotation['bigg.metabolite']=xref
                print(met.annotation)


### Check for multiplicatives among the Bigg IDs

dupl_bigg = {}
for met in model.metabolites:
    if "bigg.metabolite" not in met.annotation.keys():
        continue
    if met.id in ["cpd00261_c0", "cpd06227_c0", "cpd03572_c0", "cpd02446_c0","cpd02572_c0", "cpd01466_c0" ]:
        pass
    else:    
        if met.annotation["bigg.metabolite"] not in dupl_bigg.keys():
            dupl_bigg.setdefault(met.annotation["bigg.metabolite"],[met.id])
        else:
            dupl_bigg[met.annotation["bigg.metabolite"]].append(met.id)

for ls in dupl_bigg.values():
    if len(ls) > 1:
        print(ls)

for ls in dupl_bigg.values():
    if len(ls) > 1:
        #print(ls)
        x=[]
        y=[]
        for members in ls: 
            main_id, comp =members.split("_")
            x.append(main_id)
            y.append(comp)
        if not x[0]==x[1]:
            print(x,y)
            for pair in zip(x,y):
                print(len(model.metabolites.get_by_id(pair[0]+"_"+pair[1]).reactions))

### Change the MNXM IDs to Bigg IDs

for met in model.metabolites:
    
    if "bigg.metabolite" not in met.annotation.keys():
        continue
    if met.id in ["cpd00261_c0", "cpd06227_c0", "cpd03572_c0", "cpd02446_c0","cpd02572_c0", "cpd01466_c0" ]:
        pass
    else:    
    #if met.id.startswith('cpd'):
        split, compartment=met.id.split("_")
        met.id
    
        met.id=met.annotation['bigg.metabolite']+ '_' + compartment
        model.repair()
        print(met.id)
#Number of metabolites without Bigg ID
not_in_bigg=[]
for met in model.metabolites:
    if "bigg.metabolite" not in met.annotation:
        not_in_bigg.append(met.id)
len(not_in_bigg)  

### Save draft with mapped Bigg IDs

write_sbml_model(model, "/Users/lizrad/Dev/iVnat/iVnat.xml")