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

#skip
relative_directory = os.getcwd()
print(relative_directory)

data_dir = "/Users/lizrad/Documents/Vibrio_folder"
print("files found: ")
print(", ".join(i for i in os.listdir(data_dir) if i.startswith("GCF")))
model_orig=cobra.io.read_legacy_sbml(join(data_dir, "GCF_001456255.1_rast_metabolic_model.xml"))
#model=cobra.io.read_legacy_sbml("C:\\Users\Asus\Documents\Vibrio_project_literature\GCF_001456255.1_rast_metabolic_model.SBML\GCF_001456255.1_rast_metabolic_model.xml")
model_orig


data_dir = "/Users/lizrad/Dev/iVnat"
print("files found: ")
print(", ".join(i for i in os.listdir(data_dir) if i.startswith("iVnat")))
model=cobra.io.read_sbml_model(join(data_dir, "iVnat.xml"))
#model=cobra.io.read_legacy_sbml("C:\\Users\Asus\Documents\Vibrio_project_literature\GCF_001456255.1_rast_metabolic_model.SBML\GCF_001456255.1_rast_metabolic_model.xml")
model

print("files found: ")
reac_xref = read_csv("/Users/lizrad/Documents/Vibrio_folder/reac_xref.tsv" , sep='\t')
reac_xref

### Number of seed reactions in the model

print ('This is the total amount of reactions in the model: %i' % (len(model.reactions)))

print ('We focus on mapping those %i reactions in the model that have an SEED ID' % (len([met for met in model.reactions if met.id.startswith('rxn') ])))

for met in model.reactions:
    if not met.id.startswith('rxn'):
        print(met.id)

seed_reactions = [rxn for rxn in model.reactions if rxn.id.startswith('rxn')]

print (len(seed_reactions))
print (len(model.reactions))


### Formatting and creating SEED MNX dictionary


df_1 = reac_xref.loc[:,'XREF':'MNX_ID']
groups = df_1.groupby(['MNX_ID','XREF'])
df_1 = groups.apply(lambda x:list(x['XREF_ID']))
df_1 =df_1.unstack('XREF') 
reac_xref_MNX_1_dict = df_1.to_dict('index')
reac_xref_MNX_1_dict


reac_property_MNX = read_csv("/Users/lizrad/Documents/Vibrio_folder/reac_xref.tsv" , sep='\t')
reac_property_MNX.set_index('MNX_ID',inplace = True)
df_3 = reac_property_MNX.fillna('MISSING')  #Fill NA/NaN values using the specified method
df_3
reac_property_MNX_dict = df_3.to_dict('index')
reac_property_MNX_dict



#153644- 170637
df_2 = reac_xref.loc[153642:170637,'XREF_ID':'MNX_ID']
df_2 = df_2.set_index(['XREF_ID'])
seed_2_mnx_dict_2 = df_2.to_dict('index')
seed_2_mnx_dict_2

reac_xref_MNX_seed = read_csv("/Users/lizrad/Documents/Vibrio_folder/reac_xref.tsv" , sep='\t', low_memory=False)
df_seed = reac_xref_MNX_seed.loc[153642:170635,'XREF':'MNX_ID']
pivoted_seedd = df_seed.pivot_table(index='XREF',columns='XREF_ID',values='MNX_ID',aggfunc = lambda x: x)
#pivoted_seedd=pivoted_seed.drop(['metanetx'])
pivoted_seedd
reac_xref_MNX_seed_dict = pivoted_seedd.to_dict('index')
reac_xref_MNX_seed_dict2 = reac_xref_MNX_seed_dict['seed']
reac_xref_MNX_seed_dict2

not_in_reac_xref = []
in_reac_xref = []

for met in model.reactions:
    split=met.id.split("_")
    met_clean=split[0]
    print(met_clean)
    if not met_clean == '':
        
        if not met_clean in reac_xref_MNX_seed_dict2.keys():
            not_in_reac_xref.append(met)
        else:
            in_reac_xref.append(met)

print ('%i metabolites with SEED ID cannot be found in the MNX database' % (len(not_in_reac_xref)))


### SEED Ids to reaction annotation

for reac in model.reactions:
    if reac.id.startswith('rxn'):
        split, compartment=reac.id.split("_")
        reac.annotation['seed.compound']=split
        print(reac.annotation)
        compartment
        

for reac in model.reactions:
    data=''
    if 'seed.compound' in reac.annotation:
        data=reac.annotation['seed.compound'] 
        
        for xref in reac_xref_MNX_seed_dict2.keys():
            if xref in data:
                reac.annotation['metanetx.reaction']=reac_xref_MNX_seed_dict2[xref]
                print(reac.annotation)

### Check for multiple ids 

dupl_seed = {}
for reac in model.reactions:
    if "seed.compound" not in reac.annotation.keys():  #skip EX reactions
        #print(reac.id)
        continue
    if reac.annotation["seed.compound"] not in dupl_seed.keys():
        dupl_seed.setdefault(reac.annotation["seed.compound"],[reac.id]) #similar to get(), but will set dict[key]=default if key is not already in dict.


    else:
        dupl_seed[reac.annotation["seed.compound"]].append(reac.id)
        
dupl_seed  


model.reactions.EX_cpd00067_e0.annotation

dupl_mnxr = {}
for reac in model.reactions:
    if "metanetx.reaction" not in reac.annotation.keys():
        #print(reac.id)
        continue
        
    if reac.annotation["metanetx.reaction"] not in dupl_mnxr.keys():
            dupl_mnxr.setdefault(reac.annotation["metanetx.reaction"],[reac.id])
    else:
            dupl_mnxr[reac.annotation["metanetx.reaction"]].append(reac.id)
            
           

for ls in dupl_mnxr.values():
    if len(ls) > 1:
        print(ls)

for ls in dupl_mnxr.values():
    if len(ls) > 1:
        print(ls)
        x=[]
        y=[]
        for members in ls: 
            main_id, comp =members.split("_")
            x.append(main_id)
            y.append(comp)
        if not x[0]==x[1]:
            print(x,y)
            for pair in zip(x,y):
                print(len(model.reactions.get_by_id(pair[0]+"_"+pair[1]).reaction))




for reac in model.reactions:
    
    
    if "metanetx.reaction" not in reac.annotation.keys():
        #print(reac.id)
        continue
    if reac.id.startswith('rxn'):    
        if reac.id in ["rxn03023_c0"]:
            pass
        else:
            split, compartment=reac.id.split("_")
            reac.id
            reac.id=reac.annotation['metanetx.reaction']+ '_' + compartment
            model.repair()
            print(reac.id)
    




write_sbml_model(model, "/Users/lizrad/Dev/iVnat/iVnat.xml")


