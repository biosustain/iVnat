"""
Mapping reaction Ids from seed to metanetx then to metacyc
"""

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




### Reading the reac_prop

#reac_prop = read_csv("/Users/lizrad/Dev/iVnat/reac_prop.tsv" , sep='\t', skiprows=365)
#reac_prop_reordered = reac_prop.set_index("#MNX_ID")
#reac_prop_dict = reac_prop_reordered.to_dict('index')

### Reading the reac_xref

reac_xref = read_csv("/Users/lizrad/Documents/Vibrio_folder/reac_xref.tsv" , sep='\t')
reac_xref

### Number of seed reactions in the model

print ('This is the total amount of reactions in the model: %i' % (len(model.reactions)))

print ('We focus on mapping those %i reactions in the model that have an SEED ID' % (len([met for met in model.reactions if met.id.startswith('rxn') ])))

extracell=[]
for met in model.reactions:
    if not met.id.startswith('rxn'):
        print(met.id)
        extracell.append(met.id)
len(extracell)        

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


#reac_property_MNX = read_csv("/Users/lizrad/Documents/Vibrio_folder/reac_xref.tsv" , sep='\t')
#reac_property_MNX.set_index('MNX_ID',inplace = True)
#df_3 = reac_property_MNX.fillna('MISSING')  #Fill NA/NaN values using the specified methodss
#df_3
#reac_property_MNX_dict = df_3.to_dict('index')
#reac_property_MNX_dict

#skip
#rxn_ID_not_in_MetaNetX_2 = [rxn for rxn in rxn_nobiggid if rxn.id.startswith('MNXR') and ''.join(re.findall('(MNXR\d*)_',rxn.id)) not in reac_xref_MNX_2_dict.keys()]
#print (len(rxn_ID_not_in_MetaNetX_2))

# Get a mapping from KEGG to the new MNXR IDs, for those reactions that cannot be found in the MNX Database,
# but that have a KEGG ID in their annotation. 
#153644- 170637
df_2 = reac_xref.loc[153642:170637,'XREF_ID':'MNX_ID']
df_2 = df_2.set_index(['XREF_ID'])
seed_2_mnx_dict_2 = df_2.to_dict('index')
seed_2_mnx_dict_2

### Creating seed dict

reac_xref_MNX_seed = read_csv("/Users/lizrad/Documents/Vibrio_folder/reac_xref.tsv" , sep='\t', low_memory=False)
df_seed = reac_xref_MNX_seed.loc[153642:170635,'XREF':'MNX_ID']
pivoted_seedd = df_seed.pivot_table(index='XREF',columns='XREF_ID',values='MNX_ID',aggfunc = lambda x: x)
#pivoted_seedd=pivoted_seed.drop(['metanetx'])
pivoted_seedd
reac_xref_MNX_seed_dict = pivoted_seedd.to_dict('index')
reac_xref_MNX_seed_dict2 = reac_xref_MNX_seed_dict['seed']
reac_xref_MNX_seed_dict2

more_seed={}
for seed in reac_xref_MNX_seed_dict2.keys():
        
    if  reac_xref_MNX_seed_dict2[seed] not in more_seed.keys():
        more_seed.setdefault(reac_xref_MNX_seed_dict2[seed], [seed])
    else:
        more_seed[reac_xref_MNX_seed_dict2[seed]].append(seed)

for ls in more_seed.values():
    if len(ls) > 1:
        print(ls)

### Creating kegg dict 


#26958:46003
df_kegg = reac_xref.loc[26956:46002,'XREF':'MNX_ID']
pivoted_kegg = df_kegg.pivot_table(index='XREF',columns='XREF_ID',values='MNX_ID',aggfunc = lambda x: x)
pivoted_kegg=pivoted_kegg.drop(['metanetx'])
pivoted_kegg
reac_xref_MNX_kegg_dict = pivoted_kegg.to_dict('index')
reac_xref_MNX_kegg_dict2 = reac_xref_MNX_kegg_dict['kegg']

#for reac in reac_xref_MNX_kegg_dict2.keys():
 #   if reac.startswith("MNXR"):
  #      reac_xref_MNX_kegg_dict2.pop(reac)
        
        
print(reac_xref_MNX_kegg_dict2)

more_kegg={}
for kegg in reac_xref_MNX_kegg_dict2.keys():
        
    if  reac_xref_MNX_kegg_dict2[kegg] not in more_kegg.keys():
        more_kegg.setdefault(reac_xref_MNX_kegg_dict2[kegg], [kegg])
    else:
        more_kegg[reac_xref_MNX_kegg_dict2[kegg]].append(kegg)
    

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


for ls in more_kegg.values():
    if len(ls) > 1:
        print(ls)

### Creating Bigg dict

df_bigg = reac_xref.loc[4:26955,'XREF':'MNX_ID']
pivoted_bigg = df_bigg.pivot_table(index='XREF',columns='XREF_ID',values='MNX_ID',aggfunc = lambda x: x)
#pivoted_seedd=pivoted_seed.drop(['metanetx'])
pivoted_bigg
reac_xref_MNX_bigg_dict = pivoted_bigg.to_dict('index')
reac_xref_MNX_bigg_dict2 = reac_xref_MNX_bigg_dict['bigg']
print(reac_xref_MNX_bigg_dict2)



more_bigg={}
for bigg in reac_xref_MNX_bigg_dict2.keys():
        
    if  reac_xref_MNX_bigg_dict2[bigg] not in more_bigg.keys():
        more_bigg.setdefault(reac_xref_MNX_bigg_dict2[bigg], [bigg])
    else:
        more_bigg[reac_xref_MNX_bigg_dict2[bigg]].append(bigg)

for ls in more_bigg.values():
    if len(ls) > 1:
        print(ls)

### Creating metacyc dict

df_metacyc = reac_xref.loc[46003:70388,'XREF':'MNX_ID']
pivoted_metacyc = df_metacyc.pivot_table(index='XREF',columns='XREF_ID',values='MNX_ID',aggfunc = lambda x: x)
pivoted_metacyc
reac_xref_MNX_metacyc_dict = pivoted_metacyc.to_dict('index')
reac_xref_MNX_metacyc_dict2 = reac_xref_MNX_metacyc_dict['metacyc']
print(reac_xref_MNX_metacyc_dict2)

more_metacyc={}
for metacyc in reac_xref_MNX_metacyc_dict2.keys():
        
    if  reac_xref_MNX_metacyc_dict2[metacyc] not in more_metacyc.keys():
        more_metacyc.setdefault(reac_xref_MNX_metacyc_dict2[metacyc], [metacyc])
    else:
        more_metacyc[reac_xref_MNX_metacyc_dict2[metacyc]].append(metacyc)

for ls in more_metacyc.values():
    if len(ls) > 1:
        print(ls)

### SEED Ids to reaction annotation

for reac in model.reactions:
    if reac.id.startswith('rxn'):
        split, compartment=reac.id.split("_")
        reac.annotation['seed.compound']=split
        print(reac.annotation)
        compartment
        



### MNXR Ids to reaction annotation



for reac in model.reactions:
    data=''
    if 'seed.compound' in reac.annotation:
        data=reac.annotation['seed.compound'] 
        
        for xref in reac_xref_MNX_seed_dict2.keys():
            if xref in data:
                reac.annotation['metanetx.reaction']=reac_xref_MNX_seed_dict2[xref]
                print(reac.annotation)

### Metacyc Ids to reaction annotation

for reac in model.reactions: 
    for metacyc in reac_xref_MNX_metacyc_dict2.keys():
        if "metanetx.reaction" not in reac.annotation:
            continue
        if reac_xref_MNX_metacyc_dict2[metacyc]==reac.annotation["metanetx.reaction"]:
            if "biocyc" not in reac.annotation :
                reac.annotation["biocyc"] = []
                    
                reac.annotation["biocyc"].append(metacyc)
            elif metacyc not in met.annotation["biocyc"]:
                reac.annotation["biocyc"].append(metacyc)
            else:
                continue
            
            print(reac.annotation)

### Bigg Ids to reaction annotation

for reac in model.reactions: 
    for bigg in reac_xref_MNX_bigg_dict2.keys():
        if "metanetx.reaction" not in reac.annotation:
            continue
        if reac_xref_MNX_bigg_dict2[bigg]==reac.annotation["metanetx.reaction"]:
            if "bigg.reaction" not in reac.annotation :
                reac.annotation["bigg.reaction"] = []
                    
                reac.annotation["bigg.reaction"].append(bigg)
            elif bigg not in reac.annotation["bigg.reaction"]:
                reac.annotation["bigg.reaction"].append(bigg)
            else:
                continue
            
            print(reac.annotation)

### Kegg Ids to reaction annotation

for reac in model.reactions: 
    for kegg in reac_xref_MNX_kegg_dict2.keys():
        if "metanetx.reaction" not in reac.annotation:
            continue
        if reac_xref_MNX_kegg_dict2[kegg]==reac.annotation["metanetx.reaction"]:
            if "kegg.reaction" not in reac.annotation :
                reac.annotation["kegg.reaction"] = []
                    
                reac.annotation["kegg.reaction"].append(kegg)
            elif kegg not in reac.annotation["kegg.reaction"]:
                reac.annotation["kegg.reaction"].append(kegg)
            else:
                continue
            
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

model.reactions.rxn03023_c0.reactants
model.reactions.rxn03932_c0.reactants

model.reactions.rxn03023_c0.annotation

model.reactions.rxn03932_c0.reaction


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
            
            print(reac.id)
            model.repair()





for reac in model.reactions:
    
    
    if "biocyc" not in reac.annotation.keys():
        #print(reac.id)
        continue
    if reac.id.startswith('rxn'):    
        if reac.id in ["rxn03023_c0"]:
            pass
        else:
            split, compartment=reac.id.split("_")
            reac.id
            reac.id=reac.annotation["biocyc"][0]+ '_' + compartment           
            print(reac.id)
            model.repair()