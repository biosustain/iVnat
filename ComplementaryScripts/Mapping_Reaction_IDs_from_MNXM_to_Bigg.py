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

#skip
#data_dir = "/Users/lizrad/Dev/iVnat"
#print("files found: ")
#print(", ".join(i for i in os.listdir(data_dir) if i.startswith("iVnat")))
#model=cobra.io.read_sbml_model(join(data_dir, "iVnat_mapped.xml"))
#model=cobra.io.read_legacy_sbml("C:\\Users\Asus\Documents\Vibrio_project_literature\GCF_001456255.1_rast_metabolic_model.SBML\GCF_001456255.1_rast_metabolic_model.xml")
#model




### Reading the reac_xref

reac_xref = read_csv("/Users/lizrad/Documents/Vibrio_folder/reac_xref.tsv" , sep='\t')
reac_xref

### Number of MNXR reactions in the model

print ('This is the total amount of reactions in the model: %i' % (len(model.reactions)))

print ('We focus on mapping those %i reactions in the model that have an MNXR ID' % (len([met for met in model.reactions if met.id.startswith('MNXR') ])))

for met in model.reactions:
    if not met.id.startswith('MNXR'):
        print(met.id)

bigg_reactions = [rxn for rxn in model.reactions if rxn.id.startswith('MNXR')]

print (len(bigg_reactions))
print (len(model.reactions))

#skip
rxn_mnx2bigg = []
rxn_nobiggid = []

for rxn in [rxn for rxn in model.reactions if rxn.id.startswith('MNXR')]:
    if 'bigg' in rxn.annotation:
        reaction_id = rxn.annotation['BIGG'][0].split(', ')
        tag_split = rxn.id.split('_')
        rxn.id = str(reaction_id[0]) +'_'+ tag_split[1]
        rxn_mnx2bigg.append(reaction_id)
        cameo_model.repair()
    else:
        rxn_nobiggid.append(rxn)
        
print len(rxn_mnx2bigg)
print len(rxn_nobiggid)

### Formatting and creating Bigg MNX dictionary


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

#skip
#rxn_ID_not_in_MetaNetX_2 = [rxn for rxn in rxn_nobiggid if rxn.id.startswith('MNXR') and ''.join(re.findall('(MNXR\d*)_',rxn.id)) not in reac_xref_MNX_2_dict.keys()]
#print len(rxn_ID_not_in_MetaNetX_2)

# Get a mapping from KEGG to the new MNXR IDs, for those reactions that cannot be found in the MNX Database,
# but that have a KEGG ID in their annotation. 
#153644- 170637
df_2 = reac_xref.loc[6:26957,'XREF_ID':'MNX_ID']
df_2 = df_2.set_index(['XREF_ID'])
seed_2_mnx_dict_2 = df_2.to_dict('index')
seed_2_mnx_dict_2

reac_xref_MNXR_bigg = read_csv("/Users/lizrad/Documents/Vibrio_folder/reac_xref.tsv" , sep='\t', low_memory=False)
df_bigg = reac_xref_MNXR_bigg.loc[4:26955,'XREF':'MNX_ID']
pivoted_bigg = df_bigg.pivot_table(index='XREF',columns='XREF_ID',values='MNX_ID',aggfunc = lambda x: x)
pivoted_bigg=pivoted_bigg.drop(['metanetx'])
pivoted_bigg
reac_xref_MNX_bigg_dict = pivoted_bigg.to_dict('index')
reac_xref_MNX_bigg_dict2 = reac_xref_MNX_bigg_dict['bigg']
reac_xref_MNX_bigg_dict2

not_in_reac_xref = []
in_reac_xref = []

for met in model.reactions:
    split=met.id.split("_")
    met_clean=split[0]
    print(met_clean)
    if not met_clean == '':
        
        if not met_clean in reac_xref_MNX_bigg_dict2.keys():
            not_in_reac_xref.append(met)
        else:
            in_reac_xref.append(met)

print ('%i metabolites with MNXR ID cannot be found in the MNX database' % (len(not_in_reac_xref)))
in_reac_xref

### Bigg Ids to reaction annotation

for reac in model.reactions:
    reac.annotation

for reac in model.reactions:
    data=''
    if 'metanetx.reaction' in reac.annotation:
        data=reac.annotation['metanetx.reaction'] 
        
        for xref in reac_xref_MNX_bigg_dict2.keys():
            if reac_xref_MNX_bigg_dict2[xref]==data:
                reac.annotation['bigg.reaction']= xref
                print(reac.annotation)
                
                


### Check for multiple ids 

dupl_bigg = {}
for reac in model.reactions:
    if "bigg.reaction" not in reac.annotation.keys():  #skip EX reactions
        #print(reac.id)
        continue
    if reac.annotation["bigg.reaction"] not in dupl_bigg.keys():
        dupl_bigg.setdefault(reac.annotation["bigg.reaction"],[reac.id]) #similar to get(), but will set dict[key]=default if key is not already in dict.


    else:
        dupl_bigg[reac.annotation["bigg.reaction"]].append(reac.id)
        
dupl_bigg


for ls in dupl_bigg.values():
    if len(ls) > 1:
        print(ls)

for ls in dupl_bigg.values():
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
    
    
    if "bigg.reaction" not in reac.annotation.keys():
        #print(reac.id)
        continue
    if reac.id.startswith('MNXR'):    
        #if reac.id in ["rxn03023_c0"]:
         #   pass
        #else:
        split, compartment=reac.id.split("_")
        reac.id
        reac.id=reac.annotation['bigg.reaction']+ '_' + compartment
        model.repair()
        print(reac.id)
    

#chem_xref['#XREF'] = chem_xref['#XREF'].str.replace(':','\t')

#^[0-9A-Za-z]*:[0-9A-Za-z]*
#    (^[0-9A-Za-z]*):([0-9A-Za-z]*)
#        $1 \t $2
        



write_sbml_model(model, "/Users/lizrad/Dev/iVnat/iVnat.xml")

model.reactions.EX_cpd00001_e0.metabolites


