"""
Changing the metabolite identifiers from seed to metanetx then to metacyc

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

relative_directory = os.getcwd()
print(relative_directory)

data_dir = "/Users/lizrad/Documents/Vibrio_folder"
print("files found: ")
print(", ".join(i for i in os.listdir(data_dir) if i.startswith("GCF")))
model=cobra.io.read_legacy_sbml(join(data_dir, "GCF_001456255.1_rast_metabolic_model.xml"))
#model=cobra.io.read_legacy_sbml("C:\\Users\Asus\Documents\Vibrio_project_literature\GCF_001456255.1_rast_metabolic_model.SBML\GCF_001456255.1_rast_metabolic_model.xml")
model

# Metabolites

### Read chem_xref

chem_xref = read_csv("/Users/lizrad/Documents/Vibrio_folder/chem_xref.tsv" , sep='\t' , low_memory=False)
chem_xref


df_2 = chem_xref.loc[:,'XREF':'MNX_ID']
df_2
groups = df_2.groupby(['MNX_ID','XREF'])
df_2 = groups.apply(lambda x:list(x['XREF_ID']))
df_2
df_2 =df_2.unstack('XREF')
df_2
chem_xref_MNX_2_dict = df_2.to_dict('index')
chem_xref_MNX_2_dict




### SEED IDs to dict key


#499091 530559
chem_xref_MNX_seed = read_csv("/Users/lizrad/Documents/Vibrio_folder/chem_xref.tsv" , sep='\t', low_memory=False)
df_seed = chem_xref.loc[499088:530556,'XREF':'MNX_ID']
pivoted_seed = df_seed.pivot_table(index='XREF',columns='XREF_ID',values='MNX_ID',aggfunc = lambda x: x)
pivoted_seedd=pivoted_seed.drop(['metanetx'])
pivoted_seedd
chem_xref_MNX_seed_dict = pivoted_seedd.to_dict('index')
chem_xref_MNX_seed_dict2 = chem_xref_MNX_seed_dict['seed']


more_seed={}
for seed in chem_xref_MNX_seed_dict2.keys():
    if  chem_xref_MNX_seed_dict2[seed] not in more_seed.keys():
        more_seed.setdefault(chem_xref_MNX_seed_dict2[seed], [seed])
    else:
         more_seed[chem_xref_MNX_seed_dict2[seed]].append(seed)

for ls in more_seed.values():
    if len(ls) > 1:
        print(ls)

pivoted_seed.drop()

chem_xref_MNX_seed_dict2

### Bigg Ids to dict 

df_bigg = chem_xref.loc[2:10216,'XREF':'MNX_ID']
pivoted_bigg = df_bigg.pivot_table(index='XREF',columns='XREF_ID',values='MNX_ID',aggfunc = lambda x: x)
pivoted_bigg
#pivoted_bigg=pivoted_bigg.drop(['metanetx'])
chem_xref_MNX_bigg_dict = pivoted_bigg.to_dict('index')
chem_xref_MNX_bigg_dict2 = chem_xref_MNX_bigg_dict['bigg']
chem_xref_MNX_bigg_dict2

more_bigg={}
for bigg in chem_xref_MNX_bigg_dict2.keys():
    if  chem_xref_MNX_bigg_dict2[bigg] not in more_bigg.keys():
        more_bigg.setdefault(chem_xref_MNX_bigg_dict2[bigg], [bigg])
    else:
         more_bigg[chem_xref_MNX_bigg_dict2[bigg]].append(bigg)

for ls in more_bigg.values():
    if len(ls) > 1:
        print(ls)

### Kegg Ids to dict keys

#kegg  325803-384352

df_kegg = chem_xref.loc[325801:384349,'XREF':'MNX_ID']
pivoted_kegg = df_kegg.pivot_table(index='XREF',columns='XREF_ID',values='MNX_ID',aggfunc = lambda x: x)
#pivoted_seedd=pivoted_seed.drop(['metanetx'])
pivoted_kegg
chem_xref_MNX_kegg_dict = pivoted_kegg.to_dict('index')
chem_xref_MNX_kegg_dict2 = chem_xref_MNX_kegg_dict['kegg']
print(chem_xref_MNX_kegg_dict2)


more_kegg={}
for kegg in chem_xref_MNX_kegg_dict2.keys():
    if kegg.startswith("C"):
        
        if  chem_xref_MNX_kegg_dict2[kegg] not in more_kegg.keys():
            more_kegg.setdefault(chem_xref_MNX_kegg_dict2[kegg], [kegg])
        else:
            more_kegg[chem_xref_MNX_kegg_dict2[kegg]].append(kegg)
    
    

for ls in more_kegg.values():
    if len(ls) > 1:
        print(ls)

### Metacyc IDs to dict keys

#metacyc: 457263-481348 
chem_xref_MNX_kegg = read_csv("/Users/lizrad/Documents/Vibrio_folder/chem_xref.tsv" , sep='\t', low_memory=False)
df_metacyc = chem_xref.loc[457261:481346,'XREF':'MNX_ID']
pivoted_metacyc = df_metacyc.pivot_table(index='XREF',columns='XREF_ID',values='MNX_ID',aggfunc = lambda x: x)
#pivoted_seedd=pivoted_seed.drop(['metanetx'])
pivoted_metacyc
chem_xref_MNX_metacyc_dict = pivoted_metacyc.to_dict('index')
chem_xref_MNX_metacyc_dict2 = chem_xref_MNX_metacyc_dict['metacyc']
print(chem_xref_MNX_metacyc_dict2)


### More metacyc ids for one MNXM 

more_metacyc={}
for meta in chem_xref_MNX_metacyc_dict2.keys():
        
    if  chem_xref_MNX_metacyc_dict2[meta] not in more_metacyc.keys():
        unique=more_metacyc.setdefault(chem_xref_MNX_metacyc_dict2[meta], [meta])
    else:
        more_metacyc[chem_xref_MNX_metacyc_dict2[meta]].append(meta)

len([x for x in more_metacyc.values() if len(x) > 1])

for ls in more_metacyc.values():
    if len(ls) > 1:
        print(ls)



### Check how many Metabolites can be found in the SEED database and how many cannot be found at all. 

# Quick check-up to compare conditions before and after:
print (len(model.metabolites))
print (len([met for met in model.metabolites if met.formula == '']))  # ??? none of them has a formula
print (len([met for met in model.metabolites if met.formula == '' and met.id.startswith('cpd')]))

print (len([met for met in model.metabolites if not met.charge]))
print (len([met for met in model.metabolites if not met.charge and met.id.startswith('cpd')]))


for met in model.metabolites:
    if met.id.startswith('cpd'):
        print (met.id)

print ('This is the total amount of metabolites in the model: %i' % (len(model.metabolites)))

print ('We focus on mapping those %i metabolites in the model that have an seed ID' % (len([met for met in model.metabolites if met.id.startswith('cpd')])))

### How many SEED metabolites can we find?

not_in_chem_xref = []
in_chem_xref = []




for met in model.metabolites:
    split=met.id.split("_")
    met_clean=split[0]
    print(met_clean)
    if not met_clean == '':
        if not met_clean in chem_xref_MNX_seed_dict2.keys():
            not_in_chem_xref.append(met)
        else:
            in_chem_xref.append(met)

print ('%i metabolites with SEED ID cannot be found in the MNX database' % (len(not_in_chem_xref)))




### Seed ids to annotation field

for met in model.metabolites:
    split, compartment=met.id.split("_")
    met.annotation['seed.compound']=split
    print(met.annotation)
    compartment



### MNX ids to annotation field

for met in model.metabolites:
    data=''
    if 'seed.compound' in met.annotation:
        data=met.annotation['seed.compound'] 
        
        for xref in chem_xref_MNX_seed_dict2.keys():
            if xref==data:
                met.annotation['metanetx.chemical']=chem_xref_MNX_seed_dict2[xref]
                print(met.annotation)
          
                

### Bigg Ids to annotation field

#for met in model.metabolites:
#    data=''
#    if 'metanetx.chemical' in met.annotation:
#        data=met.annotation['metanetx.chemical'] 
        
#        for xref in chem_xref_MNX_bigg_dict2.keys():
#            if chem_xref_MNX_bigg_dict2[xref]==data:
#                met.annotation['bigg.metabolite']=xref
#                print(met.annotation)


                
for met in model.metabolites: 
    for bigg in chem_xref_MNX_bigg_dict2.keys():
        #print (kegg)
        if "metanetx.chemical" not in met.annotation:
            continue

        if chem_xref_MNX_bigg_dict2[bigg]== met.annotation["metanetx.chemical"]:
                    
            if "bigg.metabolite" not in met.annotation :
                met.annotation["bigg.metabolite"] = []
                met.annotation["bigg.metabolite"].append(bigg)
            elif bigg not in met.annotation["bigg.metabolite"]:
                met.annotation["bigg.metabolite"].append(bigg)
            else:
                continue
                    
            print(met.annotation) 

### Kegg Ids to annotation field

for met in model.metabolites: 
    for kegg in chem_xref_MNX_kegg_dict2.keys():
        #print (kegg)
        if "metanetx.chemical" not in met.annotation:
            continue
        if kegg.startswith("C"):
            if chem_xref_MNX_kegg_dict2[kegg]== met.annotation["metanetx.chemical"]:
                
                if "kegg.compound" not in met.annotation :
                    met.annotation["kegg.compound"] = []
                    met.annotation["kegg.compound"].append(kegg)
                elif kegg not in met.annotation["kegg.compound"]:
                    met.annotation["kegg.compound"].append(kegg)
                else:
                    continue
                    
                print(met.annotation)     
        
                   
                
        

for met in model.metabolites: 
    print(met.id, met.annotation)

not_in_kegg=[]
for met in model.metabolites:
    if "kegg.compound" not in met.annotation:
        not_in_kegg.append(met.id)
len(not_in_kegg)     

### Metacyc IDs to annotation field

for met in model.metabolites: 
    for metacyc in chem_xref_MNX_metacyc_dict2.keys():
        if "metanetx.chemical" not in met.annotation:
            continue
        if chem_xref_MNX_metacyc_dict2[metacyc]==met.annotation["metanetx.chemical"]:
            if "biocyc" not in met.annotation :
                met.annotation["biocyc"] = []
                    
                met.annotation["biocyc"].append(metacyc)
            elif metacyc not in met.annotation["biocyc"]:
                met.annotation["biocyc"].append(metacyc)
            else:
                continue
            
            print(met.annotation)
            

model.metabolites.cpd00049_c0.annotation



not_in_metacyc=[]
for met in model.metabolites:
    if "biocyc" not in met.annotation:
        not_in_metacyc.append(met.id)
len(not_in_metacyc)   

### Check for multiplicates among MNX IDs

dupl = {}
for met in model.metabolites:
    if met.annotation["metanetx.chemical"] not in dupl.keys():
        dupl.setdefault(met.annotation["metanetx.chemical"],[met.id])
    else:
        dupl[met.annotation["metanetx.chemical"]].append(met.id)

for ls in dupl.values():
    if len(ls) > 1:
        print(ls)

for ls in dupl.values():
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
            reaction_involved=[]
            for pair in zip(x,y):
                react_number=len(model.metabolites.get_by_id(pair[0]+"_"+pair[1]).reactions)
                reaction_involved.append(react_number)
                print (react_number)
            
        #if not y[0] == y[1]: 
        #   print("different compartment same number", x,y)

### Change the SEED IDs to MNXM IDs

 Metabolites with the following IDs are not mapped to MNXM ID: "cpd00261_c0", "cpd06227_c0", "cpd03572_c0", "cpd02446_c0","cpd02572_c0", "cpd01466_c0". They were selected by their number of involved reactions.  

for met in model.metabolites:
    split, compartment=met.id.split("_")
    if met.id.startswith('cpd'):
        met.id
    if met.id in ["cpd00261_c0", "cpd06227_c0", "cpd03572_c0", "cpd02446_c0","cpd02572_c0", "cpd01466_c0" ]:
        pass
    else:
        met.id=met.annotation['metanetx.chemical']+ '_' + compartment
        model.repair()
    print(met.id)
    



### Change the MNXM  id to Metacyc where you can

for met in model.metabolites:
    split, compartment=met.id.split("_")
    if met.id.startswith('MNXM'):
        if "biocyc" not in met.annotation:
                continue
        met.id
        if met.id in ["cpd00261_c0", "cpd06227_c0", "cpd03572_c0", "cpd02446_c0","cpd02572_c0", "cpd01466_c0" ]:
            pass
        else:
            met.id=met.annotation['biocyc'][0]+"_"+compartment
            model.repair()
        print(met.id)

### Remaining ids that have seed ids

for met in model.metabolites:
    if met.id.startswith("cpd"):
        print(met.id)
    
      
        


## Save draft with mapped metabolites

write_sbml_model(model, "/Users/lizrad/Dev/iVnat/iVnat.xml")

#Conclusion:   
#6 metabolites with seed ids, because of multiplicatives among MNXM ids, decided by number of reactions involved  
#474 with MNXM ids because they lack metacyc id

