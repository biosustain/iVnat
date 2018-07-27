"""
Changing the metabolite identifiers from seed to metanetx

"""


from __future__ import print_function
from pandas import read_csv, DataFrame
import pandas as pd
import os
import re #Both patterns and strings to be searched can be Unicode strings as well as 8-bit strings.
import math
import cobra
import cobra.test
from os.path import join


relative_directory = os.getcwd()
print(relative_directory)

data_dir = "/Users/lizrad/Documents/Vibrio_folder"
print("files found: ")
print(", ".join(i for i in os.listdir(data_dir) if i.startswith("GCF")))
model=cobra.io.read_legacy_sbml(join(data_dir, "GCF_001456255.1_rast_metabolic_model.xml"))
#model=cobra.io.read_legacy_sbml("C:\\Users\Asus\Documents\Vibrio_project_literature\GCF_001456255.1_rast_metabolic_model.SBML\GCF_001456255.1_rast_metabolic_model.xml")
model

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


chem_xref_MNX_seed_dict2

### Check how many Metabolites can be found in the SEED database and how many cannot be found at all. 

# Quick check-up to compare conditions before and after:
print (len(model.metabolites))
print (len([met for met in model.metabolites if met.formula == '']))  # ??? none of them has a formula
print (len([met for met in model.metabolites if met.formula == '' and met.id.startswith('cpd')]))

print (len([met for met in model.metabolites if not met.charge]))
print (len([met for met in model.metabolites if not met.charge and met.id.startswith('cpd')]))


for met in model.metabolites:
    print(met.id.startswith('cpd'))

print ('This is the total amount of metabolites in the model: %i' % (len(model.metabolites)))

print ('We focus on mapping those %i metabolites in the model that have an seed ID' % (len([met for met in model.metabolites if met.id.startswith('cpd')])))

### How many SEED metabolites can we find?
#splitting the seed id to get rid of the additional compartment tag
for met in model.metabolites:
    id_split=met.id.split("_")
    print(id_split[0])
    

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
len(in_chem_xref)



### Seed ids to annotation field

#met_clean, compartment=model.metabolites.cpd00001_e0.id.split("_")
    

for met in model.metabolites:
    split, compartment=met.id.split("_")
    met.annotation['seed.compound']=split
    print(met.annotation)
    compartment



### MNX ids to annotation field

for met in model.metabolites:
    print(met.id)

for met in model.metabolites:
    data=''
    if 'seed.compound' in met.annotation:
        data=met.annotation['seed.compound'] 
        
        for xref in chem_xref_MNX_seed_dict2.keys():
            if xref==data:
                met.annotation['metanetx.chemical']=chem_xref_MNX_seed_dict2[xref]
                print(met.annotation)
              
            
                


#for met in model.metabolites:
#    if met.id.startswith('MNXM'):
#        met.id=met.annotation['seed.compound']
#       print(met.id)
    
        

### Check for multiplicates among MNXM IDs

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
            for pair in zip(x,y):
                print(len(model.metabolites.get_by_id(pair[0]+"_"+pair[1]).reactions))
        #if x[0]==x[1] and y[0]==y[1]:
         #   print(x,y)
        #if not y[0] == y[1]:
         #   print("different compartment same number", x,y)

### Change the SEED IDs to MNXM IDs

 #Metabolites with the following IDs are not mapped to MNXM ID: "cpd00261_c0", "cpd06227_c0", "cpd03572_c0", "cpd02446_c0","cpd02572_c0", "cpd01466_c0". They #were selected by their number of involved reactions.  

for met in model.metabolites:
    split, compartment=met.id.split("_")
    #if met.id.startswith('cpd'):
    met.id
    if met.id in ["cpd00261_c0", "cpd06227_c0", "cpd03572_c0", "cpd02446_c0","cpd02572_c0", "cpd01466_c0" ]:
        pass
    else:
        met.id=met.annotation['metanetx.chemical']+ '_' + compartment
        model.repair()
    print(met.id)
    



## Save draft with mapped metabolites

write_sbml_model(model, "/Users/lizrad/Dev/iVnat/iVnat.xml")