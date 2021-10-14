#!/usr/bin/env python
# coding: utf-8

import sanitize
from gensim.models.doc2vec import Doc2Vec,TaggedDocument
from gensim.utils import simple_preprocess
import numpy as np
import pandas as pd
import smart_open
from scipy import spatial
from timeit import default_timer as timer
from collections import defaultdict
from sys import argv

###IMPORTANT###
#The Phenotype table must have the columns "species" and "trait". Everything else is optional.

class PhenoMatch(object):
    def __init__(self,i,description,distance,species="NA"):
        self.i = i
        self.description = description
        self.distance=distance
        self.species=species
    def __str__(self):
        pdict = {'index':self.i,
                'description':self.description,
                'distance':self.distance,
                'species':self.species}
        return(str(pdict))

def pheno_input(phenofile):
    with open(phenofile) as f:
        phenotype = f.read()
    return phenotype

def load_model(model_path):
    #Load the model (replace with actual path to model)
    model = Doc2Vec.load(model_path)
    return model

def process_traits(traitlist):
    for i, trait in enumerate(traitlist):
        if str(trait)=="nan":
            trait = "N/A"
        tokens = simple_preprocess(trait)
        yield tokens
        
def make_vectors(traitlist,model,epochs=1,species="NA"):
    vectors=[]
    for trait in traitlist:
        v = model.infer_vector(trait)
        vectors.append(([model.infer_vector(trait,epochs=epochs)],species))
    return vectors

def find_sim_phenos(query,targets,words,topn=3):
    sims = {}
    for i, trait in enumerate(targets):
        sim = spatial.distance.cosine(query,trait)
        sims[i] = sim
    sims = [(k,v) for k,v in sorted(sims.items(),key=lambda item: item[1])]
    ids = [t[0] for t in sims[:topn]]
    print(ids)
    return [words[i] for i in ids]

def get_distances(model_path,trait_list,input_trait):
    unprocessed_traits = trait_list
    traits = list(process_traits(unprocessed_traits))
    model = load_model(model_path)
    phvecs = make_vectors(traits,model,epochs=100)
    query_vec = model.infer_vector(simple_preprocess(input_trait),epochs=100)
    dists = []
    for i,vec in enumerate(phvecs):
        #vec[0] is the vector for the pheno,
        #vec[1] is the species
        sim = 1 - spatial.distance.cosine(vec[0],query_vec)
        p = PhenoMatch(i,unprocessed_traits[i],sim)
        dists.append(p)
    return dists

def get_pheno_results(in_pheno,pheno_model_file,trait_list,pheno_sim_min):
    #sanitize pheno input:
    #pheno_text = pheno_input(in_pheno)
    pheno_text = in_pheno
    pheno_text = sanitize.sanitize_text(pheno_text)

    pheno_dists = get_distances(pheno_model_file,trait_list,pheno_text)
    pheno_dist_d =[(pheno_text,p.description,p.distance,'text') for p in pheno_dists]
    pheno_table = pd.DataFrame(pheno_dist_d,columns=['input pheno','db pheno','similarity','connection'])
    pheno_table = pheno_table.assign(type1='input pheno',type2='db pheno')
    
    #only want phenos above the cutoff similarity
    pheno_table = pheno_table.loc[pheno_table['similarity'] > pheno_sim_min]
    return pheno_table


if __name__=="__main__":
    input_trait = "grain yield"
    dists=get_distances("Data/Models/dbow_wiki_2.model","Data/Phenos/All_phenotypes.csv",input_trait)
