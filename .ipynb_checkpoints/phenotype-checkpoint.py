#!/usr/bin/env python
# coding: utf-8

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
    def __init__(self,i,description,distance,species):
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
    
def load_model(model_path):
    #Load the model (replace with actual path to model)
    model = Doc2Vec.load(model_path)
    return model

def process_traits(traitlist):
    for i, trait in traitlist:
        if str(trait)=="nan":
        trait = "N/A"
        tokens = simple_preprocess(trait)
        yield tokens
        
def make_vectors(traitlist,model,epochs=1):
    vectors=[]
    for (trait,[species,i]) in traitlist:
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
    trait_df = pd.read_csv(trait_file)
    unprocessed_traits = trait_list
    traits = list(process_traits(trait_df))
    model = load_model(model_path)
    phvecs = make_vectors(traits,model,epochs=100)
    query_vec = model.infer_vector(simple_preprocess(input_trait),epochs=100)
    dists = []
    for i,vec in enumerate(phvecs):
        #vec[0] is the vector for the pheno,
        #vec[1] is the species
        sim = spatial.distance.cosine(vec[0],query_vec)
        p = PhenoMatch(i,unprocessed_traits[i],sim,vec[1])
        dists.append(p)
    return dists


if __name__=="__main__":
    input_trait = "grain yield"
    dists=get_distances("Data/Models/dbow_wiki_2.model","Data/Phenos/All_phenotypes.csv",input_trait)





