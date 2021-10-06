import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

class ScoreGraph(nx.DiGraph):
    
    #The init must contain super() for the Graph functions to work
    def __init__(self):
        super().__init__()
        #self.query_gene = query_gene
        #self.query_pheno = query_pheno
        #self.add_node(self.query_gene,source='input gene')
        #self.add_node(self.query_pheno,source='input pheno')
        self.color_dict={
            "homology" : 'green',
            "qtl" : 'blue',
            "input qtl" : 'purple',
            "text" : 'red',
            "coexpression": 'orange'
        }
        self.conn_dict={
            ('input gene','db gene'):'homology',
            ('db gene','db pheno'):'qtl',
            ('input gene','input pheno'):'input qtl',
            ('db pheno','input pheno'):'text',
            ('db gene','db gene'): 'coexpression',
        }
    
    #Must have a table (results_table) with the following header:
    #['item1','item2','weight','type1','type2','connection']
    def add_all_nodes_and_edges(self,results_table):
        result_records = results_table.to_dict(orient='records')
        def add_row_to_nodes(n1,n2,s1,s2):
            if n1 not in self:
                self.add_node(n1,source=s1)
            if n2 not in self:
                self.add_node(n2,source=s2)
        def add_row_to_edges(n1,n2,s1,s2,weight):
            if (s1,s2) not in self.conn_dict:
                try:
                    (s1,s2) = (s1,s2)[::-1]
                    assert (s1,s2) in self.conn_dict
                except AssertionError:
                    print("Error: connection type is invalid: {}".format((s1,s2)))
            conn = self.conn_dict[(s1,s2)]
            ecolor = self.color_dict[conn]
            self.add_edge(n1,n2,connection=conn,weight=weight,color=ecolor,lweight=weight*5)
        for record in result_records:
            n1 = record['item1']
            n2 = record['item2']
            s1 = record['source1']
            s2 = record['source2']
            weight = record['weight']
            add_row_to_nodes(n1,n2,s1,s2)
            add_row_to_edges(n1,n2,s1,s2,weight)
        
    #this one should only be called when all adges are added
    
    def candidate_scores(self,query_genes,query_pheno):
        cand_dict = {}
        for query_gene in query_genes:
            #the "simple_paths" method finds paths from query gene->
            #query pheno - *almost* a cycle
            gpaths = sorted(list(nx.all_simple_paths(self,query_gene,
                                                    query_pheno)))
            #print(gpaths)
            scores = []
            for path in gpaths:
                total_weight = 1
                for i in range(len(path)-1):
                    source, target = path[i], path[i+1]
                    edge = self[source][target]
                    weight = edge['weight']
                    total_weight *= weight
                scores.append(total_weight)
            #Address issue of frequent genes
            gscore = sum(scores)
            cand_dict[query_gene] = gscore
        for gene, score in cand_dict.items():
            self.add_edge(gene,query_pheno,
            connection='input qtl',
            weight=score,
            color=self.color_dict[conn],
            lweight=weight*5)
        return cand_dict

    
    
