import networkx as nx
import matplotlib.pyplot as plt

class ScoreGraph(nx.DiGraph):
    
    #The init must contain super() for the Graph functions to work
    def __init__(self,query_gene,query_pheno):
        super().__init__()
        self.query_gene = query_gene
        self.query_pheno = query_pheno
        self.add_node(self.query_gene,source='input gene')
        self.add_node(self.query_pheno,source='input pheno')
        self.color_dict={
            "homology" : 'green',
            "qtl" : 'blue',
            "in_qtl" : 'purple',
            "text" : 'red',
            "coexpression": 'orange'
        }
        self.conn_dict={
            ('input gene','db gene'):'homology',
            ('db gene','db pheno'):'qtl',
            ('input gene','input pheno'):'in_qtl',
            ('db pheno','input pheno'):'text',
            ('db gene','db gene'): 'coexpression',
        }
    
    #Must create a dict with {gene/pheno: what_kind_is_it} pairs
    def add_all_nodes(self,node_dict):
        for k,v in node_dict.items():
            if v == 'input pheno':
                k = k + " query"
            self.add_node(k,source=v)
        
    
    #this is envisioned to be called in a loop through a dataframe
    #edges MUST be added in the correct order: n1 --> n2
    def add_score_edge(self,n1,n2,weight):
        types = (self.nodes[n1]['source'],
                 self.nodes[n2]['source'])
        if types not in self.conn_dict:
            try:
                types = types[::-1]
                assert types in self.conn_dict
            except AssertionError:
                print("Error: connection type is invalid.")
        conn = self.conn_dict[types]
        ecolor = self.color_dict[conn]
        self.add_edge(n1,n2,connection=conn,weight=weight,ecolor=ecolor,lweight=weight*5)
        
    #this one should only be called when all adges are added
    def candidate_score(self):
        gpaths = sorted(list(nx.all_simple_paths(self,self.query_gene,
                                                 self.query_pheno + " query")))
        print(gpaths)
        scores = []
        for path in gpaths:
            total_weight = 1
            for i in range(len(path)-1):
                source, target = path[i], path[i+1]
                edge = self[source][target]
                weight = edge['weight']
                total_weight *= weight
            scores.append(total_weight)
        return sum(scores)

    
    