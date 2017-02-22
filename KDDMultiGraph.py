# coding: utf-8
import networkx as nx
from networkx.classes.multigraph import MultiGraph

"""
Created on 15/mar/2013

@author: "Giulio Rossetti"
@contact:    <Giulio.Rossetti@isti.cnr.it>
"""


class KDDMultiGraph(MultiGraph):
    """
    classdocs
    """

    def __init__(self):
        """
        Constructor
        """
        self.G = nx.MultiGraph()
        self.__dimensions = []
        self.graph = self.G.graph   # dictionary for graph attributes
        self.node = self.G.node    # empty node dict (created before convert)
        self.adj = self.G.adj     # empty adjacency dict
        self.edge = self.adj

    def get_dimensions(self):
        """
        Return the set of dimensions
        """
        return self.__dimensions

    def read_edgelist(self, path, weighted=False, temporal=False):
        """
        Read the multidimensional [weighted] [time stamped] graph
        """
        edges = open(path, "r")
        for e in edges:
            e = e.replace(" ","\t")
            part = e.split('\t')
            u = int(part[0])
            v = int(part[1])
 
            d = int(part[2])
            
            if not int(part[2]) in self.__dimensions:
                self.__dimensions.append(int(part[2]))
            
            if weighted and temporal:
                w = int(part[3])
                t = int(part[4])
                self.G.add_edge(u, v, dim = d, time = t, weight = w)
               
            elif weighted and (not temporal):
                w = int(part[3])
                self.G.add_edge(u, v, dim = d, weight = w)
                
            elif (not weighted) and temporal:
                t = int(part[4])
                self.G.add_edge(u, v, dim = d, time = t)
            else:
                self.G.add_edge(u, v, dim = d)

    def get_dimension_list(self):
        """
        Return the dimensions list
        """
        return self.__dimensions
        
    def extrapolate_dimension_graph(self, d):
        """
        Build the subgraph inferred on the dimension d
        """
        G = nx.Graph()
    
        for e in self.edges_iter():
            edges = self.G.get_edge_data(e[0],e[1])
            for eid in edges:
                if edges[eid]['dim'] == d: 
                    G.add_edge(e[0], e[1])
                    if "weight" in edges[eid]:
                    #if edges[eid]['weight']:
                        G.edge[e[0]][e[1]]['weight'] = edges[eid]['weight']
        return G
        
    def extrapolate_time_dimension_graph(self, d, t):
        """
        Build the subgraph inferred on the dimension d at timpestamp t
        """
        G = nx.Graph()
        
        for e in self.G.edges_iter():
            edges = self.G.get_edge_data(e[0],e[1])
            for eid in edges:
                if edges[eid]['dim'] == d and edges[eid]['time']== t: 
                    G.add_edge(e[0], e[1])
                    if edges[eid]['weight']:
                        G.edge[e[0]][e[1]]['weight'] = edges[eid]['weight']

        return G
        
    def neighbor_set(self, u, dims):
        """
        Return the neighbors of node u on reached by edges in dims 
        NeighborSet(v,D) = {u ∈ V | ∃(u,v,d) ∈ E ∧ d ∈ D}.
        """
        neigh_set = {}
        
        neigh = self.G.neighbors(u)
        for v in neigh:
            edges = self.G.get_edge_data(u,v)
            for eid in edges:
                if edges[eid]['dim'] in dims:
                    neigh_set[v] = v
                
        return neigh_set.keys()
    
    def neighbor_xor(self, u, dims):
        """
        The function NeighborsXOR : V ×P(L)→N is defined as
        NeighborsXOR(v,D)=|{u∈V|∃d∈D:(u,v,d)∈E∧!∃d′ ∈/D:(u,v,d′)∈E}| 
        
        It computes the number of neighboring nodes connected by edges belonging only
        to dimensions in D.
        """
        neigh_set = {}
        
        neigh = self.neighbors(u)
        for v in neigh:
            edges = self.G.get_edge_data(u,v)
            for eid in edges:
                for d in edges[eid]['dim']:
                    if d in dims:
                        neigh_set[v] = v
              
                    else:
                        if v in neigh_set:
                            neigh_set.pop(v)
                            break
        
        return neigh_set.keys()
    
    def dimension_relevance(self, n, dim_set):
        """
        The function DR : V × P(L) → [0,1] is defined as
        DR(v, D) = Neighbors(v,D)/Neighbors(v)
        
        It computes the ratio between the neighbors of a node v connected by edges 
        belonging to a specific set of dimensions in D and the total number of its neighbors.
        """
        all_neigh = self.G.neighbors(n)
        neigh_set = self.neighbor_set(n, dim_set)
        
        if len(all_neigh) == 0:
            return 0
        else:
            return float(len(neigh_set))/float(len(all_neigh))
    
    def dimension_relevance_xor(self, n, dim_set):
        """
        DRXOR : V ×P(L) → [0, 1] is defined as
        DRXOR(v, D) = NeighborsXOR(v,D)/Neighbors(v)
        
        It computes the fraction of neighbors directly reachable from 
        node v following edges belonging only to dimensions D.
        """
        all_neigh = self.G.neighbors(n)
        neigh_xor = self.neighbor_xor(n, dim_set)

        if len(all_neigh) == 0:
            return 0
        else:
            return float(len(neigh_xor)) / float(len(all_neigh))
        
    def dimension_relevance_weighted(self, n, dim_set):
        """
        Weighted Dimension Relevance, is defined as
        DRW (v, D) = sum_{u∈NeighborSet(v,D))} n_{u,v,d}/n_{u,v} / Neighbors(v)
      
        where: n_{u,v,d} is the number of dimensions which label the edges between 
        two nodes u and v and that belong to D; n_{u,v} is the number of dimensions 
        which label the edges between two nodes u and v.
        """
        all_neigh = self.G.neighbors(n)
        neigh_set = self.neighbor_set(n, dim_set)
        
        if len(all_neigh) == 0:
            return 0
        else:
            val = 0
            for u in neigh_set:
                val += self.__number_of_dimension_per_edge(n, u, dim_set)
            
            return float(val)/float(len(all_neigh))
    
    def lowest_redundancy_connection(self, n, dim_set):
        """
        A node v ∈ V is said to be at Lowest Redundancy Connection if each of its neighbors 
        is reachable via only one dimension, i.e.,
        
        ∀u ∈ NeighborSet(v, L) : ∃! d ∈ L (u,v,d) ∈ E. ⊓⊔
        
        Note that if a node v is LRC we have
        Degree(v, L) = Neighbors(v, L)
        """
        if len(self.neighbor_set(n, dim_set)) == self.G.degree(n):
            return True
        else:
            return False
    
    def highest_redundancy_connection(self, n):
        """
        A node v ∈ V is called Highest Redundancy Connections if each of its neighbors 
        is reachable via all the dimensions in the network, i.e.,

        ∀u ∈ NeighborSet(v,L) : ∀d ∈ L (u,v,d) ∈ E.

        Note that if a node v is HRC we have
        Degree(v, L) = N eighbors(v, L) × |L|.
        """
        if (len(self.neighbors(n))*self.__dimensions) == self.G.degree(n):
            return True
        else:
            return False
    
    def node_dimension_connectivity(self, dim):
        """
        The function NDC : L → [0,1] is defined as
        NDC(d) = |{u∈V |∃v∈V :(u,v,d)∈E}| / |V|
        
        It computes the ratio of nodes of the network that belong to the dimension d.
        """
        count = {}

        for u in self.G.nodes():
            neigh = self.G.neighbors(u)
            for v in neigh:
                edges = self.G.get_edge_data(u, v)
                for eid in edges:
                    d = edges[eid]['dim']
                    if d == dim:
                        count[u] = u
                        break
        
        return float(len(count))/float(len(self.G.nodes()))
        
    def node_exclusive_dimension_connectivity(self, dim):
        """
        The function NEDC : L → [0,1] id defined as
        
        NEDC(d) = |{u∈V |∃v∈V :(u,v,d)∈E ∧ ∀j∈L,j!=d:(u,v,j)∈/E}| / |{u∈V |∃v∈V : (u,v,d)∈E}|

        It computes the ratio of nodes belonging only to the dimension d.
        """
        excl = {}
        
        for u in self.G.nodes():
            neigh = self.G.neighbors(u)
            for v in neigh:
                edges = self.G.get_edge_data(u, v)
                for eid in edges:
                    d = edges[eid]['dim']
                    if d == dim:
                        excl[u] = u
                    else:
                        if excl.has_key(u):
                            excl.pop(u)
                        break
                        
        return float(len(excl))/float(len(self.G.nodes()))
        
    def edge_dimension_connectivity(self, dim):
        """
        The function EDC : L → [0, 1] is defined as
        
        EDC(d) = |{(u,v,d)∈E|u,v∈V }| / |E|
        
        It computes the ratio of edges of the network labeled with the dimension d.
        """
        count = 0
        
        for el in self.G.edges_iter():
            edges = self.G.get_edge_data(el[0],el[1])
            for e in edges:
                #print edges[e]['dim']
                if int(edges[e]['dim']) == dim:
                    count += 1
                    break
        
        #print str(float(count)) + " "+ str(float(len(self.G.edges())))
        return float(count)/float(len(self.G.edges()))
    
    def edge_exclusive_dimension_connectivity(self, dim):
        """
        The function EEDC : L → [0,1] is defined as

        EEDC(d) = |{(u,v,d)∈E|u,v∈V ∧ ∀j∈L,j!=d: (u,v,j)∈/E}| / |{(u,v,d)∈E|u,v∈V }|

        It computes the ratio of edges between any pair of nodes u and v labeled with 
        the dimension d such that there are no other edges between the same two nodes 
        belonging to other dimensions j != d.
        """
        excl = {}
        
        for e in self.G.edges_iter():
            edges = self.G.get_edge_data(e[0],e[1])
            for eid in edges:
                if edges[eid]['dim'] == dim:
                    excl[e] = e
                else:
                    if excl.has_key(e):
                        excl.pop(e)
                    break
                
        return float(len(excl))/float(len(self.G.edges()))
    
    def node_d_correlation(self, dim_set):
        """
        The Node D-Correlation is the function ρnodes : P(L) → [0, 1] defined as
        
        ρnodes(D) = |forall_{d∈D} intersection V_d| / |forall_{d∈D} union V_d|
        
        where V_d denotes the set of nodes belonging to dimension d.
        It computes the ratio of nodes appearing in all the dimensions in D and the 
        total number of nodes appearing in at least one dimension in D.
        """
        total = {}
        
        for u in self.G.nodes():
            neigh = self.G.neighbors(u)
            for v in neigh:
                edges = self.G.get_edge_data(u, v)
                for eid in edges:
                    d = edges[eid]['dim']
                    if d in dim_set:
                        if total.has_key(d):
                            total[d].append(u)
                        else:
                            total[d] = [u]
        
        union = []
        intersection = []
        for d in total.keys():
            union = set(union) | set(total.get(d))
            if len(intersection) == 0:
                intersection = set(total.get(d))
            else:
                intersection = set(intersection) & set(total.get(d))

        if len(union) == 0:
            return 0          
        return float(len(intersection))/float(len(union))
    
    def pair_d_correlation(self, dim_set):
        """
        The Pair D-Correlation is the function ρpairs : P(L) → [0, 1] defined as
        
        ρpairs(D) = |forall_{d∈D} intersection P_d| / |forall_{d∈D} union P_d|
        
        where P_d denotes the set of pairs of nodes (u,v) connected in dimension d. 
        It computes the ratio of pairs of nodes connected in all the dimensions in D 
        and the total number of pairs of nodes connected in at least one dimension in D.
        """
        total = {}
        
        for e in self.G.edges_iter():
            edges = self.G.get_edge_data(e[0],e[1])
            for eid in edges:
                d = edges[eid]['dim']
                if d in dim_set:
                    if total.has_key(d):
                        total[d].append(e)
                    else:
                        total[d] = [e]
                
        union = []
        intersection = []
        for d in total.keys():            
            union = set(union) | set(total.get(d))
            if len(intersection) == 0:
                intersection = set(total.get(d))
            else:
                intersection = set(intersection) & set(total.get(d))
                    
        if len(union) == 0:
            return 0
        return float(len(intersection))/float(len(union))
     
    def __number_of_dimension_per_edge(self, u, v, dims):
        """
        Computes the mean number of dimensions per edge
        """
        total = {}
        partial_dims = {}
        
        edges = self.G.get_edge_data(u, v)
        for eid in edges:
            d = edges[eid]['dim']
            total[d] = d
            if d in dims:
                partial_dims[d]=d
        
        return float(len(partial_dims))/float(len(total))
