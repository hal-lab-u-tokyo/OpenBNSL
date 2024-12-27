#!/usr/bin/env python
from collections import deque
from itertools import permutations
import random
import networkx as nx
from tqdm.auto import trange
import numpy as np

from pgmpy import config
from pgmpy.base import DAG
from pgmpy.base import PDAG
from pgmpy.base import UndirectedGraph as Graph
import sys   # これが絶対 import
sys.path.append('/workspace')

from modules.visualize_graph import display_graph_info as show
from pgmpy.estimators import (
    AICScore,
    BDeuScore,
    BDsScore,
    BicScore,
    K2Score,
    ScoreCache,
    StructureEstimator,
    StructureScore,
)
#from pgmpy.estimators.CITests import NatoriScore
from modules.CITests_fixed import NatoriScore
from collections import defaultdict
from itertools import combinations


class RAIEstimator(StructureEstimator):
    """
    Class for heuristic RAI for DAGs, to learn
    network structure from data. `estimate` attempts to find a model with optimal score.

    Parameters
    ----------
    data: pandas DataFrame object
        dataframe object where each column represents one variable.
        (If some values in the data are missing the data cells should be set to `numpy.NaN`.
        Note that pandas converts each column containing `numpy.NaN`s to dtype `float`.)

    state_names: dict (optional)
        A dict indicating, for each variable, the discrete set of states (or values)
        that the variable can take. If unspecified, the observed values in the data set
        are taken to be the only possible states.

    use_caching: boolean
        If True, uses caching of score for faster computation.
        Note: Caching only works for scoring methods which are decomposable. Can
        give wrong results in case of custom scoring methods.

    References
    ----------
    Koller & Friedman, Probabilistic Graphical Models - Principles and Techniques, 2009
    Section 18.4.3 (page 811ff)
    """

    def __init__(self, data, use_cache=True, **kwargs):
        self.use_cache = use_cache
        super(RAIEstimator, self).__init__(data, **kwargs)
    
    @staticmethod

    def orientation(self, Gs, Gall, ci_test, data, ESS):
        """
    orient edges in a PDAG to a maximally oriented graph.
    orient rules are based on rule 1~3 from Meek,C.:Causal Inference and Causal Explanation with Background Knowledge,Proc.Confon Uncertainty in Artificial Inteligence (UAl-95),p.403-410 (195)
    """
        cls_test = ci_test(data)
                                            
        #for each X-Z-Y (X and Y is not adjecent), find V-structure
        for pair in list(combinations(Gs.nodes(), 2)):
            X, Y = pair
            for Z in (
                (set(Gs.successors(X))
                & set(Gs.predecessors(X))
                & set(Gs.successors(Y)))
                | (set(Gs.successors(X))
                & set(Gs.predecessors(Y))
                & set(Gs.successors(Y)))
            ):
                if not Gs.has_edge(X, Y):
                    if not Gs.has_edge(Y, X):
                        if not cls_test.separate(X, Y, [Z], data, ESS, boolean=True):
                            if Gs.has_edge(Z, X):
                                Gs.remove_edge(Z, X)
                                Gall.remove_edge(Z, X)
                            if Gs.has_edge(Z, Y):
                                Gs.remove_edge(Z, Y)
                                Gall.remove_edge(Z, Y)
                            print("Vstructure")
                            print(X, Z, Y)
                            break
                                    
        flag = True
        while flag: # while the number of orientated edges increases
            flag = False
            # R1: X -> Y - Z, no edge between X and Z then X -> Y -> Z
            for node_y in Gs.nodes:
                for node_x in list(Gs.predecessors(node_y)):
                    if not Gs.has_edge(node_y, node_x):
                        for node_z in list(Gs.predecessors(node_y)):
                            if node_z != node_x:
                                if Gs.has_edge(node_y, node_z):
                                    if not Gs.has_edge(node_z, node_x):
                                        if not Gs.has_edge(node_x, node_z):
                                            Gall.remove_edge(node_z, node_y)
                                            Gs.remove_edge(node_z, node_y)
                                            print("nagare")
                                            print(node_x, node_y, node_z)
                                            flag = True
                                        
            # R2: X - Y and if there is a directed path from X to Y, then X -> Y
            for node_y in Gs.nodes:
                for node_x in list(Gs.predecessors(node_y)):
                    if Gs.has_edge(node_y, node_x) :                 
                        Gs.remove_edge(node_y, node_x)
                        Gs.remove_edge(node_x, node_y)
                        if self.has_directedpath(Gs, node_x, node_y, mark=[]):
                            Gs.add_edge(node_x, node_y)
                            Gall.remove_edge(node_y, node_x)
                            print("path")
                            print(node_x, node_y)
                            flag = True
                        elif self.has_directedpath(Gs, node_y, node_x, mark=[]):
                            Gs.add_edge(node_y, node_x)
                            Gall.remove_edge(node_x, node_y)
                            print("path")
                            print(node_y, node_x)
                            flag = True
                        else:
                            Gs.add_edge(node_x, node_y)
                            Gs.add_edge(node_y, node_x)
            # R3: for each X->W<-Z X-Y-Z Y-W, orient Y->W 
            for pair in list(permutations(Gs.nodes(), 2)):
                X, Z = pair
                if not Gs.has_edge(X, Z):
                    if not Gs.has_edge(Z,X):
                        for Y in (
                            set(Gs.successors(X))
                            & set(Gs.successors(Z))
                            & set(Gs.predecessors(X))
                            & set(Gs.predecessors(Z))
                        ):
                            for W in (
                                set(Gs.successors(Y))
                            & set(Gs.predecessors(Y))
                            & set(Gs.successors(X))
                            & set(Gs.successors(Z))
                            - set(Gs.predecessors(X))
                            - set(Gs.predecessors(Z))
                            ):
                                if Gs.has_edge(Y, W) and Gs.has_edge(W, Y):
                                    Gs.remove_edge(W, Y)
                                    Gall.remove_edge(W, Y)
                                    flag = True
                                    print("R3")
        return Gs, Gall
    
    def orientation_a2(self, Gs, Gall, deletededges):
        """
    orient edges in a PDAG to a maximally oriented graph.
    orient rules are based on rule 1~3 from Meek,C.:Causal Inference and Causal Explanation with Background Knowledge,Proc.Confon Uncertainty in Artificial Inteligence (UAl-95),p.403-410 (195)
    in this stage (stage A2 in Yehezkel and Lerner(2009)), only rule 1 is applied because only X -> Y - Z shape is created in stage A1 (X-Z removed).
    """
        #Rule 1: X -> Y - Z, no edge between X and Z then X -> Y -> Z
        for pair in deletededges:
            X, Z = pair
            for Y in (
                (set(Gall.successors(X))
                & set(Gall.successors(Z))
                & set(Gall.predecessors(Z))
                - set(Gall.predecessors(X)))
            ):
                Gs.remove_edge(Z, Y)
                Gall.remove_edge(Z, Y)
        return Gs, Gall
            
    def lastorientation(self, G):                            
        flag = True
        while flag: # while the number of orientated edges increases
            flag = False
            # X -> Y - Z, no edge between X and Z then X -> Y -> Z
            for node_y in G.nodes:
                for node_x in list(G.predecessors(node_y)):
                    if not G.has_edge(node_y, node_x):
                        for node_z in list(G.predecessors(node_y)):
                            if node_z != node_x:
                                if G.has_edge(node_y, node_z):
                                    if not G.has_edge(node_z, node_x):
                                        if not G.has_edge(node_x, node_z):
                                            G.remove_edge(node_z, node_y)
                                            print("nagare")
                                            print(node_x, node_y, node_z)
                                            flag = True
                                        
            #X - Y and if there is a directed path from X to Y, then X -> Y
            for node_y in G.nodes:
                for node_x in list(G.predecessors(node_y)):
                    if G.has_edge(node_y, node_x) :                 
                        G.remove_edge(node_y, node_x)
                        G.remove_edge(node_x, node_y)
                        if self.has_directedpath(G, node_x, node_y, mark=[]):
                            G.add_edge(node_x, node_y)
                            print("path")
                            print(node_x, node_y)
                            flag = True
                        elif self.has_directedpath(G, node_y, node_x, mark=[]):
                            G.add_edge(node_y, node_x)
                            print("path")
                            print(node_y, node_x)
                            flag = True
                        else:
                            G.add_edge(node_x, node_y)
                            G.add_edge(node_y, node_x)
            #for each X->Z<-Y with Z-W or Z->W, orient edges as X->W, Y->W 
            for pair in list(permutations(G.nodes(), 2)):
                X, Y = pair
                for Z in (
                    set(G.successors(X))
                    & set(G.successors(Y))
                    - set(G.predecessors(X))
                    - set(G.predecessors(Y))
                ):
                    for W in (
                        set(G.successors(Z))
                    ):
                        if G.has_edge(W, X) and G.has_edge(X, W):
                            G.remove_edge(W, X)
                            print("zibun")
                            print(X, Y, W)
                            flag = True
                        if G.has_edge(W, Y) and G.has_edge(Y, W):
                            G.remove_edge(W, Y)
                            flag = True
        return G

    def has_directedpath(self, G, node_x, node_y, mark): #if X has a path to Y such that all edges are directed
            if node_x == node_y:
                return True
            for node_z in list(G.successors(node_x)):
                if not node_z in mark:
                    if not node_x in list(G.successors(node_z)):
                        mark.append(node_z)
                        if self.has_directedpath(G, node_z, node_y, mark):
                            return True
            return False

    def _legal_operations(
        self,
        model,
        score,
        structure_score,
        tabu_list,
        max_indegree,
        black_list,
        white_list,
        fixed_edges,
    ):
        """Generates a list of legal (= not in tabu_list) graph modifications
        for a given model, together with their score changes. Possible graph modifications:
        (1) add, (2) remove, or (3) flip a single edge. For details on scoring
        see Koller & Friedman, Probabilistic Graphical Models, Section 18.4.3.3 (page 818).
        If a number `max_indegree` is provided, only modifications that keep the number
        of parents for each node below `max_indegree` are considered. A list of
        edges can optionally be passed as `black_list` or `white_list` to exclude those
        edges or to limit the search.
        """

        tabu_list = set(tabu_list)

        # Step 1: Get all legal operations for adding edges.
        potential_new_edges = (
            set(permutations(self.variables, 2))
            - set(model.edges())
            - set([(Y, X) for (X, Y) in model.edges()])
        )

        for X, Y in potential_new_edges:
            # Check if adding (X, Y) will create a cycle.
            if not nx.has_path(model, Y, X):
                operation = ("+", (X, Y))
                if (
                    (operation not in tabu_list)
                    and ((X, Y) not in black_list)
                    and ((X, Y) in white_list)
                ):
                    old_parents = model.get_parents(Y)
                    new_parents = old_parents + [X]
                    if len(new_parents) <= max_indegree:
                        score_delta = score(Y, new_parents) - score(Y, old_parents)
                        score_delta += structure_score("+")
                        yield (operation, score_delta)

        # Step 2: Get all legal operations for removing edges
        for X, Y in model.edges():
            operation = ("-", (X, Y))
            if (operation not in tabu_list) and ((X, Y) not in fixed_edges):
                old_parents = model.get_parents(Y)
                new_parents = [var for var in old_parents if var != X]
                score_delta = score(Y, new_parents) - score(Y, old_parents)
                score_delta += structure_score("-")
                yield (operation, score_delta)

        # Step 3: Get all legal operations for flipping edges
        for X, Y in model.edges():
            # Check if flipping creates any cycles
            if not any(
                map(lambda path: len(path) > 2, nx.all_simple_paths(model, X, Y))
            ):
                operation = ("flip", (X, Y))
                if (
                    ((operation not in tabu_list) and ("flip", (Y, X)) not in tabu_list)
                    and ((X, Y) not in fixed_edges)
                    and ((Y, X) not in black_list)
                    and ((Y, X) in white_list)
                ):
                    old_X_parents = model.get_parents(X)
                    old_Y_parents = model.get_parents(Y)
                    new_X_parents = old_X_parents + [Y]
                    new_Y_parents = [var for var in old_Y_parents if var != X]
                    if len(new_X_parents) <= max_indegree:
                        score_delta = (
                            score(X, new_X_parents)
                            + score(Y, new_Y_parents)
                            - score(X, old_X_parents)
                            - score(Y, old_Y_parents)
                        )
                        score_delta += structure_score("flip")
                        yield (operation, score_delta)

    def RecursiveSearch(self, Nz, Gs, Gex, Gall, ci_test, ESS): #Gexを辞書からPDAGへ Gout削除
        # Step 1: Initial checks and setup for arguments
        # print(Nz)
        # show(Gs)
        # show(Gex)
        cls_test = ci_test(self.data)
        if all(len(Gs[node]) <= Nz for node in Gs):
            #show(Gs)
            #print(Nz)
            return Gall #Gsを戻さない
        # Step 2: Do CI tests for nodes between Gex and Gs and remove edge
        # print("Gex.nodes:")
        # print(Gex.nodes)
        deletededges = []
        for node_y in Gs.nodes:
            for node_x in Gex.nodes:
                if Gall.has_edge(node_x, node_y):
                    Z = set(Gall.predecessors(node_y))# | set(Gall.successors(node_y))
                    Z = Z - {node_x}
                    #Z = Z - {node_y} #node_xは関係ない(元は実装ミス) これで正しい？
                    if len(Z) >= Nz:
                        for Z in combinations(Z, Nz):
                            if Gall.has_edge(node_x, node_y) or Gall.has_edge(node_y, node_x):
                                if cls_test.separate(node_x, node_y, Z, self.data, ESS, boolean=True):
                                    if Gall.has_edge(node_x, node_y):
                                        Gall.remove_edge(node_x, node_y)
                                    if Gall.has_edge(node_y, node_x):
                                        Gall.remove_edge(node_y, node_x)
                                    deletededges.append((node_x, node_y))
                                    break

        if Nz >= 1:
            Gs, Gall = self.orientation_a2(Gs = Gs, Gall = Gall, deletededges = deletededges)
        #show(Gall)
        checked = set()
        for node_y in Gs.nodes:
            #neighbors = list(Gs.neighbors(node_y))
            neighbors = set(Gs.predecessors(node_y)) | set(Gs.successors(node_y))
            neighbors = neighbors - checked
            checked = checked | set(node_y)
            # print("CItest:")
            # print(node_y, neighbors)
            for node_x in neighbors:
                if Nz == 0:
                    Z = []
                    if Gall.has_edge(node_x, node_y) or Gall.has_edge(node_y, node_x):
                        if cls_test.separate(node_x, node_y, Z, self.data, ESS, boolean=True):
                            if Gall.has_edge(node_x, node_y):
                                Gall.remove_edge(node_x, node_y)
                            if Gall.has_edge(node_y, node_x):
                                Gall.remove_edge(node_y, node_x)
                            if Gs.has_edge(node_x, node_y):
                                Gs.remove_edge(node_x, node_y)
                            if Gs.has_edge(node_y, node_x):
                                Gs.remove_edge(node_y, node_x)
                else:
                    #set_Pa = set(Gall.predecessors(node_y)) | set(Gall.successors(node_y))
                    set_Pa = set(Gall.predecessors(node_y)) | set(Gall.successors(node_y)) | set(Gall.predecessors(node_x)) | set(Gall.successors(node_x))
                    #set_Pa = set(Gall.predecessors(node_y)) | set(Gall.predecessors(node_x))
                    set_Pa = set_Pa - {node_y} - {node_x} 
                    num_Pa = len(set_Pa)
                    if num_Pa >= Nz:
                        for Z in combinations(set_Pa, Nz):
                            if Gall.has_edge(node_x, node_y) or Gall.has_edge(node_y, node_x):
                                if cls_test.separate(node_x, node_y, Z, self.data, ESS, boolean=True):
                                    if Gall.has_edge(node_x, node_y):
                                        Gall.remove_edge(node_x, node_y)
                                    if Gall.has_edge(node_y, node_x):
                                        Gall.remove_edge(node_y, node_x)
                                    if Gs.has_edge(node_x, node_y):
                                        Gs.remove_edge(node_x, node_y)
                                    if Gs.has_edge(node_y, node_x):
                                        Gs.remove_edge(node_y, node_x) #Frozen graph can't be modified error
                                    break

        Gs, Gall = self.orientation(self, Gs = Gs, Gall = Gall, ci_test = ci_test, data = self.data, ESS = ESS)
        Gd, g_subs = self.order_grouping(Gs)
        #show(Gd)
        Gexd = PDAG()
        #show(Gs)
        #show(Gall)
        #show(Gs)
        for subs in g_subs:
            #Gexd = nx.compose(Gexd, subs)
            Gexd.add_nodes_from(subs.nodes)
            Gexd.add_edges_from(subs.edges)
            #show(subs)
            Gall = self.RecursiveSearch(Nz + 1, subs, Gex, Gall, ci_test, ESS)     #Recursive call for structure learning of each of the ancestor sub-structures
        #show(Gexd)
        return self.RecursiveSearch(Nz + 1, Gd, Gexd, Gall, ci_test, ESS)         #Recursive call for descendants sub-structure structure learning

    def order_grouping(self, Gs):
        sccs = list(nx.strongly_connected_components(Gs))
        scc_graph = nx.condensation(Gs)
        topological_order = list(nx.topological_sort(scc_graph))
        topological_lowest_nodes = sccs[topological_order[len(scc_graph.nodes) - 1]]
        gc = self.extract_subgraph(
            topological_lowest_nodes,
            Gs
        )
        SubStructures = []
        for sccs_idx in topological_order[0 : len(scc_graph.nodes) - 1]:
            sub = self.extract_subgraph(
                sccs[sccs_idx],
                Gs
            )
            # sub = Gs.subgraph(                #うまくいかない？速いがSHD高
            #     sccs[sccs_idx]
            # ).copy()
            SubStructures.append(sub)
        return gc, SubStructures

    def extract_subgraph(self, nodes, G):
        sub = PDAG()
        sub.add_nodes_from(nodes)
        edges_to_add = [(node, sub_node) for node in nodes for sub_node in nodes if node != sub_node and G.has_edge(node, sub_node)]
        sub.add_edges_from(edges_to_add)
        return sub

    def estimate(        
        self,
        scoring_method="natoriscore",
        Gs=None,
        max_indegree=None,
        show_progress=True,
        ESS = 1
    ):

        # Step 1: Initial checks and setup for arguments
        # Step 1.1: Check scoring_method
        supported_methods = {
            "k2score": K2Score,
            "bdeuscore": BDeuScore,
            "bdsscore": BDsScore,
            "bicscore": BicScore,
            "aicscore": AICScore,
            "natoriscore": NatoriScore,
        }
        if (
            (
                isinstance(scoring_method, str)
                and (scoring_method.lower() not in supported_methods)
            )
        ) and (not isinstance(scoring_method, StructureScore)):
            raise ValueError(
                "scoring_method should either be one of k2score, bdeuscore, bicscore, bdsscore, aicscore, or an instance of StructureScore"
            )

        if isinstance(scoring_method, str):
            ci_test = supported_methods[scoring_method.lower()]
        else:
            ci_test = scoring_method

        # Step 1.2: Check the start_dag
        Gs = PDAG()
        Gs.add_nodes_from(self.variables)
        Gs.add_edges_from(
            [
                (X, Y)
                for X, Y in permutations(self.variables, 2)     #combination から変更
                if X != Y
            ]
        )
        Gall = PDAG()
        Gall.add_nodes_from(self.variables)
        Gall.add_edges_from(
            [
                (X, Y)
                for X, Y in permutations(self.variables, 2)     #combination から変更
                if X != Y
            ]
        )
        Gex = PDAG()
        Nz = 0
        #show(Gs)
        # Step 2: Define the structure_score function
        best_model = self.RecursiveSearch(
            Nz=Nz,
            Gs=Gs,
            Gex=Gex,
            Gall=Gall,
            ci_test=ci_test,
            ESS = ESS
        )
        #best_model = self.lastorientation(best_model)

        return best_model


class RAIEstimator_transitivity(StructureEstimator):
    """
    Class for heuristic RAI for DAGs using transitivity to reduce the number of CI tests.
     `estimate` attempts to find a model with optimal score.

    Parameters
    ----------
    data: pandas DataFrame object
        dataframe object where each column represents one variable.
        (If some values in the data are missing the data cells should be set to `numpy.NaN`.
        Note that pandas converts each column containing `numpy.NaN`s to dtype `float`.)

    state_names: dict (optional)
        A dict indicating, for each variable, the discrete set of states (or values)
        that the variable can take. If unspecified, the observed values in the data set
        are taken to be the only possible states.

    use_caching: boolean
        If True, uses caching of score for faster computation.
        Note: Caching only works for scoring methods which are decomposable. Can
        give wrong results in case of custom scoring methods.

    References
    ----------
    Koller & Friedman, Probabilistic Graphical Models - Principles and Techniques, 2009
    Section 18.4.3 (page 811ff)
    本田和雅, 名取和樹, 菅原聖太, 磯崎隆司, 植野真臣：推移性を利用した大規模ベイジアンネットワーク構造学習,電子情報通信学会論文誌 D, vol.J102-D, No.12, pp.796-811 (2019)
    """

    def __init__(self, data, use_cache=True, **kwargs):
        self.use_cache = use_cache
        super(RAIEstimator_transitivity, self).__init__(data, **kwargs)
    
    @staticmethod

    def transive_cut(self, Gs, Gall, ci_test, data, ESS, X, Y, Z, deletededges):
        cls_test = ci_test(data)
        A = ((set(Gall.successors(X))
             | set(Gall.predecessors(X)))
             & (set(Gall.successors(Y))
             | set(Gall.predecessors(Y)))
             - ((set(Gall.successors(X))
                & set(Gall.successors(Y)))
                | set(Z))
             )
        for a in A:
            if not ((X, a) in deletededges or (a, X) in deletededges):
                if cls_test.separate(X, a, Z, data, ESS, boolean=True):
                    if Gs.has_node(X) and Gs.has_node(a):
                        if Gs.has_edge(X, a):
                            Gs.remove_edge(X, a)
                            Gall.remove_edge(X, a)
                        if Gs.has_edge(a, X):
                            Gs.remove_edge(a, X)
                            Gall.remove_edge(a, X)
                    else:
                        if Gall.has_edge(X, a):
                            Gall.remove_edge(X, a)
                        if Gall.has_edge(a, X):
                            Gall.remove_edge(a, X)
                    deletededges.append((X, a))
                    return Gs, Gall, deletededges
            if not ((Y, a) in deletededges or (a, Y) in deletededges):
                if cls_test.separate(Y, a, Z, data, ESS, boolean=True):
                    if Gs.has_node(Y) and Gs.has_node(a):
                        if Gs.has_edge(Y, a):
                            Gs.remove_edge(Y, a)
                            Gall.remove_edge(Y, a)
                        if Gs.has_edge(a, Y):
                            Gs.remove_edge(a, Y)
                            Gall.remove_edge(a,Y)
                    else:
                        if Gall.has_edge(Y, a):
                            Gall.remove_edge(Y, a)
                        if Gall.has_edge(a, Y):
                            Gall.remove_edge(a, Y)
                    deletededges.append((Y, a))
                    return Gs, Gall, deletededges                
        return Gs, Gall, deletededges

    def orientation(self, Gs, Gall, ci_test, data, ESS):
        """
    orient edges in a PDAG to a maximally oriented graph.
    orient rules are based on rule 1~3 from Meek,C.:Causal Inference and Causal Explanation with Background Knowledge,Proc.Confon Uncertainty in Artificial Inteligence (UAl-95),p.403-410 (195)
    """
        cls_test = ci_test(data)
                                            
        #for each X-Z-Y (X and Y is not adjecent), find V-structure
        for pair in list(combinations(Gs.nodes(), 2)):
            X, Y = pair
            for Z in (
                (set(Gs.successors(X))
                & set(Gs.predecessors(X))
                & set(Gs.successors(Y)))
                | (set(Gs.successors(X))
                & set(Gs.predecessors(Y))
                & set(Gs.successors(Y)))
            ):
                if not Gs.has_edge(X, Y):
                    if not Gs.has_edge(Y, X):
                        if not cls_test.separate(X, Y, [Z], data, ESS, boolean=True):
                            if Gs.has_edge(Z, X):
                                Gs.remove_edge(Z, X)
                                Gall.remove_edge(Z, X)
                            if Gs.has_edge(Z, Y):
                                Gs.remove_edge(Z, Y)
                                Gall.remove_edge(Z, Y)
                            break
                                    
        flag = True
        while flag: # while the number of orientated edges increases
            flag = False
            # R1: X -> Y - Z, no edge between X and Z then X -> Y -> Z
            for node_y in Gs.nodes:
                for node_x in list(Gs.predecessors(node_y)):
                    if not Gs.has_edge(node_y, node_x):
                        for node_z in list(Gs.predecessors(node_y)):
                            if node_z != node_x:
                                if Gs.has_edge(node_y, node_z):
                                    if not Gs.has_edge(node_z, node_x):
                                        if not Gs.has_edge(node_x, node_z):
                                            Gall.remove_edge(node_z, node_y)
                                            Gs.remove_edge(node_z, node_y)
                                            flag = True
                                        
            # R2: X - Y and if there is a directed path from X to Y, then X -> Y
            for node_y in Gs.nodes:
                for node_x in list(Gs.predecessors(node_y)):
                    if Gs.has_edge(node_y, node_x) :                 
                        Gs.remove_edge(node_y, node_x)
                        Gs.remove_edge(node_x, node_y)
                        if self.has_directedpath(Gs, node_x, node_y, mark=[]):
                            Gs.add_edge(node_x, node_y)
                            Gall.remove_edge(node_y, node_x)
                            flag = True
                        elif self.has_directedpath(Gs, node_y, node_x, mark=[]):
                            Gs.add_edge(node_y, node_x)
                            Gall.remove_edge(node_x, node_y)
                            flag = True
                        else:
                            Gs.add_edge(node_x, node_y)
                            Gs.add_edge(node_y, node_x)
            # R3: for each X->W<-Z X-Y-Z Y-W, orient Y->W 
            for pair in list(permutations(Gs.nodes(), 2)):
                X, Z = pair
                if not Gs.has_edge(X, Z):
                    if not Gs.has_edge(Z,X):
                        for Y in (
                            set(Gs.successors(X))
                            & set(Gs.successors(Z))
                            & set(Gs.predecessors(X))
                            & set(Gs.predecessors(Z))
                        ):
                            for W in (
                                set(Gs.successors(Y))
                            & set(Gs.predecessors(Y))
                            & set(Gs.successors(X))
                            & set(Gs.successors(Z))
                            - set(Gs.predecessors(X))
                            - set(Gs.predecessors(Z))
                            ):
                                if Gs.has_edge(Y, W) and Gs.has_edge(W, Y):
                                    Gs.remove_edge(W, Y)
                                    Gall.remove_edge(W, Y)
                                    flag = True
        return Gs, Gall
    
    def orientation_a2(self, Gs, Gall, deletededges):
        """
    orient edges in a PDAG to a maximally oriented graph.
    orient rules are based on rule 1~3 from Meek,C.:Causal Inference and Causal Explanation with Background Knowledge,Proc.Confon Uncertainty in Artificial Inteligence (UAl-95),p.403-410 (195)
    in this stage (stage A2 in Yehezkel and Lerner(2009)), only rule 1 is applied because only X -> Y - Z shape is created in stage A1 (X-Z removed).
    """
        #Rule 1: X -> Y - Z, no edge between X and Z then X -> Y -> Z
        for pair in deletededges:
            X, Z = pair
            for Y in (
                (set(Gall.successors(X))
                & set(Gall.successors(Z))
                & set(Gall.predecessors(Z))
                - set(Gall.predecessors(X)))
            ):
                Gs.remove_edge(Z, Y)
                Gall.remove_edge(Z, Y)
        return Gs, Gall
    
    def lastorientation(self, G):                            
        flag = True
        while flag: # while the number of orientated edges increases
            flag = False
            # X -> Y - Z, no edge between X and Z then X -> Y -> Z
            for node_y in G.nodes:
                for node_x in list(G.predecessors(node_y)):
                    if not G.has_edge(node_y, node_x):
                        for node_z in list(G.predecessors(node_y)):
                            if node_z != node_x:
                                if G.has_edge(node_y, node_z):
                                    if not G.has_edge(node_z, node_x):
                                        if not G.has_edge(node_x, node_z):
                                            G.remove_edge(node_z, node_y)
                                            flag = True
                                        
            #X - Y and if there is a directed path from X to Y, then X -> Y
            for node_y in G.nodes:
                for node_x in list(G.predecessors(node_y)):
                    if G.has_edge(node_y, node_x) :                 
                        G.remove_edge(node_y, node_x)
                        G.remove_edge(node_x, node_y)
                        if self.has_directedpath(G, node_x, node_y, mark=[]):
                            G.add_edge(node_x, node_y)
                            flag = True
                        elif self.has_directedpath(G, node_y, node_x, mark=[]):
                            G.add_edge(node_y, node_x)
                            flag = True
                        else:
                            G.add_edge(node_x, node_y)
                            G.add_edge(node_y, node_x)
            #for each X->Z<-Y with Z-W or Z->W, orient edges as X->W, Y->W 
            for pair in list(permutations(G.nodes(), 2)):
                X, Y = pair
                for Z in (
                    set(G.successors(X))
                    & set(G.successors(Y))
                    - set(G.predecessors(X))
                    - set(G.predecessors(Y))
                ):
                    for W in (
                        set(G.successors(Z))
                    ):
                        if G.has_edge(W, X) and G.has_edge(X, W):
                            G.remove_edge(W, X)
                            flag = True
                        if G.has_edge(W, Y) and G.has_edge(Y, W):
                            G.remove_edge(W, Y)
                            flag = True
        return G

    def has_directedpath(self, G, node_x, node_y, mark): #if X has a path to Y such that all edges are directed
            if node_x == node_y:
                return True
            for node_z in list(G.successors(node_x)):
                if not node_z in mark:
                    if not node_x in list(G.successors(node_z)):
                        mark.append(node_z)
                        if self.has_directedpath(G, node_z, node_y, mark):
                            return True
            return False

    def _legal_operations(
        self,
        model,
        score,
        structure_score,
        tabu_list,
        max_indegree,
        black_list,
        white_list,
        fixed_edges,
    ):
        """Generates a list of legal (= not in tabu_list) graph modifications
        for a given model, together with their score changes. Possible graph modifications:
        (1) add, (2) remove, or (3) flip a single edge. For details on scoring
        see Koller & Friedman, Probabilistic Graphical Models, Section 18.4.3.3 (page 818).
        If a number `max_indegree` is provided, only modifications that keep the number
        of parents for each node below `max_indegree` are considered. A list of
        edges can optionally be passed as `black_list` or `white_list` to exclude those
        edges or to limit the search.
        """

        tabu_list = set(tabu_list)

        # Step 1: Get all legal operations for adding edges.
        potential_new_edges = (
            set(permutations(self.variables, 2))
            - set(model.edges())
            - set([(Y, X) for (X, Y) in model.edges()])
        )

        for X, Y in potential_new_edges:
            # Check if adding (X, Y) will create a cycle.
            if not nx.has_path(model, Y, X):
                operation = ("+", (X, Y))
                if (
                    (operation not in tabu_list)
                    and ((X, Y) not in black_list)
                    and ((X, Y) in white_list)
                ):
                    old_parents = model.get_parents(Y)
                    new_parents = old_parents + [X]
                    if len(new_parents) <= max_indegree:
                        score_delta = score(Y, new_parents) - score(Y, old_parents)
                        score_delta += structure_score("+")
                        yield (operation, score_delta)

        # Step 2: Get all legal operations for removing edges
        for X, Y in model.edges():
            operation = ("-", (X, Y))
            if (operation not in tabu_list) and ((X, Y) not in fixed_edges):
                old_parents = model.get_parents(Y)
                new_parents = [var for var in old_parents if var != X]
                score_delta = score(Y, new_parents) - score(Y, old_parents)
                score_delta += structure_score("-")
                yield (operation, score_delta)

        # Step 3: Get all legal operations for flipping edges
        for X, Y in model.edges():
            # Check if flipping creates any cycles
            if not any(
                map(lambda path: len(path) > 2, nx.all_simple_paths(model, X, Y))
            ):
                operation = ("flip", (X, Y))
                if (
                    ((operation not in tabu_list) and ("flip", (Y, X)) not in tabu_list)
                    and ((X, Y) not in fixed_edges)
                    and ((Y, X) not in black_list)
                    and ((Y, X) in white_list)
                ):
                    old_X_parents = model.get_parents(X)
                    old_Y_parents = model.get_parents(Y)
                    new_X_parents = old_X_parents + [Y]
                    new_Y_parents = [var for var in old_Y_parents if var != X]
                    if len(new_X_parents) <= max_indegree:
                        score_delta = (
                            score(X, new_X_parents)
                            + score(Y, new_Y_parents)
                            - score(X, old_X_parents)
                            - score(Y, old_Y_parents)
                        )
                        score_delta += structure_score("flip")
                        yield (operation, score_delta)

    def RecursiveSearch(self, Nz, Gs, Gex, Gall, ci_test, ESS):
        # Step 1: Initial checks and setup for arguments
        cls_test = ci_test(self.data)
        if all(len(Gs[node]) <= Nz for node in Gs):
            return Gall
        # Step 2: Do CI tests for nodes between Gex and Gs and remove edge
        deletededges = []
        for node_y in Gs.nodes:
            for node_x in Gex.nodes:
                if Gall.has_edge(node_x, node_y):
                    if not ((node_x, node_y) in deletededges or (node_y, node_x) in deletededges):
                        Z = set(Gall.predecessors(node_y))
                        Z = Z - {node_x}
                        if len(Z) >= Nz:
                            if Gall.has_edge(node_x, node_y) or Gall.has_edge(node_y, node_x):
                                for Z in combinations(Z, Nz):
                                    if cls_test.separate(node_x, node_y, Z, self.data, ESS, boolean=True):
                                        if Gall.has_edge(node_x, node_y):
                                            Gall.remove_edge(node_x, node_y)
                                        if Gall.has_edge(node_y, node_x):
                                            Gall.remove_edge(node_y, node_x)
                                        Gs, Gall, deletededges = self.transive_cut(self, Gs, Gall, ci_test, self.data, ESS, node_x, node_y, Z, deletededges)
                                        deletededges.append((node_x, node_y))
                                        break
        if Nz >= 1:
            Gs, Gall = self.orientation_a2(Gs = Gs, Gall = Gall, deletededges = deletededges)
        checked = set()
        deletededges = []
        for node_y in Gs.nodes:
            neighbors = set(Gs.predecessors(node_y)) | set(Gs.successors(node_y))
            neighbors = neighbors - checked
            checked = checked | set(node_y)
            for node_x in neighbors:
                if not ((node_x, node_y) in deletededges or (node_y, node_x) in deletededges):
                    if Nz == 0:
                        Z = []
                        if Gall.has_edge(node_x, node_y) or Gall.has_edge(node_y, node_x):
                            if cls_test.separate(node_x, node_y, Z, self.data, ESS, boolean=True):
                                if Gall.has_edge(node_x, node_y):
                                    Gall.remove_edge(node_x, node_y)
                                if Gall.has_edge(node_y, node_x):
                                    Gall.remove_edge(node_y, node_x)
                                if Gs.has_edge(node_x, node_y):
                                    Gs.remove_edge(node_x, node_y)
                                if Gs.has_edge(node_y, node_x):
                                    Gs.remove_edge(node_y, node_x)
                                deletededges.append((node_x, node_y))
                                Gs, Gall, deletededges = self.transive_cut(self, Gs, Gall, ci_test, self.data, ESS, node_x, node_y, Z, deletededges)
                    else:
                        set_Pa = set(Gall.predecessors(node_y)) | set(Gall.successors(node_y)) | set(Gall.predecessors(node_x)) | set(Gall.successors(node_x))
                        set_Pa = set_Pa - {node_y} - {node_x} 
                        num_Pa = len(set_Pa)
                        if num_Pa >= Nz:
                            for Z in combinations(set_Pa, Nz):
                                if Gall.has_edge(node_x, node_y) or Gall.has_edge(node_y, node_x):
                                    if cls_test.separate(node_x, node_y, Z, self.data, ESS, boolean=True):
                                        if Gall.has_edge(node_x, node_y):
                                            Gall.remove_edge(node_x, node_y)
                                        if Gall.has_edge(node_y, node_x):
                                            Gall.remove_edge(node_y, node_x)
                                        if Gs.has_edge(node_x, node_y):
                                            Gs.remove_edge(node_x, node_y)
                                        if Gs.has_edge(node_y, node_x):
                                            Gs.remove_edge(node_y, node_x) 
                                        deletededges.append((node_x, node_y))
                                        Gs, Gall, deletededges = self.transive_cut(self, Gs, Gall, ci_test, self.data, ESS, node_x, node_y, Z, deletededges)
                                        break

        Gs, Gall = self.orientation(Gs, Gall = Gall, ci_test = ci_test, data = self.data, ESS = ESS)
        Gd, g_subs = self.order_grouping(Gs)
        Gexd = PDAG()
        for subs in g_subs:
            Gexd.add_nodes_from(subs.nodes)
            Gexd.add_edges_from(subs.edges)
            Gall = self.RecursiveSearch(Nz + 1, subs, Gex, Gall, ci_test, ESS)     # Recursive call for structure learning of each of the ancestor sub-structures
        return self.RecursiveSearch(Nz + 1, Gd, Gexd, Gall, ci_test, ESS)        # Recursive call for descendants sub-structure structure learning

    def order_grouping(self, Gs):
        sccs = list(nx.strongly_connected_components(Gs))
        scc_graph = nx.condensation(Gs)
        topological_order = list(nx.topological_sort(scc_graph))
        topological_lowest_nodes = sccs[topological_order[len(scc_graph.nodes) - 1]]
        gc = self.extract_subgraph(
            topological_lowest_nodes,
            Gs
        )
        SubStructures = []
        for sccs_idx in topological_order[0 : len(scc_graph.nodes) - 1]:
            sub = self.extract_subgraph(
                sccs[sccs_idx],
                Gs
            )
            SubStructures.append(sub)
        return gc, SubStructures

    def extract_subgraph(self, nodes, G):
        sub = PDAG()
        sub.add_nodes_from(nodes)
        edges_to_add = [(node, sub_node) for node in nodes for sub_node in nodes if node != sub_node and G.has_edge(node, sub_node)]
        sub.add_edges_from(edges_to_add)
        return sub

    def estimate(               #全てをPDAGで定義
        self,
        scoring_method="natoriscore",
        Gs=None,
        max_indegree=None,
        show_progress=True,
        ESS = 1
    ):

        # Step 1: Initial checks and setup for arguments
        # Step 1.1: Check scoring_method
        supported_methods = {
            "k2score": K2Score,
            "bdeuscore": BDeuScore,
            "bdsscore": BDsScore,
            "bicscore": BicScore,
            "aicscore": AICScore,
            "natoriscore": NatoriScore,
        }
        if (
            (
                isinstance(scoring_method, str)
                and (scoring_method.lower() not in supported_methods)
            )
        ) and (not isinstance(scoring_method, StructureScore)):
            raise ValueError(
                "scoring_method should either be one of k2score, bdeuscore, bicscore, bdsscore, aicscore, or an instance of StructureScore"
            )

        if isinstance(scoring_method, str):
            ci_test = supported_methods[scoring_method.lower()]
        else:
            ci_test = scoring_method

        # Step 1.2: Check the start_dag
        Gs = PDAG()
        Gs.add_nodes_from(self.variables)
        Gs.add_edges_from(
            [
                (X, Y)
                for X, Y in permutations(self.variables, 2)     #combination から変更
                if X != Y
            ]
        )
        Gall = PDAG()
        Gall.add_nodes_from(self.variables)
        Gall.add_edges_from(
            [
                (X, Y)
                for X, Y in permutations(self.variables, 2)     #combination から変更
                if X != Y
            ]
        )
        Gex = PDAG()
        Nz = 0
        #show(Gs)
        # Step 2: Define the structure_score function
        best_model = self.RecursiveSearch(
            Nz=Nz,
            Gs=Gs,
            Gex=Gex,
            Gall=Gall,
            ci_test=ci_test,
            ESS = ESS
        )
        #best_model = self.lastorientation(best_model)

        return best_model