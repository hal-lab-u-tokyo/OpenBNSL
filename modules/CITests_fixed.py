import numpy as np

from scipy.special import gammaln
from pgmpy.estimators import BDeuScore
from pgmpy.models import BayesianNetwork


class NatoriScore():
    def __init__(self, data, equivalent_sample_size=100, **kwargs):
        self.equivalent_sample_size = equivalent_sample_size


    def calc_g1(self, x_table, y_table, alpha):
        x_nijk = np.sum(x_table, axis=1)
        y_nijk = np.sum(y_table, axis=1)
        # n_ij = np.concatenate([x_nij, y_nij])
        # n_ijk = np.concatenate([x_table, y_table])
        n_ijk = np.hstack((x_table, y_table))
        x_nij = np.sum(x_nijk, axis=-1)
        y_nij = np.sum(y_nijk, axis=-1)
        r1 = x_table.shape[-1]
        r2 = y_table.shape[-1]
        g1_term1 = np.sum(gammaln(r1 * alpha) - gammaln(r1 * alpha + x_nij)) \
                    + np.sum(gammaln(r2 * alpha) - gammaln(r2 * alpha + y_nij))
        # g1_term1 = np.sum(gammaln(r1 * alpha) - gammaln(r1 * alpha + n_ij[r1:]))
        g1_term2 = np.sum(gammaln(alpha + n_ijk) - gammaln(alpha))
        return g1_term1 + g1_term2


    def calc_g3(self, n_xyz, r1, r2, alpha):
        n_j = np.sum(n_xyz, axis=1)
        g3_term1 = np.sum(gammaln(r1 * r2 * alpha) - gammaln(r1 * r2 * alpha + n_j))
        g3_term2 = np.sum(gammaln(alpha + n_xyz) - gammaln(alpha))
        return g3_term1 + g3_term2

    def bayes_factor(self, X, Y, Z, data, ESS):
        if Z:
            g1 = BayesianNetwork()#従属モデル
            g1.add_node(X)
            g1.add_node(Y)
            g1.add_nodes_from(Z)
            for node_z in Z:
                g1.add_edge(node_z, X)
                g1.add_edge(node_z, Y)
            g1.add_edge(X, Y)
            g3 = BayesianNetwork()#独立モデル
            g3.add_node(X)
            g3.add_node(Y)
            g3.add_nodes_from(Z)
            for node_z in Z:
                g3.add_edge(node_z, X)
                g3.add_edge(node_z, Y)
        else:
            g1 = BayesianNetwork()#従属モデル
            g1.add_node(X)
            g1.add_node(Y)
            g1.add_edge(X, Y)
            g3 = BayesianNetwork()#独立モデル
            g3.add_node(X)
            g3.add_node(Y)
        g1_score = BDeuScore(data, equivalent_sample_size=ESS).score(g1)
        #g1_score = BicScore(data).score(g1)
        # print("dep")
        print(g1_score)
        g3_score = BDeuScore(data, equivalent_sample_size=ESS).score(g3)
        #g3_score = BicScore(data).score(g3)
        # print("in")
        print(g3_score)
        if g3_score > g1_score:
            print("independent")
            return True #independent
        else:
            print("dependent") 
            return False
            

    def separate(self, X, Y, Z, data, ESS, boolean=True):
        if hasattr(Z, "__iter__"):
            Z = list(Z)
        else:
            raise (f"Z must be an iterable. Got object type: {type(Z)}")

        if (X in Z) or (Y in Z):
            raise ValueError(
                f"The variables X or Y can't be in Z. Found {X if X in Z else Y} in Z."
            )

        # Step 2: Do a simple contingency test if there are no conditional variables.
        independent = self.bayes_factor(
            X, Y, Z, data, ESS
        )
        if independent:
            #self.separating_sets[frozenset([X,Y])] = set(Z)
            return True
        else: return False
