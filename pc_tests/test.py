import os, sys, signal

import multiprocessing
from multiprocessing.sharedctypes import Value, Array

import numpy

from sklearn.covariance import graph_lasso
from sklearn import linear_model
from sklearn.cross_validation import KFold
from sklearn import preprocessing

import networkx as nx

NTHREADS = 12

class ThreadSafeFile( file ):
    def __init__( *args ):
        file.__init__( *args )
        args[0].lock = multiprocessing.Lock()

    def write( self, string ):
        with self.lock:
            file.write( self, string )
            self.flush()


rep1_cols = numpy.array((0,2,4))
rep2_cols = numpy.array((1,3,5))
def cv_gen():
    yield rep1_cols, rep2_cols
    yield rep2_cols, rep1_cols
    return


def load_data():
    samples = None
    genes = []
    expression = numpy.zeros((103501, 10), dtype=float)
    with open("all_quant.txt") as fp:
        for i, line in enumerate(fp):
            data = line.split()
            if i == 0:
               samples = data[1:]
               continue
            genes.append(data[0])
            expression[i-1,:] = map(float, data[1:])
    
    #res = numpy.zeros((103501, 103501))
    #cov = graph_lasso(expression.dot(expression.T), 100)
    #cov = expression.dot(expression.T)
    #numpy.save("cov", cov)
    #cov = numpy.load("cov.npy")
    cov = None
    return samples, genes, expression[:,(3,4,0,2,5,7)], cov

def estimate_covariates(pred_expression, resp_expression, i):
    #if ( pred_expression[i,(0,2,4)].sum() < 1 #1e-12 
    #     or pred_expression[i,(1,3,5)].sum() < 1): # 1e-12):
    #    nonzero_coefs = []
    #    alpha = 0
    #else:
    response = resp_expression[i,:]
    pred_expression[i,:] = 0
    clf = linear_model.LassoCV(
       fit_intercept=False, cv=cv_gen(), verbose=0, n_jobs=1)
    clf.fit(pred_expression.T, response.T)
    alpha = clf.alpha_
    nonzero_coefs = numpy.nonzero(numpy.abs(clf.coef_) > 1e-12)[0].tolist()
    
    return alpha, nonzero_coefs

def fit_regression_models(expression, expression_indices):
    pred_expression = expression[:,(0,1,2,3,4,5)]
    resp_expression = expression[:,(0,1,2,3,4,5)]
    
    shared_index = Value('i', 0)
    pids = []
    with ThreadSafeFile("output_noshift.txt", "w") as ofp:
        for p_index in xrange(NTHREADS):
            pid = os.fork()
            if pid == 0:
                while True:
                    with shared_index.get_lock():
                        i = shared_index.value
                        if i >= len(expression): break
                        shared_index.value += 1
                    alpha, nonzero_coefs = estimate_covariates(
                        pred_expression, resp_expression, i)
                    output_str = "{}\t{}\t{}\n".format(
                        expression_indices[i], alpha, 
                        "\t".join(str(expression_indices[x]) 
                                  for x in nonzero_coefs))
                    print output_str,
                    ofp.write(output_str)
                sys._exit()
            else:
                pids.append(pid)
        try:
            for pid in pids:
                os.waitpid(pid, 0)
        except:
            for pid in pids:
                os.kill(pid, signal.SIGTERM)
            raise

def build_graph(edge_fname, expression, genes):
    G = nx.Graph()
    for node_id, gene in enumerate(genes):
        G.add_node(node_id, gene_id=genes[node_id], alpha=0.0)
    
    # add all of the node alphas and all edges
    with open(edge_fname) as fp:
        for line in fp:
            data = line.split()
            node_id = int(data[0])
            alpha = float(data[1])
            G.node[node_id]['alpha'] = alpha
            for edge_target in (int(x) for x in data[2:]):
                G.add_edge(node_id, edge_target)
    

    for x in nx.connected_component_subgraphs(G):
        if len(x.nodes()) <= 10:
            for node in x.nodes():
                G.remove_node(node)
    return G

def main():
    sample, genes, raw_expression, cov = load_data()
    expression = raw_expression[raw_expression.sum(1) > 6]
    expression_indices = numpy.nonzero(
        raw_expression.sum(1) > 6)[0].tolist()
    expression = numpy.log10(expression + 1)
    #fit_regression_models(expression, expression_indices)
    G = build_graph(sys.argv[1], expression, genes)
    #nx.write_gml(G, edge_fname + ".gml")
    print G
    return

def test_main():
    G = nx.balanced_tree(0.1, 3)
    nx.draw(G)
    return
    import matplotlib.pyplot as plt
    pos=nx.graphviz_layout(G,prog="neato")
    
    #nx.write_gml(G, edge_fname + ".gml")
    print G
    return

    
main()









# tmp
