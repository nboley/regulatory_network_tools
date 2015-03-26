import os, sys

import math

from collections import defaultdict, OrderedDict

from itertools import combinations

import multiprocessing
import signal

import networkx as nx
import matplotlib.pyplot as plt

import random
import numpy.random
import numpy

from numpy.linalg import pinv

import numpy as np
from scipy import stats, linalg

from sklearn import linear_model, preprocessing
from sklearn.feature_selection import f_regression
from scipy.linalg import lstsq
from scipy.stats import f, norm, pearsonr

N_THREADS = 32

#print numpy.random.seed()
try: 
    numpy.random.seed(int(sys.argv[1]))
except:
    seed = random.randrange(100)
    print "SEED:", seed
    numpy.random.seed(seed)
# 82 is a good seed 

DEBUG_VERBOSE = False 
VERBOSE = True 
REVERSE_CAUSALITY = True

ALPHA = 0.01
MAX_ORDER = 8
MIN_TPM = 2

def partial_corr(C):
    inv_cov = pinv(numpy.cov(C))
    normalization_mat = numpy.sqrt(
        numpy.outer(numpy.diag(inv_cov), numpy.diag(inv_cov)))
    return inv_cov/normalization_mat


def build_causal_graph(n_children, max_height):
    def add_children(G, parent, num):
        child_level = G.node[parent]['level'] + 1
        children_ids = []
        for i in xrange(num):
            id = "%i_%i_%i" % (child_level, G.node[parent]['order'], i)
            G.add_node(id, level=child_level, order=i, parent_order=G.node[parent]['order'])
            if not REVERSE_CAUSALITY:
                G.add_edge(parent, id)
            else:
                G.add_edge(id, parent)
            yield id
        return
    
    G = nx.DiGraph()
    G.add_node("0", level=0, order=0)
    pending_nodes = ["0",]
    while len(pending_nodes) > 0:
        curr_node = pending_nodes.pop()
        if G.node[curr_node]['level'] < max_height:
            pending_nodes.extend(list(add_children(G, curr_node, n_children)))
    return G

def simulate_causal_graph(depth, n_children):
    # build a causal graph
    return build_causal_graph(n_children, depth)


def simulate_data_from_causal_graph(G, n_timepoints, corr=0.9):
    # find all nodes without successors
    assert nx.is_directed_acyclic_graph(G)
    nodes_stack = [node for node in G.nodes() 
                   if len(G.predecessors(node)) == 0]
    expression_values = dict(
        (node, numpy.random.randn(n_timepoints)) 
        for node in nodes_stack )
    
    # deal with the successors
    while len(nodes_stack) > 0:
        parent = nodes_stack.pop()
        for child in G.successors(parent):
            val = ( (1-corr)*numpy.random.randn(n_timepoints) 
                    + corr*expression_values[parent] )
            val = val/len(G.predecessors(child))
            if child not in expression_values:
                expression_values[child] = val
            else:
                expression_values[child] += val
            nodes_stack.append(child)

    return ( sorted(expression_values), 
             numpy.vstack(y for x, y in sorted(expression_values.items())))

def estimate_covariates(sample1, sample2, resp_index, alpha=0.50):
    def cv_gen():
        n_tps = sample1.shape[1]
        yield range(n_tps), range(n_tps, 2*n_tps)
        yield range(n_tps, 2*n_tps), range(n_tps)
        return
    
    merged_samples = numpy.hstack((sample1, sample2))
    # normalize to sum 1
    merged_samples = ((merged_samples.T)/(merged_samples.sum(1))).T

    pred_mask = numpy.ones(merged_samples.shape[0], dtype=bool)
    pred_mask[resp_index] = False

    # filter marginally uncorrelated covariates
    response = merged_samples[resp_index,:]    
    best_i, min_pval = None, 10
    for i, row in enumerate(merged_samples):
        if i == resp_index: continue
        corr_coef, p_value = pearsonr(row, response)
        if p_value > alpha: pred_mask[i] = False
        if p_value < min_pval: best_i = i
    pred_mask[best_i] = True
    
    clf = linear_model.LassoCV(
       fit_intercept=True, cv=cv_gen(), verbose=0, max_iter=100000)
    clf.fit(merged_samples[pred_mask,:].T, response.T)
    
    alpha = clf.alpha_
    
    clf = linear_model.Lasso(
       alpha=alpha, max_iter=100000)
    clf.fit(merged_samples[pred_mask,:].T, response.T)

    coefs = clf.coef_.tolist()
    for i in sorted(int(x) for x in (1-pred_mask).nonzero()[0]):
        coefs.insert(i, 0)
    return alpha, coefs


def estimate_skeleton_from_samples(sample1, sample2, labels, thresh_ratio=1000):
    G = nx.DiGraph()
    for i in xrange(sample1.shape[0]):
        G.add_node(i, label=labels[i])
        alpha, regression_coefs = estimate_covariates(sample1, sample2, i)
        max_coef = max(numpy.abs(regression_coefs))
        for j, val in enumerate(regression_coefs):
            if abs(val) > 1e-6 and abs(val) >= max_coef/thresh_ratio:
                G.add_edge(i, j, weight=val)
    
    # remove non bi-directed edges and set the edge weights
    for a, b in list(G.edges()):
        if not G.has_edge(b, a): 
            G.remove_edge(a, b)
        else:
            G[a][b]['weight'] = min(G[a][b]['weight'], G[b][a]['weight'])
    
    return G.to_undirected()

def estimate_marginal_correlations(merged_samples, resp_index, alpha=ALPHA,
                                   min_num_neighbors=0, max_num_neighbors=10000 ):
    sig_samples = []
    response = merged_samples[resp_index,:]
    rep1 = numpy.array((0,1,2))
    rep2 = numpy.array((3,4,5))
    for i, row in enumerate(merged_samples):
        if i <= resp_index: continue
        corr_coef1, p_value1 = pearsonr(row[rep1], response[rep2])
        corr_coef2, p_value2 = pearsonr(row[rep2], response[rep1])
        if p_value1 < p_value2:
            p_value = p_value1
            corr_coef = corr_coef1
        else:
            p_value = p_value2
            corr_coef = corr_coef2
        #print p_value, corr_coef
        if p_value < alpha or len(sig_samples) < min_num_neighbors: 
            sig_samples.append((p_value, corr_coef, i))
    sig_samples.sort(key=lambda x: (x[0], -x[1]))
    while (len(sig_samples) > min_num_neighbors 
           and sig_samples[-1][0] > alpha)/len(merged_samples):
        del sig_samples[-1]
    return sig_samples[:max_num_neighbors]

def estimate_initial_skeleton(normalized_data, labels, alpha):
    print "Estimating marginal independence relationships"
    # initialize a graph storing the covariance structure
    manager = multiprocessing.Manager()
    edges_and_data = manager.list()
    curr_node = multiprocessing.Value('i', 0)
    n_nodes = normalized_data.shape[0]
    pids = []
    for p_index in xrange(N_THREADS):
        pid = os.fork()
        if pid == 0:
            while True:
                with curr_node.get_lock():
                    node = curr_node.value
                    if node >= n_nodes: break
                    curr_node.value += 1
                if node%100 == 0: 
                    print "O0: Finished processing %i/%i nodes" % (
                        node, n_nodes)
                marginal_dep = estimate_marginal_correlations(
                    normalized_data, node, alpha)
                edges_and_data.append((node, marginal_dep ))
            os._exit(0)
        else:
            pids.append(pid)
    try:
        for pid in pids:
            os.waitpid(pid, 0)
    except:
        for pid in pids:
            os.kill(pid, signal.SIGTERM)
        raise
    
    print "Building O0 graph skeleton from marginal independence data"
    G = nx.Graph()
    for i, data in edges_and_data:
        if not G.has_node(i): G.add_node(i, label=labels[i])
        for p, corr, j in data:
            if not G.has_node(j): G.add_node(j, label=labels[j])
            G.add_edge(i, j, corr=corr, marginal_p=p)
    
    manager.shutdown()
    return G


def find_unoriented_edges(G):
    unoriented_edges = set()
    for start, stop in G.edges():
        if start > stop: continue
        if G.has_edge(stop, start):
            unoriented_edges.add((start, stop))
    return sorted(unoriented_edges)

def iter_unoriented_v_structures(G):
    for node, node_data in G.nodes(data=True):
        neighbors = [neighbor for neighbor in nx.all_neighbors(G, node)
                     if G.has_edge(neighbor, node) and G.has_edge(node, neighbor)]
        for n1 in neighbors:
            for n2 in neighbors:
                if n1 >= n2: 
                    continue
                if G.has_edge(n1, n2) or G.has_edge(n2, n1): 
                    continue
                yield n1, node, n2
    return

def iter_colliders(G):
    for node, node_data in G.nodes(data=True):
        predecessors = [predecessor for predecessor in G.predecessors(node)
                        if not G.has_edge(node, predecessor) 
                        and G.has_edge(predecessor, node)]
        for n1 in predecessors:
            assert G.has_edge(n1, node)
            assert not G.has_edge(node, n1)
            for n2 in predecessors:
                if n1 >= n2: continue
                yield n1, node, n2
    return

def orient_single_v_structure(G, data):
    for a, c, b in iter_unoriented_v_structures(G):
        if not are_cond_indep(a, b, c, data):
            # remove the edges point from c to a and b
            try: G.remove_edge(c, a)
            except: pass
            try: G.remove_edge(c, b)
            except: pass
            
            if VERBOSE:
                print "Orienting %s->%s<-%s" % (
                    G.node[a]['label'], G.node[c]['label'], G.node[b]['label'])
            return True
    
    return False

def orient_v_structures(G, ci_sets):
    for a, c, b in iter_unoriented_v_structures(G):
        # the a and b are marginally independent or they
        # are independent conditionally on c orient the edges
        if (a,b) not in ci_sets or c not in ci_sets[(a,b)]:
            # remove the edges point from c to a and b
            try: G.remove_edge(c, a)
            except: pass
            try: G.remove_edge(c, b)
            except: pass
    return

def apply_rule_1(b, c, G):
    """Check is there is a directed edge a->b such that a and c are not adjacent
    """
    for a in G.predecessors(b):
        # skip bi-directed edges
        if G.has_edge(b, a): 
            continue
        # skip adjacent edges
        if G.has_edge(a, c) or G.has_edge(c, a):
            continue
        return True
    return False

def apply_rule_2(a, b, G):
    """Check is there is a directed edge a->b such that a and b are not adjacent
    """
    for c in G.successors(a):
        # skip bi-directed edges
        if G.has_edge(c, a): 
            continue
        # if there also exists a directed edge from c,b 
        # then rule 2 applies
        if G.has_edge(c, b) and not G.has_edge(b, c):
            return True
    return False

def apply_rule_3(a, b, G):
    """Check is there are two chains a--c->b and a--d->b such that 
       c and d are not adjacent.
    """
    intermediate_nodes = set()
    # find all such chains
    for c in G.successors(a):
        # skip bi-directed edges
        if G.has_edge(c, a) and G.has_edge(c, b) and not G.has_edge(b, c): 
            intermediate_nodes.add(c)

    # try to find a pair of non adjacent nodes
    for c in intermediate_nodes:
        for d in intermediate_nodes:
            if c == d: continue
            if not G.has_edge(c, d) and not G.has_edge(d, c):
                return True
    return False

def apply_rule_4(a, b, G):
    """Check is there is a chain a--c->d->b where c and b are non adjacent.
    """
    # find all such chains
    for c in G.successors(a):
        # skip non bi-directed edges
        if not G.has_edge(c, a):
            continue
        # skip c that are adjacent to b
        if G.has_edge(c, b) or G.has_edge(b, c):
            continue
        for d in G.successors(c):
            # skip bi directed edges
            if G.has_edge(d, c): continue
            if G.has_edge(d, b) and not G.has_edge(b, d):
                return True
    return False

def apply_IC_rules(G):
    unoriented_edges = find_unoriented_edges(G)
    for a, b in unoriented_edges:
        if apply_rule_1(a, b, G):
            if DEBUG_VERBOSE: print "Applying Rule 1:", a, b
            G.remove_edge(b ,a)
            return True
        elif apply_rule_2(a, b, G):
            if DEBUG_VERBOSE: print "Applying Rule 2:", a, b
            G.remove_edge(b ,a)
            return True
        elif apply_rule_3(a, b, G):
            if DEBUG_VERBOSE: print "Applying Rule 3:", a, b
            G.remove_edge(b ,a)
            return True
        elif apply_rule_4(a, b, G):
            if DEBUG_VERBOSE: print "Applying Rule 4:", a, b
            G.remove_edge(b ,a)
            return True
    return False


def break_cycles(G):
    G = G.copy()
    cycles = nx.cycle_basis(G)
    while len(cycles) > 0:
        # break cycles
        edge_cnts = defaultdict(int)
        edge_weights = {}
        for cycle in nx.cycle_basis(G):
            for a,b in zip(cycle[:-1], cycle[1:]) + [(cycle[-1],cycle[0]),]:
                a, b = min(a,b), max(a,b)
                edge_cnts[(a,b)] += 1
                edge_weights[(a,b)] = G[a][b]['weight']
        edge = sorted(
            edge_cnts, key=lambda x: (-edge_cnts[x], abs(edge_weights[x])))[0]
        #print edge, edge_cnts[edge], abs(edge_weights[edge])
        #print "Removing edge:", edge
        G.remove_edge(*edge)
        cycles = nx.cycle_basis(G)
    return G

def test_for_CI(G, n1, n2, normalized_data, order, alpha):
    """Test if n1 and n2 are conditionally independent. 

    If they are not return None, else return the conditional independence set.
    """
    N = normalized_data.shape[1]
    ones = numpy.ones((N,1), dtype=float)
    
    n1_resp = normalized_data[n1,:]    
    n1_neighbors = set(G.neighbors(n1))    
    
    n2_resp = normalized_data[n2,:]
    n2_neighbors = set(G.neighbors(n2))
    
    common_neighbors = n1_neighbors.intersection(n2_neighbors) - set((n1, n2))
    # if there aren't enough neighbors common to n1 and n2, return none
    if len(common_neighbors) < order: 
        return None
    
    min_score = 1e100
    best_p_val = None
    best_neighbors = None
    n_common_neighbors = 0
    for covariates in combinations(common_neighbors, order):
        n_common_neighbors += 1
        predictors = numpy.hstack(
            (ones, normalized_data[numpy.array(covariates),:].T))
        # test if node is independent of neighbors given for some subset
        rv1, _, _, _ = lstsq(predictors, n1_resp)
        rv2, _, _, _ = lstsq(predictors, n2_resp)
        cor, pval =  pearsonr(n1_resp - rv1.dot(predictors.T), 
                              n2_resp - rv2.dot(predictors.T))
        if abs(cor) < min_score:
            min_score = abs(cor)
            best_neighbors = covariates
            best_p_val = pval
        # make the multiple testing correction /n_common_neighbors
        if best_p_val < alpha/n_common_neighbors:
            return None
        #score = math.sqrt(N-order-3)*0.5*math.log((1+cor)/(1-cor))
        #print abs(score),  norm.isf(alpha/(len(neighbors)*2)), cor, pval
        #if abs(score) < norm.isf(alpha/(len(neighbors)*2)): 
    
    # make the multiple testing correction /n_common_neighbors
    if best_p_val < alpha/n_common_neighbors:
        return None
    else:
        return best_neighbors

def apply_pc_iteration_serial(G, normalized_data, order, alpha=ALPHA):    
    cond_independence_sets = defaultdict(set)
    
    ## we can't estimate higher order interactiosn than we have samples
    #N = normalized_data.shape[1]
    #if N - order - 3 <= 0:
    #    return cond_independence_sets

    num_nodes = len(G.nodes())
    for n1 in G.nodes():    
        n1_neighbors = sorted(
            G.neighbors(n1), key=lambda n2: G[n1][n2]['marginal_p'])
        num_neighbors = len(n1_neighbors)
        if DEBUG_VERBOSE:
            print "Test O%i: %i/%i %i neighbors ... " % (
                order, n1, num_nodes, num_neighbors),

        for n2 in n1_neighbors:
            if n2 <= n1: continue
            are_CI = test_for_CI(G, n1, n2, normalized_data, order, alpha)
            if are_CI == None:
                if DEBUG_VERBOSE: print "%i NOT CI of %i" % (n1, n2)
            else:
                if DEBUG_VERBOSE:
                    print "%i IS CI %i | %s" % (
                        n1, n2, ",".join(str(x) for x in are_CI))
                cond_independence_sets[(n1, n2)].update(are_CI)
                G.remove_edge(n1, n2)
                num_neighbors -= 1
        if VERBOSE: 
            print "%i remain (%i removed)" % (
                num_neighbors, len(n1_neighbors) - num_neighbors)
    
    return cond_independence_sets

def remove_edges_for_single_node(G, n1, normalized_data, order, alpha):
    # the edges removed, and their corresponding conditional independence sets
    cond_independence_sets = defaultdict(set)
    
    # n1's neighbors, sorted by decreasing marginal p value. Since the algorithm
    # is order depedent, all else being equal we prefer to remove edges that 
    # have the lowest marginal correlation first
    n1_neighbors = sorted(
        G.neighbors(n1), key=lambda n2: G[n1][n2]['marginal_p'])
    num_neighbors = len(n1_neighbors)

    for n2 in n1_neighbors:
        if n2 <= n1: continue
        if num_neighbors-1 <= order: break
        are_CI = test_for_CI(G, n1, n2, normalized_data, order, alpha)
        if are_CI == None:
            if DEBUG_VERBOSE: print "%i NOT CI of %i" % (n1, n2)
        else:
            if DEBUG_VERBOSE:
                print "%i IS CI %i | %s" % (
                    n1, n2, ",".join(str(x) for x in are_CI))
            cond_independence_sets[(n1, n2)].update(are_CI)
            num_neighbors -= 1
    if DEBUG_VERBOSE: 
        print "Test O%i: %i/%i %i neighbors ... %i remain (%i removed)" % (
            order, n1, len(G), len(n1_neighbors), 
            num_neighbors, len(n1_neighbors) - num_neighbors)
    
    return cond_independence_sets

def find_edges_to_remove(G, normalized_data, order, alpha=ALPHA):    
    cond_independence_sets = {}
    
    num_nodes = len(G.nodes())
    for n1 in G.nodes():    
        edges_to_remove = find_edges_to_remove(
            G, n1, normalized_data, order, alpha)
        for edge, ci_sets in cond_independence_sets:
            assert edge not in cond_independence_sets
            cond_independence_sets[edge] = ci_sets
    
    return cond_independence_sets

def partition_nodes_into_nonadjacent_sets(G):
    grouped_nodes = set()
    nodes_sets = []
    degree_sorted_nodes = sorted(G.degree().items(), key=lambda x:-x[1])
    while len(grouped_nodes) < len(G):
        nodes_sets.append( [] )
        remaining_nodes = set(G.nodes())
        for node, degree in degree_sorted_nodes:
            # skip nodes htat we've already grouped
            if node in grouped_nodes: continue
            # skip nodes that are neighbors to nodes already in the set
            if node not in remaining_nodes: continue
            # add the node to the latest set, and mark it as having
            # been added
            nodes_sets[-1].append(node)
            grouped_nodes.add(node)
            # remove neighbots of this node from the set of nodes that
            # can still be considered
            remaining_nodes.remove(node)
            remaining_nodes.difference_update(G.neighbors(node))
    return nodes_sets

def find_CI_relationships_in_subprocess(
        normalized_data, skeleton, ind_order, alpha,
        node_set_queue, shared_edges_to_remove):
    # make a thread local copy of the graph (shouldnt be necessary on
    # a fork, but to be safe... )
    curr_skeleton = skeleton.copy()
    while True:
        node = node_set_queue.get()
        # if there are no nodes left, break
        if node == None: break
        node_CI_sets = remove_edges_for_single_node(
            curr_skeleton, node, normalized_data, ind_order, alpha)
        for edge, CI_set in node_CI_sets.iteritems():
            shared_edges_to_remove.put((edge, CI_set))
    return

def estimate_pdag(sample1, sample2, labels, alpha=ALPHA):
    # combine and normalize the samples
    normalized_data = numpy.hstack((sample1, sample2))
    normalized_data = ((normalized_data.T)/(normalized_data.sum(1))).T
    
    skeleton = estimate_initial_skeleton(
        normalized_data, labels, alpha=alpha)
    nx.write_gml(skeleton, "skeleton_O0.gml")

    cond_independence_sets = {}
    for ind_order in xrange(1,min(MAX_ORDER,min(normalized_data.shape)-2+1)):
        nonadjacent_node_sets = partition_nodes_into_nonadjacent_sets(skeleton)
        print "O%i: Processing cluster %i/%i: %i/%i nodes remain" % (
            ind_order, 1, len(nonadjacent_node_sets),
            sum(len(x) for i, x in enumerate(nonadjacent_node_sets) 
                if i >= 0),
            sum(len(x) for x in nonadjacent_node_sets))

        for cluster_i, nonadjacent_nodes in enumerate(nonadjacent_node_sets):
            new_cond_independence_sets = multiprocessing.Queue()
            # populate the queue
            n_threads_to_spawn = min(N_THREADS, len(nonadjacent_nodes))
            nodes_queue = multiprocessing.Queue()
            for node in nonadjacent_nodes: nodes_queue.put(node)
            for i in xrange(n_threads_to_spawn): nodes_queue.put(None)

            # spaqwn worker processes
            pids = []
            for i in xrange(n_threads_to_spawn):
                pid = os.fork()
                if pid == 0:
                    find_CI_relationships_in_subprocess(
                        normalized_data, skeleton, ind_order, alpha, 
                        nodes_queue, new_cond_independence_sets)
                    os._exit(0)
                else:
                    pids.append(pid)
            try:
                for pid in pids:
                    os.waitpid(pid, 0)
            except:
                for pid in pids:
                    os.kill(pid, signal.SIGTERM)
                raise
            print "O%i: Finished processing cluster %i/%i: %i/%i nodes remain" % (
                ind_order, cluster_i+1, len(nonadjacent_node_sets),
                sum(len(x) for i, x in enumerate(nonadjacent_node_sets) 
                    if i > cluster_i),
                sum(len(x) for x in nonadjacent_node_sets))
            
            # update the global cond_independence_sets
            while not new_cond_independence_sets.empty():
                edge, CI_set = new_cond_independence_sets.get()
                assert edge not in cond_independence_sets
                cond_independence_sets[edge] = CI_set
                skeleton.remove_edge(*edge)

        print "Writing O%i skeleton to disk" % ind_order        
        nx.write_gml(skeleton, "skeleton_O%i.gml" % ind_order)

    
    est_G = skeleton.to_directed()
    orient_v_structures(est_G, cond_independence_sets)
    applied_rule = True
    while apply_IC_rules(est_G): pass
    return est_G, cond_independence_sets

def hierarchical_layout(real_G):
    level_grouped_nodes = defaultdict(list)
    for node, data in real_G.nodes(data=True):
        level_grouped_nodes[data['level']].append(node)
    pos = {}
    for level, nodes in level_grouped_nodes.iteritems():
        y_pos = 1 - float(level+1)/(len(level_grouped_nodes)+1)
        for i, node in enumerate(sorted(
                nodes, key=lambda id: [int(x) for x in id.split("_")])):
            x_pos = float(i+1)/(len(nodes)+1)
            pos[node] = numpy.array((x_pos, y_pos+(x_pos/1.7)**2))
    
    return pos

def iter_all_subsets_of_siblings(G, node):
    siblings = set(nx.all_neighbors(G, node))
    for i in xrange(1, len(siblings)+1):
        for subset in combinations(siblings, i):
            yield subset
    return

def brute_force_find_all_consistent_dags(pdag, data):
    def get_undirected_edges(G):
        undirected_edges = []
        for a, b in G.edges():
            if G.has_edge(b, a):
                undirected_edges.append( (a,b) )
        return undirected_edges
    
    def orient_edge_and_propogate_changes(G, a, b):
        G = G.copy()
        G.remove_edge(b, a)
        #orient_v_structures(G, data)
        while apply_IC_rules(G): pass
        return G

    # the dags that we've found
    dags = []
    # stack containing both edges that we have already oriented, and 
    # the current pdag after those edges have been oriented
    stack = [ [[(a,b),], orient_edge_and_propogate_changes(pdag, a, b)]
              for a, b in get_undirected_edges(pdag) ]
    while len(stack) > 0:
        oriented_edges, curr_pdag = stack.pop()
        edges_to_orient = get_undirected_edges(curr_pdag)
        # if there are no more edges to orient, then make sure that we haven't 
        # already seen it and that it is acyclic, and add it
        if len(edges_to_orient) == 0: 
            if (not any(set(curr_pdag.edges()) == set(x.edges()) for x in dags)
                    and nx.is_directed_acyclic_graph(curr_pdag)):
                dags.append(curr_pdag)
        else:
            for a, b in edges_to_orient:
                new_pdag = orient_edge_and_propogate_changes(curr_pdag, a, b)
                stack.append([oriented_edges + [(a, b),], new_pdag])
    
    return dags


def find_all_consistent_dags(pdag, cond_independence_sets, real_G):
    initial_colliders = set(iter_colliders(pdag))
    
    def get_undirected_edges(G):
        undirected_edges = []
        for a, b in G.edges():
            if G.has_edge(b, a):
                undirected_edges.append( (a,b) )
        return undirected_edges
    
    def orient_edge_and_propogate_changes(G, a, b):
        G = G.copy()
        G.remove_edge(b, a)
        orient_v_structures(G, cond_independence_sets)
        while apply_IC_rules(G): pass
        # make sure that the resulting graph is acyclic
        #if not nx.is_directed_acyclic_graph(G): return None
        if any( collider not in initial_colliders 
                for collider in  iter_colliders(G) ): 
            return None
        return G

    #unoriented_v_structures = list(iter_unoriented_v_structures(pdag))
    # the dags that we've found
    dags = []
    # stack containing both edges that we have already oriented, and 
    # the current pdag after those edges have been oriented
    stack = [ [[(a,b),], orient_edge_and_propogate_changes(pdag, a, b)]
              for a, b in get_undirected_edges(pdag) ]
    stack = [(x, g) for x, g in stack if g != None]
    while len(stack) > 0:
        oriented_edges, curr_pdag = stack.pop()
        edges_to_orient = get_undirected_edges(curr_pdag)
        # if there are no more edges to orient, then make sure that we haven't 
        # already seen it and that it is acyclic, and add it
        if len(edges_to_orient) == 0: 
            if (not any(set(curr_pdag.edges()) == set(x.edges()) for x in dags)
                    and nx.is_directed_acyclic_graph(curr_pdag)):
                #plot_pdag(curr_pdag, real_G)
                dags.append(curr_pdag)
        else:
            for a, b in edges_to_orient:
                new_pdag = orient_edge_and_propogate_changes(curr_pdag, a, b)
                if new_pdag != None:
                    stack.append([oriented_edges + [(a, b),], new_pdag])
    
    return dags


def plot_pdag(pdag, real_G):
    real_G_layout = hierarchical_layout(real_G)
    #nx.draw(est_G, nx.graphviz_layout(est_G,prog='twopi',args=''))
    labels = dict((id, data['label']) for id, data in pdag.nodes(data=True))
    pos = dict((id, real_G_layout[data['label']]) 
               for id, data in pdag.nodes(data=True))
    nx.draw(pdag, pos, labels=labels, node_size=1500, node_color='white')
    plt.show()
    return


def load_data():
    samples = None
    genes = []
    expression = []
    with open("all_quant.txt") as fp:
        for i, line in enumerate(fp):
            data = line.split()
            if i == 0:
               samples = data[1:]
               continue
            # append the gene name
            gene_expression = numpy.array(map(float, data[1:]))
            if numpy.max(gene_expression) < MIN_TPM: continue
            genes.append(data[0])
            # add random normal noise to the gene expressionvalues to prevent 
            # high correlation artifacts due to rounding error, etc. 
            gene_expression += (numpy.random.randn(len(gene_expression)))**2
            expression.append( gene_expression )
    
    expression = numpy.vstack(expression)
    s1 = expression[:,(3,0,5)]
    s2 = expression[:,(4,2,7)]
    return genes, s1, s2

def main():
    genes, sample1, sample2 = load_data()
    pdag, cond_independence_sets = estimate_pdag(sample1, sample2, genes)
    nx.write_gml(pdag, "expression_GT_%i_pdag.gml" % MIN_TPM)
    print pdag
    return

def test():
    real_G = simulate_causal_graph(5, 5)
    #nx.draw(real_G, hierarchical_layout(real_G), node_size=1500, node_color='blue')
    #plt.show()
    #return

    labels, sample1 = simulate_data_from_causal_graph(real_G, 30, 0.5)
    labels, sample2 = simulate_data_from_causal_graph(real_G, 30, 0.5)
    #print partial_corr(preprocessing.scale(numpy.hstack((sample1, sample2)).T).T)
    #return
    pdag, cond_independence_sets = estimate_pdag(sample1, sample2, labels)
    plot_pdag(pdag, real_G)
    #for x in find_all_consistent_dags(
    #        pdag, cond_independence_sets, real_G):
    #    print x.edges()
    #    plot_pdag(x, real_G)
    
    return


if __name__ == '__main__':
    main()
    #test()
