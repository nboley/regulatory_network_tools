import os, sys, signal

import multiprocessing
from multiprocessing.sharedctypes import Value, Array

import numpy

from sklearn.covariance import GraphLasso, GraphLassoCV
from sklearn import linear_model
from sklearn.cross_validation import KFold
from sklearn import preprocessing

NTHREADS = 24

class ThreadSafeFile( file ):
    def __init__( *args ):
        file.__init__( *args )
        args[0].lock = multiprocessing.Lock()

    def write( self, string ):
        with self.lock:
            file.write( self, string )
            self.flush()


rep1_cols = numpy.array((0,2))
rep2_cols = numpy.array((1,3))
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
    return samples, genes, expression, cov

def estimate_covariates(pred_expression, resp_expression, i):
    if ( pred_expression[i,(0,2)].sum() < 1 #1e-12 
         or pred_expression[i,(1,3)].sum() < 1): # 1e-12):
        nonzero_coefs = []
        alpha = 0
    else:
        response = resp_expression[i,:]
        pred_expression[i,:] = 0
        clf = linear_model.LassoCV(
            fit_intercept=False, cv=cv_gen(), verbose=0, n_jobs=1)
        clf.fit(pred_expression.T, response.T)
        alpha = clf.alpha_
        nonzero_coefs = numpy.nonzero(numpy.abs(clf.coef_) > 1e-12)[0].tolist()
    
    return alpha, nonzero_coefs

def main():
    sample, genes, raw_expression, cov = load_data()
    expression = raw_expression[raw_expression.min(1) > 100]
    expression_indices = numpy.nonzero(raw_expression.sum(1) > 6)[0].tolist()
    
    ## reorder and filter data
    #rep1_cols = numpy.array((3,0,5)) # 8 is co culture
    #rep2_cols = numpy.array((4,2,7)) # 9 is MRC5
    expression = expression[:,(3,4,0,2,5,7)]

    # log data
    expression = numpy.log10(expression + 1)[1:100,]
    cov = expression.dot(expression.T)
    print cov.shape
    #mo = GraphLasso(alpha=95, mode='lars', verbose=True) #, cv=KFold(3,2), n_jobs=24)
    mo = GraphLassoCV(mode='lars', verbose=True, cv=KFold(3,2), n_jobs=24)
    sparse_cov = mo.fit(cov)
    print( numpy.nonzero(sparse_cov)[0].sum() )
    return

main()









# tmp
