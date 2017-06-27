#!/usr/local/bin/python3

# Copyright 2016 Maxwell Shinn (mws41@cam.ac.uk)
# Available under the GNU GPLv3.

# If you use this code, please cite the associated paper.

# This script contains a number of useful functions for analyzing
# brain networks, especially functional networks.  It is written in
# Python3, and can only be guaranteed to work there.  

# Dependencies:

# - networkx
# - scipy/numpy/matplotlib
# - bctpy: The module "bct" is bctpy, a port of the Brain Connectivity
#   Toolbox to Python.  The latest version supports Python3.
#   https://github.com/aestrivex/bctpy or "pip install bcpty".  If you
#   don't want to install bctpy, it should be pretty easy to modify
#   this code to remove the dependency.

# NOTE: Functions are documented with inline comments, not in docstrings.

# Here is a quick example to get you started:
#
#     import networkx
#     from versatility import *
#     G = networkx.karate_club_graph()
#     find_nodal_mean_versatility(G, find_communities_louvain, processors=2)
#     find_nodal_versatility(G, find_communities_louvain, algargs={"gamma" : 0.5})

# You've probably noticed that there are a lot of assert statements.
# This code will run a bit faster if you run "python -O", which will
# skip the asserts.  I can't say I recommend this, but if you really
# want to do so then be my guest...

# # The following line may make this run in Python2, though your
# # results may vary.
# from __future__ import print_function, unicode_literals, absolute_import, division


import bct
import networkx
import matplotlib.pyplot as plt
import numpy
import scipy.stats

# ▌ ▌            ▐  ▗▜ ▗▐             ▐     ▗       
# ▚▗▘▞▀▖▙▀▖▞▀▘▝▀▖▜▀ ▄▐ ▄▜▀ ▌ ▌ ▛▚▀▖▞▀▖▜▀ ▙▀▖▄ ▞▀▖▞▀▘
# ▝▞ ▛▀ ▌  ▝▀▖▞▀▌▐ ▖▐▐ ▐▐ ▖▚▄▌ ▌▐ ▌▛▀ ▐ ▖▌  ▐ ▌ ▖▝▀▖
#  ▘ ▝▀▘▘  ▀▀ ▝▀▘ ▀ ▀▘▘▀▘▀ ▗▄▘ ▘▝ ▘▝▀▘ ▀ ▘  ▀▘▝▀ ▀▀ 

# find_nodal_versatility computes the versatility of each node.  This
# replaces find_nodal_uncertainty and find_nodal_uncertaintyplusplus
# as measures of how difficult it is to group a node into a module.
# `g` should be a directed/undirected, weighted/unweighted graph, and
# `alg` should be a function compatible with the "alg" argument of the
# "consensus_matrix" function. The definition of nodal versatility is:
# 
# $U_n = ∑_{i∈N} sin(π a_{n,i})/(∑_{i∈N} a_{n,i})$
# 
# where $a_{i,j}$ is an element of the association matrix, i.e. each
# cell is the probability that algorithm `alg` will classify those
# nodes in the same module.  This assigns the property "[algname]vers"
# to each node, where "[algname]" is `algname`.  It also assigns the
# mean uncertainty in all nodes to the graph property
# "[algname]meanvers". This returns a dictionary of nodal versatility
# indexed by node.
#
# The optional `algargs` parameter is a dictionary that specifies what
# extra arguments should be passed to `alg`.  (This makes it so that
# you don't have to use partial functions for something as simple as a
# resolution parameter.)  For example, to specify the resolution
# parameter for the louvain algorithm to 1.5, set algargs={"gamma" :
# 1.5}.
#
# The `it` argument is the number of times to run the modularity
# algorithm when estimating the association matrix.  It should
# generally not be lower than 100, and there is little need to make it
# higher than 1000.  The default of 200 is a good tradeoff between
# precision and speed.
#
# The optional argument `processors` should be set to an integer
# greater than 0.  If it is greater than 1, consensus_matrix_par will
# be used instead of consensus_matrix.
def find_nodal_versatility(g, alg, algname="", processors=1, algargs={}, it=200):
    assert type(g) == networkx.classes.graph.Graph, "Not a graph"
    assert type(processors) == int and processors > 0, "Invalid number of processors"
    # We don't check `alg`, `it`, or `algargs` because
    # consensus_matrix function does that for us, and we just pass it
    # directly to there unmodified.
    if processors == 1:
        C = consensus_matrix(g, alg, it=it, algargs=algargs)
    else:
        C = consensus_matrix_par(g, alg, it=it, processors=processors, algargs=algargs)
    g.graph['%sconsmatrix' % algname] = C.astype('float16')
    Cs = numpy.sin(numpy.pi*C)
    assert numpy.all(C == C.T) and numpy.all(Cs == Cs.T), "Assocation matrix or versatility matrix not symmetric"
    assert type(C) == numpy.ndarray and type(Cs) == numpy.ndarray, "Not ndarrays" # Not numpy.matrix 
    versatility = numpy.sum(Cs, axis=0)/numpy.sum(C, axis=0)
    versatility[versatility<1e-10] = 0 # Prevent really small values
    versatilitydict = {g.nodes()[i] : versatility[i] for i in range(0, len(g.nodes()))}
    networkx.set_node_attributes(g, "%svers" % algname, versatilitydict)
    g.graph["%svers" % algname] = numpy.mean(versatility)
    return versatilitydict


# find_nodal_mean_versatility computes the versatility across a
# spectrum of different parameters (e.g. gamma in the louvain
# algorithm) and computes the mean for each node.  `g`, `alg`,
# `algname`, `processors`, and `it` are as documented in
# `find_nodal_versatility`.  `argname` is the argument which we are to
# vary.  `argvals` are the values across which the mean is to be
# taken.  This assigns the value "[algname]meanvers" to each node in
# `g` corresponding to the mean versatility across the spectrum of
# parameters.  It also returns these values as a dictionary indexed by
# node.  Furthermore, it computes the mean of these mean versatility
# values and assigns it as a graph property named "[algname]meanvers".
_argvalsm = numpy.array(range(4, 25), dtype=float)/10
def find_nodal_mean_versatility(g, alg, algname="", processors=1, argname="gamma", argvals=_argvalsm, it=100):
    gc = g.copy()
    for v in argvals:
        find_nodal_versatility(gc, alg, algname=str(v), processors=processors, algargs={argname : v}, it=it)
        print(v)
    means = { n : numpy.mean([gc.node[n][str(v)+"vers"] for v in argvals]) for n in gc.nodes() }
    allvals = { n : dict(zip(argvals, [gc.node[n][str(v)+"vers"] for v in argvals])) for n in gc.nodes() }
    networkx.set_node_attributes(g, "%smeanvers" % algname, means)
    networkx.set_node_attributes(g, "%smeanversvals" % algname, allvals)
    g.graph["%smeanvers" % algname] = numpy.mean(list(means.values()))
    return means

# find_optimal_gamma_curve plots the mean network versatility across a
# spectrum of gammas.  `g` should be a networkx network.  All other
# arguments are optional, and may be specified as in
# find_nodal_versatility.  The algorithm `alg` should also take an
# extra argument, specified by `argname`.  This should be a string
# describing an algorithm argument which can vary across a spectrum.
# It will vary according to the argument `argvals`.  These arguments
# are similar to those in find_nodal_mean_versatility, but instead of
# taking the average, it will plot them.  This returns a list of the
# versatility values.
_argvalsc = numpy.asarray(list(range(0, 40)))/10+.1
def find_optimal_gamma_curve(G, alg, algarg="gamma", argvals=_argvalsc, it=100, show=True, **kwargs):
    import scipy.stats
    import sys
    gs = argvals
    vs = []
    sems = []
    for g in gs:
        v = find_nodal_versatility(G, alg=alg, algargs={algarg : g}, it=it, **kwargs)
        vs.append(numpy.mean(list(v.values())))
        sems.append(scipy.stats.sem(list(v.values())))
        print(g, end=" ")
        sys.stdout.flush()
    print("\n")
    if show == True:
        plt.errorbar(gs, vs, yerr=sems)
        plt.title("Versatility across different values of %s" % algarg)
        plt.xlabel(algarg)
        plt.ylabel("Versatility")
        plt.show()
    return (gs,vs,sems)


# ▞▀▖                                  ▐  ▗          
# ▌  ▞▀▖▛▀▖▞▀▘▞▀▖▛▀▖▞▀▘▌ ▌▞▀▘ ▙▀▖▞▀▖▌ ▌▜▀ ▄ ▛▀▖▞▀▖▞▀▘
# ▌ ▖▌ ▌▌ ▌▝▀▖▛▀ ▌ ▌▝▀▖▌ ▌▝▀▖ ▌  ▌ ▌▌ ▌▐ ▖▐ ▌ ▌▛▀ ▝▀▖
# ▝▀ ▝▀ ▘ ▘▀▀ ▝▀▘▘ ▘▀▀ ▝▀▘▀▀  ▘  ▝▀ ▝▀▘ ▀ ▀▘▘ ▘▝▀▘▀▀ 


# concensus_matrix finds a probability matrix of nodes i and j being
# in the same communitiy.  `g` should be a networkx graph with N
# nodes.  `algorithm` should be a function that takes a networkx graph
# as its input and gives a dictionary as its output, where the index
# is the node and the value is an identifier representing the
# community.  This function runs `algorithm` `it` times (where `it` ∈
# ℕ^+) and returns a NxN array, where the (i,j)-th cell is the
# probability that node i is in the same community as node j.  The
# rows and columns of the matrix are sorted by the list g.nodes().
# Optionally, `algargs` is a list of arguments to pass to the
# modularity algorithm.
def consensus_matrix(g, algorithm, it=500, algargs={}):
    assert type(g) == networkx.classes.graph.Graph, "Not a graph"
    assert callable(algorithm), "f is not a function"
    assert type(it) == int, "Non-integer iterations"
    assert it > 0, "It needs to be greater than 0"
    consensus = numpy.zeros((len(g), len(g)))
    for i in range(0, it):
        p = algorithm(g, **algargs)
        assert type(p) == dict, "Wrong algorithm return type"
        assert len(p) == len(g), "Wrong algorithm return length"
        assert set(list(p.keys())) == set(g.nodes()), "Keys not nodes"
        consensus += numpy.array([[p[k] == p[j] for k in g.nodes()] for j in g.nodes()])
    consensus /= it
    return consensus

# consensus_matrix_par is the same as consensus_matrix, except it is
# parallelized.  It has the extra argument `processors`, which is the
# number of processors on which to run the parallelization.
def consensus_matrix_par(g, algorithm, it=500, processors=2, algargs={}):
    import multiprocessing, functools
    assert type(g) == networkx.classes.graph.Graph, "Not a graph"
    assert callable(algorithm), "f is not a function"
    assert type(it) == int, "Non-integer iterations"
    assert it > 0, "It needs to be greater than 0"
    assert type(processors) == int and processors > 0, "Invalid number of processors"
    # The idea here is that, since a consensus matrix is just a bunch
    # of matrices averaged together, we do "`it` divided by
    # `processors`" iterations on each processor and then average
    # together the results.  It sometimes does slightly more
    # iterations if it doesn't divide perfectly.
    # 
    # To do this, we make `processors` copies of the graph, and then
    # use a parallel pool map to send it to `processors` processors.
    itadj = int(numpy.ceil(it/processors))
    gs = [g for i in range(0, processors)] # This doesn't actually use any more memory since g is a reference
    # Use try-except catch-all to make sure we close the processes.
    # That way when we ctrl+c, we don't have 20 python processes
    # permanently running on our computer.
    try:
        p = multiprocessing.Pool(processors)
        cs = p.map(functools.partial(consensus_matrix, algorithm=algorithm, it=itadj, algargs=algargs), gs)
    finally:
        p.terminate()
    return numpy.mean(cs, axis=0)


# ▞▀▖                 ▗▐         ▜          ▗▐  ▌         
# ▌  ▞▀▖▛▚▀▖▛▚▀▖▌ ▌▛▀▖▄▜▀ ▌ ▌ ▝▀▖▐ ▞▀▌▞▀▖▙▀▖▄▜▀ ▛▀▖▛▚▀▖▞▀▘
# ▌ ▖▌ ▌▌▐ ▌▌▐ ▌▌ ▌▌ ▌▐▐ ▖▚▄▌ ▞▀▌▐ ▚▄▌▌ ▌▌  ▐▐ ▖▌ ▌▌▐ ▌▝▀▖
# ▝▀ ▝▀ ▘▝ ▘▘▝ ▘▝▀▘▘ ▘▀▘▀ ▗▄▘ ▝▀▘ ▘▗▄▘▝▀ ▘  ▀▘▀ ▘ ▘▘▝ ▘▀▀ 


# find_communities_louvain uses the louvain algorithm to detect
# community (module) structure.  The input should be a graph in
# networkx graph format.  This function adds the "louvain" property to
# each of the nodes, describing which module (community) the node
# belongs to.  It returns a dictionary, where each index is a node and
# the value is the value is an integer representing the community
# index that node is a part of.
def find_communities_louvain(g, gamma=1):
    assert type(g) == networkx.classes.graph.Graph, "Not a graph"
    #assert networkx.is_connected(g), "Graph not connected"
    lmodule = bct.community_louvain(numpy.asarray(networkx.to_numpy_matrix(g)), gamma=gamma)
    c = dict(zip(g.nodes(), list(map(int, lmodule[0]))))
    networkx.set_node_attributes(g, 'louvain', c)
    return c


