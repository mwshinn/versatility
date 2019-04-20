# Versatility

This package implements versatility [(Shinn et al.,
2017)](https://www.nature.com/articles/s41598-017-03394-5), which
describes how closely affiliated a node is with a network community
structure.  It is written in Python3, and can only be guaranteed to
work there.  (This MAY work in Python2 if you `import __future__` but
this is untested... see code for details.)

Install with:

    pip3 install versatility

Alternatively, clone the git repo and install with:

    python3 setup.py install

Dependencies:

- Python3
- [networkx](https://networkx.github.io/)
- [Scipy](https://www.scipy.org/) (including numpy and matplotlib)
- [bctpy](https://github.com/aestrivex/bctpy): The module "bct" is
  bctpy, a port of the Brain Connectivity Toolbox to Python.  The
  latest version supports Python3, and can be installed most easily
  with "`pip install bcpty`".

See function help for full documentation, but the most useful
functions are:

- `find_nodal_versatility` - Compute the versatility of each node in a
  graph using a specific community detection algorithm.
- `find_nodal_mean_versatility` - Compute the versatility of each node
  across a spectrum of community detection algorithm parameters (most
  notably the resolution parameter) and find the average.
- `find_optimal_gamma_curve` - Find the mean and standard error of
  versatility across a spectrum of resolution parameters and
  (optionally) plot the result.  This is most useful for finding the
  best resolution parameter, e.g. in [Figure 3c of the original
  paper](https://www.nature.com/articles/s41598-017-03394-5/figures/3).

Here is a quick example to get you started:

    import networkx
    from versatility import *
    G = networkx.karate_club_graph()
    find_nodal_mean_versatility(G, find_communities_louvain, processors=2)
    find_nodal_versatility(G, find_communities_louvain, algargs={"gamma" : 0.5})

If you use this code, please cite:

    Shinn, M., Romero-Garcia, R., Seidlitz, J., Vasa, F., Vertes, P.,
    Bullmore, E. (2017). Versatility of nodal affiliation to
    communities. Scientific Reports 7: 4273.
    doi:10.1038/s41598-017-03394-5

Copyright 2016-2019 Maxwell Shinn (maxwell.shinn@yale.edu)
Available under the GNU GPLv3.
