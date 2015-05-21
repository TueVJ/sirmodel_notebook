import networkx as nx
import numpy as np
import pandas as pd
import random
import collections


def single_p_infection_history(
        G, p, tau=4, initial_node=None, immunize=set(),
        nodelist=None, iterations=1, samplesize=1):
    N = G.number_of_nodes()
    nodeset = set(G.nodes_iter())
    immunize = immunize & nodeset
    if N == 0:
        raise ValueError('empty graph')
    if initial_node is None:
        initial_node = np.random.choice(list(set(G.nodes()) - immunize), samplesize)
    elif not isinstance(initial_node, collections.Iterable):
        # Assume you've given me a single node
        initial_node = [initial_node]
    if p < 0 or p > 1:
        raise ValueError('p must be between 0 and 1')
    slist = []
    ilist = []
    rlist = []
    infectcount = {n: 0 for n in nodeset}
    n_iters = 0
    for ni in initial_node:
        for _ in xrange(iterations):
            infected = set([ni])
            susceptible = set(G.nodes_iter()) - immunize - infected
            immune = set()
            times = {n: 0 for n in G.nodes_iter()}
            lslist = [len(susceptible)]
            lilist = [len(infected)]
            lrlist = [0]
            for n in infected:
                times[n] = 4
            while len(infected) > 0:
                # MAIN LOOP
                # Spread infection
                to_infect = set()
                for n in infected:
                    for nn in (set(G[n]) & susceptible):
                        if random.random() < p:
                            to_infect.add(nn)
                # Remove nodes with time 0
                for n in list(infected):
                    if times[n] <= 0:
                        infected.remove(n)
                        immune.add(n)
                # Reduce times
                times = {k: max(v-1, 0) for k, v in times.iteritems()}
                # Add new infected
                infected.update(to_infect)
                susceptible.difference_update(to_infect)
                for n in to_infect:
                    times[n] = tau
                # Log populations
                lslist.append(len(susceptible) + len(immunize))
                lilist.append(len(infected))
                lrlist.append(len(immune))
            slist.append(lslist)
            ilist.append(lilist)
            rlist.append(lrlist)
            for n in immune:
                infectcount[n] += 1
            n_iters += 1
    infectcount = {k: float(v)/n_iters for k, v in infectcount.iteritems()}
    # BUILD DATAFRAME
    outdata = pd.DataFrame({
        'Healthy': pd.DataFrame(slist).T.fillna(method='pad').mean(axis=1),
        'Sick': pd.DataFrame(ilist).T.fillna(0).mean(axis=1),
        'Immune': pd.DataFrame(rlist).T.fillna(method='pad').mean(axis=1)})
    return outdata, infectcount
