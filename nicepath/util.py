import networkx as nx
import numpy as np
import os
from itertools import islice

def calculateDistance(_conserved_atom_ratio, _operator):
    lower_thresholds = {'sqrt': 0.001, '':0.01, 'exp': 0.1} # longest distance possible for sqrt: 31.6, none: 100, exp: 8103.1
    assert lower_thresholds.get(_operator), "No lower threshold defined for operator%s"%_operator
    if _conserved_atom_ratio < lower_thresholds[_operator]:
        distance = transform_operator(1.0 / lower_thresholds[_operator] , _operator)
    else:
        distance = transform_operator(1.0 / _conserved_atom_ratio, _operator)
        distance = round(distance, 2)

    return distance


def collectStatistics(Statistics, pathway, diff):
    pw_length = len(pathway)-1
    if not Statistics.get(pw_length):
        Statistics[pw_length] = diff
    else:
        Statistics[pw_length] += diff
    return Statistics


def k_shortest_paths(Graph, source, target, k, weight):
    return list(islice(nx.shortest_simple_paths(Graph, source, target, weight=weight), k))


def getBoolean(_string):
    if _string in ['False', 'false', '0']:
        return False
    elif _string in ['True', 'true', '1']:
        return True
    else:
        raise ValueError('ERROR: Invalid value', _string, 'for conversion to boolean.')

def transform_operator(x, _operator):
    if _operator == 'sqrt':
        return np.sqrt(x)
    elif _operator == 'exp':
        return np.exp(x)/np.e
    else:
        return x
