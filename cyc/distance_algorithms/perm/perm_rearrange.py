#!/usr/bin/python3
from Elias.elias2005 import *
import sys

def to_perm(g1,g2):
    m = {}
    for i,a in enumerate(g2):
        m[a] = i+1
    perm = []
    for a in g1:
        perm.append(m[a])
    return perm

if __name__ == "__main__" :
    g1 = eval("[%s]" % sys.argv[1])
    g2 = eval("[%s]" % sys.argv[2])
    perm = to_perm(g1, g2)

    algo = elias2005(perm)
    dist = algo.sort()
    sys.stdout.write(str(dist) + "\n")
    sys.stdout.flush()
