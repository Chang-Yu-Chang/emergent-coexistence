import igraph
from igraph import *

# Read all file from 
#read("../data/temp/motif1.txt").feedback_arc_set()


g = read("../data/temp/networks/C1R4.txt")
print(g)

g = read("../data/temp/networks/motif1.txt")
g = g.as_undirected()
print(g)

g.feedback_arc_set()
g.es.is_mutual()

cb = g.cohesive_blocks()
cb.cohesions()

read("../data/temp/networks/motif2.txt").feedback_arc_set()  
