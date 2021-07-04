#!/usr/bin/env python

import json
import logging

import numpy
import numpy as np

##

logging.basicConfig()
logger=logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


## Create a test graph network
##
## see the graph here: https://en.wikipedia.org/wiki/File:Graph_with_Chordless_and_Chorded_Cycles.svg
##
if False:
    import networkx as nx

    #
    # { "node" : {"connected nodes"},
    #   ...
    # }
    dict_of_dicts = {
      "A" : {"B", "F"},
      "B" : {"A", "C", "G"},
      "C" : {"B", "D", "G", "L"},
      "D" : {"C", "E", "K"},
      "E" : {"D", "F"},
      "F" : {"A", "E"},
      "G" : {"B", "C", "H", "L"},
      "H" : {"G", "I"},
      "I" : {"H", "J", "K"},
      "J" : {"I", "K"},
      "K" : {"D", "I", "J", "L"},
      "L" : {"C", "G", "K"},
    }

    G = nx.Graph(dict_of_dicts)
    G.edges["A","B"]["color"] = "blue"
    print(list(G.nodes))
    print(list(G.edges))


    # make it a directed bipartite graph
    ## generate reactions
    graph_data_dict = {}
    reaction_numbers = []
    for reaction_number, edge in enumerate(G.edges):
        # because every edge is expanded and we can only construct in-coming edges,
        # split the data:
        # A - B
        # becomes:
        # A -> 1 ; 1 -> B
        assert len(edge) == 2
        first_node  = edge[0]
        second_node = edge[1]
        reaction_numbers.append(reaction_number)

        #print(first_node, reaction_number, second_node)

        ## add outgoing links if not existent, else create
        if first_node in graph_data_dict:
            graph_data_dict[first_node].append(reaction_number)
        else:
            graph_data_dict[first_node] = [ reaction_number ]

        # reaction_number is unique, can be created always safely
        graph_data_dict.update( {  reaction_number : [ second_node ] } )

        #print(graph_data_dict)

    B = nx.DiGraph(graph_data_dict)
    # label reactions
    for reaction_number in reaction_numbers:
        B.nodes[reaction_number]["bipartite"] = 1
        B.nodes[reaction_number]["type"]      = "reaction"
    # label compounds
    for node_key in B.nodes:
        if B.nodes[node_key] == {}:
            B.nodes[node_key]["bipartite"] = 0
            B.nodes[node_key]["type"]      = "compound"

    ## generate stoichiometry matrix
    reaction_numbers = { reaction_number for reaction_number, data in B.nodes(data=True) if data["type"] == "reaction" }
    compound_numbers = { compound_number for compound_number, data in B.nodes(data=True) if data["type"] == "compound" }


    ## separate nodes
    #top_nodes = {n for n, d in B.nodes(data=True) if d["bipartite"] == 0}
    #bottom_nodes = set(B) - top_nodes

    print(nx.neighbors(B, 5))
    #for r in reaction_numbers:
    #    print(B["D"])

    ## plot bipartite
    import matplotlib.pyplot as plt
    ax = plt.subplot(121)
    pos = nx.multipartite_layout(B, subset_key="bipartite")
    nx.draw(B, pos=pos, ax=ax, with_labels=True)
    ax = plt.subplot(222)
    nx.draw_kamada_kawai(G, ax=ax, with_labels=True)
    ax = plt.subplot(224)
    nx.draw_kamada_kawai(B, ax=ax, with_labels=True)
    plt.show()


## template:
if False:
    ## figure x) name / description of test
    ##
    ## indicate reactions:
    #  rxn0 : A = B
    #  rxn1 : B = C
    #  rxn...
    ##
    ## indicate designed cycles:
    #  cycle I:  0 = rxn0 + rxn1 - rxn2  == A-B-C-A
    #  cycle II: ...
    ##
    ## indicate known solutions describing the whole figure:
    #  solution s1: cycle I + cycle II
    ##
    #  comments if needed
    ##
    #  number of nodes (mols) in the figure:     x
    #  number of edges (rxns) in the figure:     y
    #  linear dependent reactions of the system: z
    ##
    stoich_matrix = np.array(
    [                                         ## mol ...
    # rxn ... --->                            ##     |
    # rxn 0   1   2   3   4   5   6   7   8   ##
        [-1, 00, 00, 00, 00, 00, 00, 00, 00], ## mol A
        [+1, 00, 00, 00, 00, 00, 00, 00, 00], ##     B
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##     C
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##     D
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##     E
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##     F
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##     G
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##     H
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##     I
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##     J
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##     K
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##     L
    ])

## one-substrate one-product scenarios:

if False:
    ## figure a) minimal trivial cycle of length 3
    ##
    #  rxn0 : A = B
    #  rxn1 : A = C
    #  rxn2 : B = C
    ##
    #  cycle I: 0 = rxn0 + rxn2 - rxn1  == A-B-C-A
    ##
    #  number of nodes (mols) in the figure:     3
    #  number of edges (rxns) in the figure:     3
    #  linear dependent reactions of the system: 1
    ##
    stoich_matrix = np.array(
    [                 ## mol ...
    # rxn ... --->    ##     |
    # rxn 0   1   2   ##     |
        [-1, -1 ,00], ## mol A
        [+1, 00, -1], ##     B
        [00, +1, +1], ##     C
    ])
    ##
    # rref = np.array(
    # [[1, 0, -1],
    #  [0, 1,  1],
    #  [0, 0,  0]])



if False:
    ## figure b) minimal trivial cycle of length 3 + linear path away
    ##
    #  rxn0 : A = B
    #  rxn1 : A = C
    #  rxn2 : B = C
    #  rxn3 : C = D
    #  rxn4 : D = E
    ##
    #  cycle I: 0 = rxn0 + rxn2 - rxn1  == A-B-C-A
    ##
    #  number of nodes (mols) in the figure:     5
    #  number of edges (rxns) in the figure:     5
    #  linear dependent reactions of the system: 1
    ##
    stoich_matrix = np.array(
    [                         ## mol ...
    # rxn ... --->            ##     |
    # rxn 0   1   2,  3   4   ##     |
        [-1, -1 ,00, 00, 00], ## mol A
        [+1, 00, -1, 00, 00], ##     B
        [00, +1, +1, -1, 00], ##     C
        [00, 00, 00, +1, -1], ##     D
        [00, 00, 00, 00, +1], ##     E
    ])
    ##
    # rref = np.array(
    # [[1, 0, -1, 0, 0],
    #  [0, 1,  1, 0, 0],
    #  [0, 0,  0, 1, 0],
    #  [0, 0,  0, 0, 1],
    #  [0, 0,  0, 0, 0]])


if False:
    ## figure c) two minimal trivial cycles of length 3 without shared edge
    ##
    #  rxn0 : A = B
    #  rxn1 : A = C
    #  rxn2 : B = C
    #  rxn3 : D = E
    #  rxn4 : D = F
    #  rxn5 : E = F
    ##
    #  cycle I:  0 = rxn0 + rxn2 - rxn1   == A-B-C-A
    #  cycle II: 0 = rxn3 + rxn5 - rxn4   == D-E-F-D
    ##
    #  solution s1: cycle I + cycle II
    ##
    #  number of nodes (mols) in the figure:     6
    #  number of edges (rxns) in the figure:     6
    #  linear dependent reactions of the system: 2
    ##
    stoich_matrix = np.array(
    [                             ## mol ...
    # rxn ... --->                ##     |
    # rxn 0   1   2,  3   4,  5   ##     |
        [-1, -1 ,00, 00, 00, 00], ## mol A
        [+1, 00, -1, 00, 00, 00], ##     B
        [00, +1, +1, 00, 00, 00], ##     C
        [00, 00, 00, -1, -1, 00], ##     D
        [00, 00, 00, +1, 00, -1], ##     E
        [00, 00, 00, 00, +1, +1], ##     F
    ])
    ##
    # rref = np.array(
    # [[1, 0, -1, 0, 0,  0],
    #  [0, 1,  1, 0, 0,  0],
    #  [0, 0,  0, 1, 0, -1],
    #  [0, 0,  0, 0, 1,  1],
    #  [0, 0,  0, 0, 0,  0],
    #  [0, 0,  0, 0, 0,  0]])

if False:
    ## figure d) two minimal trivial cycles of length 3 with one shared edge
    ##
    # rxn0 : A = B
    # rxn1 : A = C
    # rxn2 : B = C
    # rxn3 : B = D
    # rxn4 : C = D
    ##
    # cycle I:  0 = rxn0 + rxn2 - rxn1   == A-B-C-A
    # cycle II: 0 = rxn2 + rxn4 - rxn3   == B-C-D-B
    ##
    # solution s1: cycle I + cycle II
    ##
    # the solution s1 involves an edge that appears twice, in each cycle.
    # There are more complex solutions, but those are clearly inferior, e.g.:
    # cycle III: 0 = rxn0 + rxn3 - rxn4 - rxn1  == A-B-C-D-A
    # solution s2: cycle I + cycle III
    # This solution would have two edges shared between its cycles (whereas
    # solution s1 has only one edge shared between its two cycles), and
    # cycle III has length 4 instead of the cycles of length 3 from solution s1.
    #
    #  |--------------------------------------------------------------------------------------------------------|
    #  |            | number of |   length of   | total edges | total nodes |     edges in    |    nodes in     | 
    #  | solution # |   cycles  | longest cycle | in solution | in solution | multiple cycles | multiple cycles | 
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #  |     s1     |     2     |      3        |      6      |      6      |      1 in 2     |     2 in 2      |
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #  |     s2     |     2     |      4        |      7      |      7      |      2 in 2     |     3 in 2      |
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #
    ##
    #  number of nodes (mols) in the figure:     6
    #  number of edges (rxns) in the figure:     5
    #  linear dependent reactions of the system: 2
    ##
    stoich_matrix = np.array(
    [                         ## mol ...
    # rxn ... --->            ##     |
    # rxn 0   1   2,  3   4   ##     |
        [-1, -1 ,00, 00, 00], ## mol A
        [+1, 00, -1, -1, 00], ##     B
        [00, +1, +1, 00, -1], ##     C
        [00, 00, 00, +1, +1], ##     D
    ])
    ##
    # rref = np.array(
    # [[1, 0, -1, 0,  1],
    #  [0, 1,  1, 0, -1],
    #  [0, 0,  0, 1,  1],
    #  [0, 0,  0, 0,  0]])

if False:
    ## figure e) multiple distinct equally-well solutions: two cycles of length 4 each
    ##
    #  rxn0 : A = B
    #  rxn1 : A = D
    #  rxn2 : B = C
    #  rxn3 : B = E
    #  rxn4 : C = D
    #  rxn5 : D = E
    ##
    #  cycle I:   0 = rxn0 + rxn2 + rxn4 - rxn1  == A-B-C-D-A
    #  cycle II:  0 = rxn0 + rxn3 - rxn5 - rxn1  == A-B-E-D-A
    #  cycle III: 0 = rxn2 + rxn4 + rxn5 - rxn3  == B-C-D-E-B
    ##
    #  solution s1: cycle I  + cycle II
    #  solution s2: cycle I  + cycle III
    #  solution s3: cycle II + cycle III
    ##
    #  When starting with the path D-A-B, there is no difference
    #  in whether to choose C or E as the "bridge" to close the cycle.
    #  One can then describe the figure by chosing two cycles, each with a
    #  different bridge (solution s1: cycle I + cycle II).
    #  But: one can also construct cycle III.
    #  Is it equally well to use solution s2: cycle I + cycle III, or
    #  solution s3: cycle II + cycle III?
    #  All solutions have two edges which are present in two
    #  cycles; there is no "measurable" difference between them.
    #
    #  |--------------------------------------------------------------------------------------------------------|
    #  |            | number of |   length of   | total edges | total nodes |     edges in    |    nodes in     | 
    #  | solution # |   cycles  | longest cycle | in solution | in solution | multiple cycles | multiple cycles | 
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #  |     s1     |     2     |      4        |      8      |      8      |      2 in 2     |     3 in 2      |
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #  |     s2     |     2     |      4        |      8      |      8      |      2 in 2     |     3 in 2      |
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #  |     s3     |     2     |      4        |      8      |      8      |      2 in 2     |     3 in 2      |
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #
    ##
    #  number of nodes (mols) in the figure:     5
    #  number of edges (rxns) in the figure:     6
    #  linear dependent reactions of the system: 2
    ##
    stoich_matrix = np.array(
    [                             ## mol ...
    # rxn ... --->                ##     |
    # rxn 0   1   2   3   4   5   ##
        [-1, -1, 00, 00, 00, 00], ## mol A
        [+1, 00, -1, -1, 00, 00], ##     B
        [00, 00, +1, 00, -1, 00], ##     C
        [00, +1, 00, 00, +1, -1], ##     D
        [00, 00, 00, +1, 00, +1], ##     E
    ])
    ##
    # rref = np.array(
    # [[1, 0, 0, 0, -1,  1],
    #  [0, 1, 0, 0,  1, -1],
    #  [0, 0, 1, 0, -1,  0],
    #  [0, 0, 0, 1,  0,  1],
    #  [0, 0, 0, 0,  0,  0]]


if False:
    ## figure f) multiple distinct equally-well solutions: three cycles of length 3 each
    ##
    #  rxn0 : A = B
    #  rxn1 : A = C
    #  rxn2 : A = D
    #  rxn3 : B = C
    #  rxn4 : B = D
    #  rxn5 : C = D
    ##
    #  cycle I:   0 = rxn0 + rxn3 - rxn1    == A-B-C-A
    #  cycle II:  0 = rxn3 + rxn5 - rxn4    == B-C-D-B
    #  cycle III: 0 = rxn2 - rxn5 - rxn1    == A-D-C-A
    #  cycle IV:  0 = rxn0 + rxn4 - rxn2    == A-B-D-A
    ##
    #  solution s1: cycle I  + cycle II  + cycle III
    #  solution s2: cycle I  + cycle II  + cycle IV
    #  solution s3: cycle I  + cycle III + cycle IV
    #  solution s4: cycle II + cycle III + cycle IV
    ##
    #  After constructing, e.g. cycle I and cycle II, is it equally well to
    #  construct cycle III or cycle IV? Both include the final missing edge (rxn2),
    #  but they choose different nodes / edges, compared to
    #  the already existing cycles. Illustrating the problem explicitly:
    #  If having chosen cycle I and cycle II already,
    #  and chosing now: 1) cycle III, the node C is present in three cycles,
    #  and the edges rxn1, rxn3, and rxn5 are present in two cycles each;
    #  whereas when chosing now: 2) cycle IV, the node B is present in
    #  three cycles, and the edges rxn0, rxn3, and rxn4
    #  are present in two cycles each.
    #  The solutions describe the figure equally well, but are distinct.
    #
    #  |--------------------------------------------------------------------------------------------------------|
    #  |            | number of |   length of   | total edges | total nodes |     edges in    |    nodes in     | 
    #  | solution # |   cycles  | longest cycle | in solution | in solution | multiple cycles | multiple cycles | 
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #  |     s1     |     3     |      3        |      9      |      9      |      3 in 2     | 1 in 3, 3 in 2  |
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #  |     s2     |     3     |      3        |      9      |      9      |      3 in 2     | 1 in 3, 3 in 2  |
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #  |     s3     |     3     |      3        |      9      |      9      |      3 in 2     | 1 in 3, 3 in 2  |
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #  |     s4     |     3     |      3        |      9      |      9      |      3 in 2     | 1 in 3, 3 in 2  |
    #  |--------------------------------------------------------------------------------------------------------|
    #
    ##
    #  number of nodes (mols) in the figure:     4
    #  number of edges (rxns) in the figure:     6
    #  linear dependent reactions of the system: 3
    ##
    stoich_matrix = np.array(
    [                             ## mol ...
    # rxn ... --->                ##     |
    # rxn 0   1   2   3   4   5   ##
        [-1, -1, -1, 00, 00, 00], ## mol A
        [+1, 00, 00, -1, -1, 00], ##     B
        [00, +1, 00, +1, 00, -1], ##     C
        [00, 00, +1, 00, +1, +1], ##     D
    ])
    ##
    # rref = np.array(
    # [[1, 0, 0, -1, -1,  0],
    #  [0, 1, 0,  1,  0, -1],
    #  [0, 0, 1,  0,  1,  1],
    #  [0, 0, 0,  0,  0,  0]])


if False:
    ## figure g) trade-off between objectives: two cycles of length 4, 3; or three cycles of length 3, 3, 3?
    ##
    #  rxn0 : A = B
    #  rxn1 : A = C
    #  rxn2 : A = D
    #  rxn3 : A = E
    #  rxn4 : B = C
    #  rxn5 : C = D
    #  rxn6 : C = E
    ##
    #  cycle I:   0 = rxn0 + rxn4 + rxn5 - rxn2     == A-B-C-D-A
    #  cycle II:  0 = rxn1 + rxn6 - rxn3            == A-C-E-A
    #  cycle III: 0 = rxn0 + rxn4 - rxn1            == A-B-C-A
    #  cycle IV:  0 = rxn1 + rxn5 - rxn2            == A-C-D-A
    #  cycle V:   0 = rxn0 + rxn4 - rxn6 - rxn3     == A-B-C-E-A
    #  cycle VI:  0 = rxn2 - rxn5 - rxn4 - rxn0     == A-D-C-E-A
    ##
    #  solution s1: cycle I   + cycle II
    #  solution s2: cycle II  + cycle III + cycle IV
    #  solution s3: cycle IV  + cycle V
    #  solution s4: cycle III + cycle VI
    ##
    #  A complete description of the figure can be achieved e.g. with
    #  solution s1, which is distinct from solution s2. Depending on the
    #  objectives one sets, one solution is better than the other, although
    #  both describe the figure completely.
    #  solution s1 contains only two cycles, but one has length 4;
    #  two nodes (A and C) are present in two cycles.
    #  solution s2 contains three cycles of length 3 each;
    #  two nodes (A and C) are present in all three cycles; also
    #  one edge (rxn1) is present in two cycles, whereas solution s1 does
    #  not contain any edge present in multiple cycles.
    #  solution s3 and solution s4 are similar to solution s1, but this is
    #  a "multiple distinct equally-well solutions" problem not to be
    #  considered here. The question really is: Which solution is "better",
    #  solution s1 or solution s2?
    #
    #  |--------------------------------------------------------------------------------------------------------|
    #  |            | number of |   length of   | total edges | total nodes |     edges in    |    nodes in     | 
    #  | solution # |   cycles  | longest cycle | in solution | in solution | multiple cycles | multiple cycles | 
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #  |     s1     |     2     |      4        |      7      |      7      |        0        |     2 in 2      |
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #  |     s2     |     3     |      3        |      9      |      9      |      1 in 3     |     2 in 3      |
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #  |     s3     |     2     |      4        |      7      |      7      |        0        |     2 in 2      |
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #  |     s4     |     2     |      4        |      7      |      7      |        0        |     2 in 2      |
    #  |--------------------------------------------------------------------------------------------------------|
    #
    ##
    #  number of nodes (mols) in the figure:     5
    #  number of edges (rxns) in the figure:     7
    #  linear dependent reactions of the system: 3
    ##
    stoich_matrix = np.array(
    [                                 ## mol ...
    # rxn ... --->                    ##     |
    # rxn 0   1   2   3   4   5   6   ##
        [-1, -1, -1, -1, 00, 00, 00], ## mol A
        [+1, 00, 00, 00, -1, 00, 00], ##     B
        [00, +1, 00, 00, +1, -1, -1], ##     C
        [00, 00, +1, 00, 00, +1, 00], ##     D
        [00, 00, 00, +1, 00, 00, +1], ##     E
    ])
    ##
    # rref = np.array(
    # [[1, 0, 0, 0, -1,  0,  0],
    #  [0, 1, 0, 0,  1, -1, -1],
    #  [0, 0, 1, 0,  0,  1,  0],
    #  [0, 0, 0, 1,  0,  0,  1],
    #  [0, 0, 0, 0,  0,  0,  0]])


if False:
    ## figure h) trade-off between objectives: two cycles of length 5, 4; or three cycles of length 4, 4, 3?
    ##
    #  rxn0 : A = B
    #  rxn1 : A = F
    #  rxn2 : B = C
    #  rxn3 : B = E
    #  rxn4 : B = G
    #  rxn5 : C = D
    #  rxn6 : D = E
    #  rxn7 : E = F
    #  rxn8 : E = G
    ##
    #  cycle I:   0 = rxn0 + rxn4 - rxn8 + rxn7 - rxn1          == A-B-G-E-F-A
    #  cycle II:  0 = rxn2 + rxn5 + rxn6 - rxn3                 == B-C-D-E-B
    #  cycle III: 0 = rxn0 + rxn3 + rxn7 - rxn1                 == A-B-E-F-A
    #  cycle IV:  0 = rxn4 - rxn8 - rxn3                        == B-G-E-B
    #  cycle V:   0 = rxn0 + rxn2 + rxn5 + rxn6 + rxn7 - rxn1   == A-B-C-D-E-F-A
    ##
    #  solution s1: cycle I  + cycle II
    #  solution s2: cycle II + cycle III + cycle IV
    #  solution s3: cycle IV + cycle V
    ##
    #  Here, solution s1 has actually a bigger length of longest cycle compared to
    #  solution s2 (5 vs 4), still it has a smaller number of cycles and even a smaller
    #  number of total edges (and nodes) in the solution.
    #  Solution s3 is strange in that most numbers are optimal except for the length
    #  of the longest cycle.
    #  Note also that cycle I is the XOR combination of cycle III and cycle IV, and
    #  thus might not be expected to be involved in an optimal solution.
    #
    #  |--------------------------------------------------------------------------------------------------------|
    #  |            | number of |   length of   | total edges | total nodes |     edges in    |    nodes in     | 
    #  | solution # |   cycles  | longest cycle | in solution | in solution | multiple cycles | multiple cycles | 
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #  |     s1     |     2     |      5        |      9      |      9      |        0        |     2 in 2      |
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #  |     s2     |     3     |      4        |     11      |     11      |      1 in 3     |     2 in 3      |
    #  |------------|-----------|---------------|-------------|-------------|-----------------|-----------------|
    #  |     s3     |     2     |      6        |      9      |      9      |        0        |     2 in 2      |
    #  |--------------------------------------------------------------------------------------------------------|
    #
    ##
    #  number of nodes (mols) in the figure:     7
    #  number of edges (rxns) in the figure:     9
    #  linear dependent reactions of the system: 3
    ##
    stoich_matrix = np.array(
    [                                         ## mol ...
    # rxn ... --->                            ##     |
    # rxn 0   1   2   3   4   5   6   7   8   ##
        [-1, -1, 00, 00, 00, 00, 00, 00, 00], ## mol A
        [+1, 00, -1, -1, -1, 00, 00, 00, 00], ##     B
        [00, 00, +1, 00, 00, -1, 00, 00, 00], ##     C
        [00, 00, 00, 00, 00, +1, -1, 00, 00], ##     D
        [00, 00, 00, +1, 00, 00, +1, -1, -1], ##     E
        [00, +1, 00, 00, 00, 00, 00, +1, 00], ##     F
        [00, 00, 00, 00, +1, 00, 00, 00, +1], ##     G
    ])
    ##
    # rref = np.array(
    # [[1, 0, 0, 0, 0, 0,  0, -1,  0],
    #  [0, 1, 0, 0, 0, 0,  0,  1,  0],
    #  [0, 0, 1, 0, 0, 0, -1,  0,  0],
    #  [0, 0, 0, 1, 0, 0,  1, -1, -1],
    #  [0, 0, 0, 0, 1, 0,  0,  0,  1],
    #  [0, 0, 0, 0, 0, 1, -1,  0,  0],
    #  [0, 0, 0, 0, 0, 0,  0,  0,  0]])


if False:
    ## figure i) why are redundant cycles needed?
    ##
    #  rxn0 : A = B
    #  rxn1 : A = C
    #  rxn2 : A = D
    #  rxn3 : A = E
    #  rxn4 : B = F
    #  rxn5 : C = F
    #  rxn6 : D = F
    #  rxn7 : E = F
    ##
    #  cycle I:   0 = rxn0 + rxn4 - rxn5 - rxn1  == A-B-F-C-A
    #  cycle II:  0 = rxn3 + rxn7 - rxn6 - rxn2  == A-E-F-D-A
    #  cycle III: 0 = rxn1 + rxn5 - rxn6 - rxn2  == A-C-F-D-A
    #  cycle IV:  0 = rxn0 + rxn4 - rxn7 - rxn3  == A-B-F-E-A
    #  cycle V:   0 = rxn0 + rxn4 - rxn6 - rxn2  == A-B-F-D-A
    #  cycle VI:  0 = rxn3 + rxn7 - rxn5 - rxn1  == A-E-F-C-A
    ##
    #  solution s1: cycle I  + cycle II  (NB: not really a solution)
    #  solution s2: cycle I  + cycle II + cycle III
    #  solution s3: cycle IV + cycle V  + cycle VI
    ##
    #  Why is solution s1 not sufficient? It covers the whole figure!
    #  But: from the cycles of solution s1 it is not possible to construct
    #  all other possible cycles with XOR operations, e.g. it is not possible
    #  to construct cycle IV.
    #  In contrast, with solution s2, which adds a seemingly redundant cycle (cycle III),
    #  it becomes now possible to construct the others: cycle IV, V, VI, ...
    ##
    #  number of nodes (mols) in the figure:     6
    #  number of edges (rxns) in the figure:     8
    #  linear dependent reactions of the system: 3
    ##
    stoich_matrix = np.array(
    [                                         ## mol ...
    # rxn ... --->                            ##     |
    # rxn 0   1   2   3   4   5   6   7   ##
        [-1, -1, -1, -1, 00, 00, 00, 00], ## mol A
        [+1, 00, 00, 00, -1, 00, 00, 00], ##     B
        [00, +1, 00, 00, 00, -1, 00, 00], ##     C
        [00, 00, +1, 00, 00, 00, -1, 00], ##     D
        [00, 00, 00, +1, 00, 00, 00, -1], ##     E
        [00, 00, 00, 00, +1, +1, +1, +1], ##     F
    ])
    # rref = np.array(
    # [[1, 0, 0, 0, 0,  1,  1,  1],
    #  [0, 1, 0, 0, 0, -1,  0,  0],
    #  [0, 0, 1, 0, 0,  0, -1,  0],
    #  [0, 0, 0, 1, 0,  0,  0, -1],
    #  [0, 0, 0, 0, 1,  1,  1,  1],
    #  [0, 0, 0, 0, 0,  0,  0,  0]])



## multi-substrate multi-product scenarios:

if False:
    ## mixed one/two-substrate/product cycle of length 3
    ##
    # rxn0 : A = B
    # rxn1 : C = D
    # rxn2 : A + C = B + D
    ##
    # 0 = rxn0 + rxn1 - rxn2
    ##
    stoich_matrix = np.array(
    [                 ## mol ...
    # rxn ... --->    ##     |
    # rxn 0   1   2   ##     |
        [+1, 00, +1], ## mol A
        [-1, 00, -1], ##     B
        [00, +1, +1], ##     C
        [00, -1, -1], ##     D
    ])

if False:
    ## minimal trivial two-substrate two-product cycle of length 3
    ##
    # rxn0 : A + B = C + D
    # rxn1 : C + D = E + F
    # rxn2 : A + B = E + F
    ##
    # 0 = rxn0 - rxn1 - rxn2
    ##
    stoich_matrix = np.array(
    [                 ## mol ...
    # rxn ... --->    ##     |
    # rxn 0   1   2   ##     |
        [+1, 00, +1], ## mol A
        [+1, 00, +1], ##     B
        [-1, +1, 00], ##     C
        [-1, +1, 00], ##     D
        [00, -1, -1], ##     E
        [00, -1, -1], ##     F
    ])



if False:
    ## what is this?
    stoich_matrix = np.array(
    [                                         ## mol ...
    # rxn ... --->                            ##     |
    # rxn 0   1   2   3   4   5   6   7   8   ##
        [+1, +1 ,00 ,00, 00, 00, 00, -1, 00], ## mol A
        [+1, 00, +2, 00, 00, 00, 00, -1, 00], ##     B
        [-1, -1, 00, 00, 00, 00, 00, 00, 00], ##     C
        [-1, 00, -2, 00, 00, 00, 00, 00, 00], ##     D
        [00, +1, -2, 00, 00, 00, 00, 00, 00], ##     E
        [00, -1, +2, 00, 00, 00, 00, 00, 00], ##     F
        [00, 00, 00, -1, 00, -1, -1, +1, 00], ##     G
        [00, 00, 00, +1, -1, 00, 00, 00, 00], ##     H
        [00, 00, 00, +1, -1, 00, 00, 00, 00], ##     I
        [00, 00, 00, 00, +1, +1, 00, 00, 00], ##     J
        [00, 00, 00, 00, 00, 00, +1, 00, 00], ##     K
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##     L
    ])


collected_mols = list( [f"mol_{chr(65+i)}" for i in range(stoich_matrix.shape[0])])
collected_rxns = list( [f"rxn_{i}" for i in range(stoich_matrix.shape[1])])


# "prove" rref equivalence for different orderings of the stoich. matrix
if False:

    import itertools
    import sympy

    number_of_rows    = stoich_matrix.shape[0]
    number_of_columns = stoich_matrix.shape[1]
    original_rref, original_inds = sympy.Matrix(stoich_matrix).rref()
    number_of_dependent_columns = number_of_columns - len(original_inds)
    logger.debug(original_rref)
    logger.debug(original_inds)
    logger.debug(number_of_dependent_columns)

    ## the reduced row echolon form (rref) is NOT the same for permutated columns (i.e. rxn ordering)
    for column_ordering in itertools.permutations(range(number_of_columns)):
        #logger.debug(f"order of columns to investigate: {column_ordering}")
        reordered_stoich_matrix = stoich_matrix[:,column_ordering]
        _, inds = sympy.Matrix(reordered_stoich_matrix).rref()
        #logger.debug(inds)
        #logger.debug(_)
        dependent_columns = sorted(list([column_ordering[i] for i in range(number_of_columns) if i not in inds]))
        logger.debug(dependent_columns)

        assert len(dependent_columns) == number_of_dependent_columns


    ## the reduced row echolon form (rref) is the same for permutated rows (i.e. mol ordering)
    for row_ordering in itertools.permutations(range(number_of_rows)):
        break
        logger.debug(f"order of rows to investigate: {row_ordering}")
        reordered_stoich_matrix = stoich_matrix[row_ordering,:]
        _, inds = sympy.Matrix(reordered_stoich_matrix).rref()
        assert original_rref == _
        assert original_inds == inds



## calculate rref
if True:
    import time
    import sympy

    ## time the rref calculation for subsets
    if False:
        stoich_matrix2 = stoich_matrix.copy()
        for i in []: #10]: #,100,300,500]:
            stoich_matrix = stoich_matrix2[:i]
            print(stoich_matrix.shape)

            t1 = time.time()
            print("calculating rref...")
            _, inds = sympy.Matrix(stoich_matrix).rref()
            print("...done.")
            t2 = time.time()
            print(f"t: {t2-t1}")

        # reset the stoich_matrix
        stoich_matrix = stoich_matrix2

    ## time the actual calculation
    t1 = time.time()
    print("calculating rref...")
    _, inds = sympy.Matrix(stoich_matrix).rref()
    print("...done.")
    t2 = time.time()
    print(f"t: {t2-t1}")

    print("rref result:")
    print(_)
    print(inds)

    print("saving rref's indxs..")
    with open("./inds.npz", "wb") as f:
        numpy.savez(f, inds=inds)

else:
    print("loading rref's indxs..")
    with open("./inds.npz", "rb") as f:
        loaded = numpy.load(f)
        inds = loaded["inds"]
print("...done.")

## this gave us all linearly independent columns, but we are interested in the others:
linear_dependent_rxn_indices = [x for x in range(len(collected_rxns)) if not x in inds]
linear_dependent_rxns = [ collected_rxns[i] for i in linear_dependent_rxn_indices ]
print(f"The number of linear dependent reactions is: {len(linear_dependent_rxns)}")


if False:
    ### reduce matrix so it contains linearly independent columns only
    matrix_complete = stoich_matrix.copy()
    matrix_reduced  = stoich_matrix.copy() #[:, [row_idx for row_idx in range(len(collected_rxns)) if not row_idx in linear_dependent_rxn_indices]]
    collected_rxns_reduced = [r for r in collected_rxns] # if not r in linear_dependent_rxns]
    print(f"collected_rxns_reduced : {collected_rxns_reduced}")
    print("---")

    ## alternative to integer Linear Programming: least squares solution:
    # try:
    #     x, residuals, ranks, s = numpy.linalg.lstsq(matrix_reduced, matrix_complete[:, solve_for_rxn_id])
    # except numpy.LinAlgError as e:
    #     x, residuals, ranks, s = None, None, None, None
    #     print(e)



    ## if desired, only do a calculation for a specific index of the dependent reactions
    import json
    import pulp
    ## these variables allow to spread jobs
    linear_dependent_rxn_indices_to_calculate = [x for x in range(len(collected_rxns))] # linear_dependent_rxn_indices
    _filename_extra = ""
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("number_of_job", type=int, nargs='?')
    args = parser.parse_args()
    print(f"job number: {args.number_of_job}")
    if args.number_of_job or args.number_of_job == 0:
        linear_dependent_rxn_indices_to_calculate = [ linear_dependent_rxn_indices[args.number_of_job], ]
        _filename_extra = f"_job{args.number_of_job}"


    output_collector={}
    for solve_for_rxn_id in linear_dependent_rxn_indices_to_calculate:
        print(collected_rxns[solve_for_rxn_id])

        matrix_reduced  = stoich_matrix.copy()[:,[c for c in range(stoich_matrix.shape[1]) if not c == solve_for_rxn_id]]
        #print(matrix_reduced)
        print("setting up the problem...")
        rxn_variables  = list([pulp.LpVariable(f"rxn_{rxn_id}", -3, 3, pulp.const.LpInteger) for rxn_id in collected_rxns_reduced if not rxn_id == collected_rxns[solve_for_rxn_id]])
        prob = pulp.LpProblem(f"solve_for_rxn_{collected_rxns[solve_for_rxn_id]}", pulp.const.LpMinimize)
        for row_idx in range(matrix_reduced.shape[0]):
            constraint_for_row = 0
            for column_idx, value in enumerate(matrix_reduced[row_idx]):
                constraint_for_row += value * rxn_variables[column_idx]
            prob += constraint_for_row == matrix_complete[row_idx, solve_for_rxn_id], f"row_{row_idx}"

        ## objective function is actually supposed to be sum(abs(rxn_vars)), but needs workaround
        abs_of_rxn_variables  = list([pulp.LpVariable(f"abs_of_rxn_{rxn_id}") for rxn_id in collected_rxns_reduced]) #, -3, 3, pulp.const.LpInteger
        objective_function = pulp.lpSum( [abs_of_rxn_var for abs_of_rxn_var in abs_of_rxn_variables] ) # if not abs_of_rxn_var is abs_of_rxn_variables[solve_for_rxn_id]] )
        ## workaround formula for abs, by replacing actual variable by dummy and defining constraints of actual<=>dummy:
        for abs_of_rxn_var, rxn_var in zip(abs_of_rxn_variables, rxn_variables):
            prob += rxn_var  <= abs_of_rxn_var
            prob += -rxn_var <= abs_of_rxn_var
        prob += objective_function
        print("...done.")

        #print(prob)

        print("Solving the problem...")
        #prob.solve()
        prob.solve(pulp.apis.PULP_CBC_CMD(msg=False))
        print("...done.")

        print(collected_rxns[solve_for_rxn_id])
        abs_objective = 0
        for r in rxn_variables:
            if r.value() == 0: continue
            print(f"{r} : {r.value()}")
            abs_objective += abs(r.value())
        print(abs_objective)

        ## check the result, just to be sure
        matrix = matrix_reduced.copy()
        for colidx, rxn in zip(range(len(rxn_variables)),rxn_variables):
            if rxn.value() is None:
                matrix[:, colidx] = 0
                continue
            else:
                matrix[:, colidx] = rxn.value()*matrix[:, colidx]
        resulting_overall_reaction = matrix.sum(axis=1)
        expected_overall_reaction  = matrix_complete[:,solve_for_rxn_id]
        test = expected_overall_reaction == resulting_overall_reaction
        assert test.all(), "Obtained result is actually not valid!"
        #print(test.all())

        output_dict = { collected_rxns[solve_for_rxn_id]: [] }
        for r, rxn in zip(rxn_variables, collected_rxns_reduced):
            if r.value() == 0: continue
            output_dict[collected_rxns[solve_for_rxn_id]].append( {rxn : r.value()} )

        output_collector.update(output_dict)

        with open(f"./metabolic_thermodynamic_cycles{_filename_extra}.json", "w") as json_file:
            json.dump(output_collector, json_file, indent=2)
