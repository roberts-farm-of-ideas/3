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
    ## name / description of test
    ##
    ## indicate reactions:
    #  rxn0 : A = B
    #  rxn1 : B = C
    #  rxn...
    ##
    ## indicate designed cycles:
    #  0 = rxn0 + rxn1 - rxn2
    ##
    #  comments if needed
    ##
    stoich_matrix = np.array(
    [                                         ## mol_...
    # rxn_... --->                            ##     |
    # rxn_0   1   2   3   4   5   6   7   8   ##
        [+1, 00, 00, 00, 00, 00, 00, 00, 00], ## mol_A
        [-1, 00, 00, 00, 00, 00, 00, 00, 00], ##    _B
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##    _C
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##    _D
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##    _E
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##    _F
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##    _G
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##    _H
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##    _I
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##    _J
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##    _K
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##    _L
    ])

## one-substrate one-product scenarios:

if False:
    ## minimal trivial cycle of length 3
    ##
    #  rxn0 : A = B
    #  rxn1 : B = C
    #  rxn2 : A = C
    ##
    #  0 = rxn0 + rxn1 - rxn2
    ##
    #  the cycle is A-B-C-A
    #
    stoich_matrix = np.array(
    [                 ## mol_...
    # rxn_... --->    ##     |
    # rxn_0   1   2   ##     |
        [+1, 00 ,+1], ## mol_A
        [-1, +1, 00], ##    _B
        [00, -1, -1], ##    _C
    ])


if False:
    ## minimal trivial cycle of length 3 + linear path away
    ##
    #  rxn0 : A = B
    #  rxn1 : B = C
    #  rxn2 : A = C
    #  rxn3 : C = D
    #  rxn4 : D = E
    ##
    #  0 = rxn0 + rxn1 - rxn2
    ##
    #  the cycle is A-B-C-A
    #
    stoich_matrix = np.array(
    [                         ## mol_...
    # rxn_... --->            ##     |
    # rxn_0   1   2,  3   4   ##     |
        [+1, 00 ,+1, 00, 00], ## mol_A
        [-1, +1, 00, 00, 00], ##    _B
        [00, -1, -1, +1, 00], ##    _C
        [00, 00, 00, -1, +1], ##    _D
        [00, 00, 00, 00, -1], ##    _E
        [00, 00, 00, 00, 00], ##    _F
    ])


if False:
    ## two minimal trivial cycles of length 3 without shared edge
    ##
    #  rxn0 : A = B
    #  rxn1 : B = C
    #  rxn2 : A = C
    #  rxn3 : D = E
    #  rxn4 : E = F
    #  rxn5 : D = F
    ##
    #  0 = rxn0 + rxn1 - rxn2
    #  0 = rxn3 + rxn4 - rxn5
    ##
    #  the cycles are: A-B-C-A and D-E-F-D
    ##
    stoich_matrix = np.array(
    [                             ## mol_...
    # rxn_... --->                ##     |
    # rxn_0   1   2,  3   4,  5   ##     |
        [+1, 00 ,+1, 00, 00, 00], ## mol_A
        [-1, +1, 00, 00, 00, 00], ##    _B
        [00, -1, -1, 00, 00, 00], ##    _C
        [00, 00, 00, +1, 00, +1], ##    _D
        [00, 00, 00, -1, +1, 00], ##    _E
        [00, 00, 00, 00, -1, -1], ##    _F
    ])


if False:
    ## two minimal trivial cycles of length 3 with one shared edge
    ##
    # rxn0 : A = B
    # rxn1 : B = C
    # rxn2 : A = C
    # rxn3 : C = D
    # rxn4 : B = D
    ##
    # 0 = rxn0 + rxn1 - rxn2
    # 0 = rxn1 + rxn3 - rxn4
    ##
    #
    #  the cycles are:      A-B-C-A and D-B-C-D ,
    #  the shared edge is:    B-C
    ##
    stoich_matrix = np.array(
    [                         ## mol_...
    # rxn_... --->            ##     |
    # rxn_0   1   2,  3   4   ##     |
        [+1, 00 ,+1, 00, 00], ## mol_A
        [-1, +1, 00, 00, +1], ##    _B
        [00, -1, -1, +1, 00], ##    _C
        [00, 00, 00, -1, -1], ##    _D
    ])



## decision problems: which node to choose?

if False:
    ## two cycles of length 4 each, which 3 shared nodes
    ## -- the remaining two nodes close the cycles equally well
    ##
    #  rxn0 : A = B
    #  rxn1 : B = C
    #  rxn2 : C = D
    #  rxn3 : A = D
    #  rxn4 : B = E
    #  rxn5 : D = E
    ##
    ## indicate designed cycles:
    #  0 = rxn0 + rxn1 + rxn2 - rxn3
    #  0 = rxn0 + rxn4 - rxn5 - rxn3
    ##
    #  the cycles are: A-B-C-D-A and A-B-E-D-A
    #  Both paths are equivalent, there is no difference in whether to choose
    #  C or E as the bridge.
    ##
    stoich_matrix = np.array(
    [                             ## mol_...
    # rxn_... --->                ##     |
    # rxn_0   1   2   3   4   5   ##
        [+1, 00, 00, +1, 00, 00], ## mol_A
        [-1, +1, 00, 00, +1, 00], ##    _B
        [00, -1, +1, 00, 00, 00], ##    _C
        [00, 00, -1, -1, 00, +1], ##    _D
        [00, 00, 00, 00, -1, -1], ##    _E
    ])




## decision problems: which edges to choose?

if False:
    ## three cycles of length 3, with only one unique edge each -- after closing
    ## two cycles, the third can be constructed in two ways
    ##
    #  rxn0 : A = B
    #  rxn1 : B = C
    #  rxn2 : C = D
    #  rxn3 : A = D
    #  rxn4 : A = C
    #  rxn5 : B = D
    ##
    ## indicate designed cycles:
    #  0 = rxn0 - rxn3 + rxn5    == A-B-D-A
    #  0 = rxn1 + rxn2 - rxn5    == B-C-D-B
    #  0 =-rxn2 + rxn3 - rxn4    == C-A-D-C
    ##
    #  a further cycle is: A-B-C-A.
    #  After constructing, e.g. A-B-D-A and B-C-D-B, is it equally well to
    #  construct C-A-D-C or A-B-C-A? Both include the final missing edge (A-C),
    #  but they choose different nodes / edges again. in the C-A-D-C case,
    #  the node D is present in three cycles, whereas in the A-B-C-A case,
    #  the edges A-B and B-C are part of two cycles each.
    ##
    stoich_matrix = np.array(
    [                             ## mol_...
    # rxn_... --->                ##     |
    # rxn_0   1   2   3   4   5   ##
        [+1, 00, 00, +1, +1, 00], ## mol_A
        [-1, +1, 00, 00, 00, +1], ##    _B
        [00, -1, +1, 00, -1, 00], ##    _C
        [00, 00, -1, -1, 00, -1], ##    _D
    ])


## decision problems: what is the optimal solution -- more paths, or smaller ones?

if False:
    ## two cycles of length 4, 3 == three cycles of length 3, 3, 3
    ##
    #  rxn0 : A = B
    #  rxn1 : B = C
    #  rxn2 : C = D
    #  rxn3 : A = D
    #  rxn4 : A = C
    #  rxn5 : A = E
    #  rxn6 : C = E
    ##
    ## indicate designed cycles:
    #  0 = rxn0 + rxn1 + rxn2 - rxn3
    #  0 = rxn4 - rxn5 + rxn6
    #  0 = rxn0 - rxn1 - rxn4
    #  0 = rxn2 - rxn3 + rxn4
    ##
    #  the cycles are: A-B-C-D-A and A-C-E-A; or:
    #  could also be: A-B-C-A, A-C-D-A, A-C-E-A.
    #  In both, A-C-E-A is a shared cycle.
    #  Is it now better to duplicate an edge (A-C) to get smaller solutions,
    #  or is it better to have a larger solution (cycle of length 4)?
    ##
    stoich_matrix = np.array(
    [                                 ## mol_...
    # rxn_... --->                    ##     |
    # rxn_0   1   2   3   4   5   6   ##
        [+1, 00, 00, +1, +1, +1, 00], ## mol_A
        [-1, +1, 00, 00, 00, 00, 00], ##    _B
        [00, -1, +1, 00, -1, 00, +1], ##    _C
        [00, 00, -1, -1, 00, 00, 00], ##    _D
        [00, 00, 00, 00, 00, -1, -1], ##    _E
    ])


if True:
    ## two cycles of length 5, 4 == three cycles of length 4, 4, 3
    ##
    ## indicate reactions:
    #  rxn0 : A = B
    #  rxn1 : B = C
    #  rxn2 : C = D
    #  rxn3 : D = E
    #  rxn4 : E = F
    #  rxn5 : A = F
    #  rxn6 : B = G
    #  rxn7 : B = E
    #  rxn8 : E = G
    ##
    ## indicate designed cycles:
    #  0 = rxn0 + rxn6 - rxn8 + rxn4 + rxn5     == A-B-G-E-F-A
    #  0 = rxn1 + rxn2 + rxn3 - rxn7            == B-C-D-E-B
    ##
    #  possible cycles are also: B-E-G-B and A-B-E-F-A. This allows to construct
    #  e.g. (A-B-E-F-A, B-C-D-E-B, B-E-G-B) as another set of cycles of
    #  lengths (4,4,3).
    ##
    stoich_matrix = np.array(
    [                                         ## mol_...
    # rxn_... --->                            ##     |
    # rxn_0   1   2   3   4   5   6   7   8   ##
        [+1, 00, 00, 00, 00, +1, 00, 00, 00], ## mol_A
        [-1, +1, 00, 00, 00, 00, +1, +1, 00], ##    _B
        [00, -1, +1, 00, 00, 00, 00, 00, 00], ##    _C
        [00, 00, -1, +1, 00, 00, 00, 00, 00], ##    _D
        [00, 00, 00, -1, +1, 00, 00, -1, +1], ##    _E
        [00, 00, 00, 00, -1, -1, 00, 00, 00], ##    _F
        [00, 00, 00, 00, 00, 00, -1, 00, -1], ##    _G
    ])




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
    [                 ## mol_...
    # rxn_... --->    ##     |
    # rxn_0   1   2   ##     |
        [+1, 00, +1], ## mol_A
        [-1, 00, -1], ##    _B
        [00, +1, +1], ##    _C
        [00, -1, -1], ##    _D
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
    [                 ## mol_...
    # rxn_... --->    ##     |
    # rxn_0   1   2   ##     |
        [+1, 00, +1], ## mol_A
        [+1, 00, +1], ##    _B
        [-1, +1, 00], ##    _C
        [-1, +1, 00], ##    _D
        [00, -1, -1], ##    _E
        [00, -1, -1], ##    _F
    ])



if False:
    ## what is this?
    stoich_matrix = np.array(
    [                                         ## mol_...
    # rxn_... --->                            ##     |
    # rxn_0   1   2   3   4   5   6   7   8   ##
        [+1, +1 ,00 ,00, 00, 00, 00, -1, 00], ## mol_A
        [+1, 00, +2, 00, 00, 00, 00, -1, 00], ##    _B
        [-1, -1, 00, 00, 00, 00, 00, 00, 00], ##    _C
        [-1, 00, -2, 00, 00, 00, 00, 00, 00], ##    _D
        [00, +1, -2, 00, 00, 00, 00, 00, 00], ##    _E
        [00, -1, +2, 00, 00, 00, 00, 00, 00], ##    _F
        [00, 00, 00, -1, 00, -1, -1, +1, 00], ##    _G
        [00, 00, 00, +1, -1, 00, 00, 00, 00], ##    _H
        [00, 00, 00, +1, -1, 00, 00, 00, 00], ##    _I
        [00, 00, 00, 00, +1, +1, 00, 00, 00], ##    _J
        [00, 00, 00, 00, 00, 00, +1, 00, 00], ##    _K
        [00, 00, 00, 00, 00, 00, 00, 00, 00], ##    _L
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

    ## the reduced row echolon form (rref) is not the same for permutated columns
    for column_ordering in itertools.permutations(range(number_of_columns)):
        #logger.debug(f"order of columns to investigate: {column_ordering}")
        reordered_stoich_matrix = stoich_matrix[:,column_ordering]
        _, inds = sympy.Matrix(reordered_stoich_matrix).rref()
        #logger.debug(inds)
        #logger.debug(_)
        dependent_columns = sorted(list([column_ordering[i] for i in range(number_of_columns) if i not in inds]))
        logger.debug(dependent_columns)

        assert len(dependent_columns) == number_of_dependent_columns


    ## the reduced row echolon form (rref) is the same for permutated rows
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


### reduce matrix so it contains linearly independent columns only
matrix_complete = stoich_matrix.copy()
matrix_reduced  = stoich_matrix.copy()[:, [row_idx for row_idx in range(len(collected_rxns)) if not row_idx in linear_dependent_rxn_indices]]
collected_rxns_reduced = [r for r in collected_rxns if not r in linear_dependent_rxns]


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
linear_dependent_rxn_indices_to_calculate = linear_dependent_rxn_indices
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

    print("setting up the problem...")
    rxn_variables  = list([pulp.LpVariable(f"rxn_{rxn_id}", -3, 3, pulp.const.LpInteger) for rxn_id in collected_rxns_reduced])
    prob = pulp.LpProblem(f"solve_for_rxn_{collected_rxns[solve_for_rxn_id]}", pulp.const.LpMinimize)
    for row_idx in range(matrix_reduced.shape[0]):
        constraint_for_row = 0
        for column_idx, value in enumerate(matrix_reduced[row_idx]):
            constraint_for_row += value * rxn_variables[column_idx]
        prob += constraint_for_row == matrix_complete[row_idx, solve_for_rxn_id], f"row_{row_idx}"

    ## objective function is actually supposed to be abs(sum(rxn_vars)), but needs workaround
    abs_of_rxn_variables  = list([pulp.LpVariable(f"abs_of_rxn_{rxn_id}") for rxn_id in collected_rxns_reduced]) #, -3, 3, pulp.const.LpInteger
    objective_function = pulp.lpSum( [abs_of_rxn_var for abs_of_rxn_var in abs_of_rxn_variables] ) # if not abs_of_rxn_var is abs_of_rxn_variables[solve_for_rxn_id]] )
    ## workaround formula for abs, by replacing actual variable by dummy and defining constraints of actual<=>dummy:
    for abs_of_rxn_var, rxn_var in zip(abs_of_rxn_variables, rxn_variables):
        prob += rxn_var  <= abs_of_rxn_var
        prob += -rxn_var <= abs_of_rxn_var
    prob += objective_function
    print("...done.")

    print("Solving the problem...")
    prob.solve()
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
