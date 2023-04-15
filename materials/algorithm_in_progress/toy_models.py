#!/usr/bin/env python

import numpy
import numpy as np

##



## template:
def template():
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
    ##
    return stoich_matrix



## one-substrate one-product scenarios:

def a():
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
    ##
    return stoich_matrix



def b():
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
    ##
    return stoich_matrix



def c():
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
    ##
    return stoich_matrix



def d():
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
    ##
    return stoich_matrix



def e():
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
    ##
    return stoich_matrix



def f():
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
    ##
    return stoich_matrix



def g():
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
    ##
    return stoich_matrix



def h():
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
    ##
    return stoich_matrix



def i():
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
    ##
    return stoich_matrix

pass