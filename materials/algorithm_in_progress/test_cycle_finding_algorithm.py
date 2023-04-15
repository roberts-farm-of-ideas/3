#!/usr/bin/env python

import collections
import json
import logging
import os
import pickle
import re
import shutil
import time


import numpy
import numpy as np
import ortools
from ortools.sat.python import cp_model
import pulp
import sympy

import toy_models

##

def sort_variables(foobar_tuple):
    what_to_split_with = ("_", "RXN")
    for splitter in what_to_split_with:
        match = [int(s) for s in foobar_tuple[0].split(splitter) if s.isdigit()]
        if len(match)==1:
            return match[0]
    return foobar_tuple

def variables_to_coefficient_array(selected_variables_for_this_solution, collected_rxns):
    coefficient_array = numpy.zeros((1,len(collected_rxns)),dtype="int")
    for i,r in enumerate(collected_rxns):
        selected_coefficients_for_this_reaction = [c for _,c in selected_variables_for_this_solution if _==r]
        assert len(selected_coefficients_for_this_reaction) <= 1
        if len(selected_coefficients_for_this_reaction) == 1:
            coefficient_array[0,i] = selected_coefficients_for_this_reaction[0]
    return coefficient_array

class PlaybackSolver():
    def __init__(self, selected_variables, status_to_return=cp_model.OPTIMAL):
        self.__variables = selected_variables
        self.__status_to_return = status_to_return

    def Value(self, variable):
        if variable in self.__variables:
            return 1 
        else:
            return 0
    
    def Solve(self, *args, **kwargs):
        print("The following status comes from the PlaybackSolver!")
        return self.__status_to_return 


class VarArraySolutionPrinter(cp_model.CpSolverSolutionCallback):
    """Print intermediate solutions."""

    def __init__(self, variables):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.__variables = variables
        self.__solution_count = 0
        self.__solution_storage_in_order_of_discovery = list()

    def on_solution_callback(self):
        self.__solution_count += 1
        foo = sorted([key for key in self.__variables if self.Value(self.__variables[key]) != 0], key=sort_variables)
        
        solution_was_already_discovered = False
        padders = [" "," "]
        if foo in self.__solution_storage_in_order_of_discovery:
            padders = ["(",")"]
            solution_was_already_discovered = True

        if not solution_was_already_discovered:
            self.__solution_storage_in_order_of_discovery.append(foo)


        print(f'{padders[0]}{str(self.__solution_count).rjust(3)}{padders[1]}: ', end="")
        #for v in self.__variables:
        #    print(f'{str(v)[:5]} = {str(self.Value(v)).rjust(2)}', end='  ')
        #print(" ".join([str(v)[4:5] for v in self.__variables if self.Value(v) != 0]))
        print(len([key for key in self.__variables if self.Value(self.__variables[key]) != 0]))
        #print()

        write_out_all_solutions_every_nth_solutions = 1
        if self.__solution_count % write_out_all_solutions_every_nth_solutions == 0:
            #with open("solution_storage_max.txt", "w") as f:
            #    f.write(repr(self.__solution_storage_in_order_of_discovery))
            pass


    def solution_count(self):
        return self.__solution_count

    def solution_storage_in_order_of_discovery(self):
        return self.__solution_storage_in_order_of_discovery
    
###


logging.basicConfig()
logger=logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

##TODO check ./sat/sat_parameters_pb2.pyi
from ortools.init import pywrapinit
#pywrapinit.CppBridge.InitLogging('integer_programming.py') 
cpp_flags = pywrapinit.CppFlags() 
cpp_flags.cp_model_dump_lns = False
cpp_flags.cp_model_dump_models =False
cpp_flags.cp_model_dump_prefix = ""
cpp_flags.cp_model_dump_response = False
cpp_flags.log_prefix = False
cpp_flags.stderrthreshold = 0#2
pywrapinit.CppBridge.SetFlags(cpp_flags)


if False:
    stoich_matrix = toy_models.a()
    #stoich_matrix = toy_models.b()
    #stoich_matrix = toy_models.c()
    #stoich_matrix = toy_models.d()
    #stoich_matrix = toy_models.e()
    #stoich_matrix = toy_models.f()
    #stoich_matrix = toy_models.g()
    #stoich_matrix = toy_models.h()
    stoich_matrix = toy_models.i()
    collected_mols = list( [f"mol_{chr(65+i)}" for i in range(stoich_matrix.shape[0])])
    collected_rxns = list( [f"rxn_{i}" for i in range(stoich_matrix.shape[1])])
    
else:
    ## load data from chemtunes biopath
    with open("./stoich_matrix_and_indices.npz", "rb") as f:
        loaded = numpy.load(f)
        stoich_matrix = loaded["stoich_matrix"]
        collected_mols = loaded["collected_mols"]
        collected_rxns = loaded["collected_rxns"]


## blacklist trivial reactions (back- and forward)
blacklist = []
import itertools
while False:
#for i,j in itertools.combinations(range(len(collected_rxns)),2):
    print(i,j)

    
    ##DEBUG
    ## in chemtunes biopath, rxn272 and rxn198 are a trivial pair
    #j=272-1
    #i=198-1

    ## quick version with numpy
    diff1 = +stoich_matrix[:,i] - stoich_matrix[:,j]
    diff2 = -stoich_matrix[:,i] + stoich_matrix[:,j]
    diff3 = +stoich_matrix[:,i] + stoich_matrix[:,j]
    all_zeros = not np.any(diff1) or not np.any(diff2) or not np.any(diff3)

    ## TEST!
    calc1 = stoich_matrix[:,i] / stoich_matrix[:,j]
    vals = [v for v in list(calc1) if v!=0]
    unique_vals = list(vals).unique()
    assert len(unique_vals)==1
    #all_zeros = not np.any(diff1) or not np.any(diff2) or not np.any(diff3)


    ## stop rref at 83 562 .... 2023-03-23 15.33

    ##DEBUG
    #print("---")
    #print(list(stoich_matrix[:,i]))
    #print(list(stoich_matrix[:,j]))
    #print("+++")
    #print(list(diff1))
    #print(all_zeros)

    if all_zeros:
        blacklist.append(frozenset([i,j]))
        print("BLACKLISTED:")
        print(collected_rxns[[i,j]])
        
        ##DEBUG
        #import pdb; pdb.set_trace()
        #print("exit this by entering on pdb command line:")
        #print("pdb.set_trace = lambda: None")
        #print("c")
        #break
        with open("blacklist.txt", "w") as f:
            f.write(repr(blacklist))

    continue 

    ##TIMING
    # 1min35s for reaching 2, 0
    # 2min08s for reaching 2, 1034 (first trivial pair)
    # 5min25s for reaching 6, 0 (~6*1545)
    # 7min30s for reaching 9, 0 (~10*1545)
    # ...         reaching 9, 1156 (second trivial pair)

    ### there are 1 192 740 combinations for 2 out of 1545 ->
    
    mat = sympy.Matrix(stoich_matrix[:,[i,j]])
    _,counts = mat.rref()
    #print(mat,counts)
    #print(counts)

    if counts == (0,):
        blacklist.append(frozenset([i,j]))
        print("BLACKLISTED:")
        print(collected_rxns[[i,j]])
        
        ##DEBUG
        #import pdb; pdb.set_trace()
        #print("exit this by entering on pdb command line:")
        #print("pdb.set_trace = lambda: None")
        #print("c")
        #break
        with open("blacklist.txt", "w") as f:
            f.write(repr(blacklist))

##DEBUG
#blacklist.append(frozenset([1,2]))
blacklist = set(blacklist)

## TODO: what about trivial triples, i.e. the blacklist would be [(1,2), (1,3)]; we would like to remove "all the others"
## WORKAROUND: can be checked by another run after first blacklist-removal
blacklist_indices = [list(i)[1] for i in blacklist]

with open("blacklist.txt", "w") as f:
    f.write(repr(blacklist))

print(f"I will reduce the stoich_matrix due to blacklisting by {len(blacklist)} columns.")
for foo in blacklist:
    print(collected_rxns[list(foo)])

stoich_matrix = numpy.delete(stoich_matrix,blacklist_indices,axis=1)
collected_rxns = numpy.delete(collected_rxns,blacklist_indices,axis=0)


print("saving stoichiometry matrix and indices...")
with open("./stoich_matrix_and_indices_after_blacklisting.npz", "wb") as f:
    numpy.savez(f, stoich_matrix=stoich_matrix, collected_rxns=collected_rxns, collected_mols=collected_mols)


## re-calculate names
#collected_mols = list( [f"mol_{chr(65+i)}" for i in range(stoich_matrix.shape[0])])
#collected_rxns = list( [f"rxn_{i}" for i in range(stoich_matrix.shape[1])])
#print(collected_mols[:2])

##DEBUG
#quit()



## ortools
from ortools.sat.python import cp_model

model = cp_model.CpModel()

coefficients = [-1,1]
## do NOT include 0 as a possible coefficient!

# choices[(r,c)]: choice (yes/no) for coefficient 'c' for reaction with id 'r'
#  
choices = {}
for r in collected_rxns:
    for c in coefficients:
        choices[(r,c)] = model.NewBoolVar(f"choice_{r}_{c}")


print("adding objective function...")
## try to find cycles with minimal numbers of reactions
model.Minimize(sum(choices[r,c] for r in collected_rxns for c in coefficients))


print("adding constraints...")

print(stoich_matrix.shape)
## only one coefficient can be active for each reaction
for i,r in enumerate(collected_rxns):
    model.AddAtMostOne(choices[(r,c)] for c in coefficients)

## for any thermodynamic cycle, all metabolites' stoichiometry coefficients must sum to 0
for metabolite_idx,_ in enumerate(collected_mols):
    constraint_for_metabolite = 0
    for counter, metabolite_stoichiometry_in_this_reaction in enumerate(stoich_matrix[metabolite_idx]):
        r = collected_rxns[counter]
        for c in coefficients:
            constraint_for_metabolite += metabolite_stoichiometry_in_this_reaction * c * choices[(r,c)]
    #print(constraint_for_metabolite)
    model.Add(constraint_for_metabolite == 0) #, f"constraint_for_metabolite_{metabolite_idx}"
    ##progress
    print(metabolite_idx, end=".",flush=True)

## exclude zero solution
print("adding zero solution constraint...")
model.AddAtLeastOne([choices[r,c] for r in collected_rxns for c in coefficients])


print("solve...", flush=True)
solver = cp_model.CpSolver()
solution_printer = VarArraySolutionPrinter(choices)

##TODO is this clever???
solver.parameters.enumerate_all_solutions = True

optimal_solution_storage_in_order_of_discovery = list()
combinatorially_obtained_solution_storage = list()
optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix = numpy.ndarray((0,len(collected_rxns)),dtype="int")

## Allows to resume at previous stage of solution
total_playback_mode = False
playback_mode = True
if playback_mode:
    print("Preparing playback mode...")
    with open("optimal_solution_storage.txt") as f:
        playback_optimal_solution_storage_in_order_of_discovery = eval(f.read())
        print(f"Will playback {len(playback_optimal_solution_storage_in_order_of_discovery)} solutions.")
    shutil.copyfile("optimal_solution_storage.txt", f"optimal_solution_storage.txt.{int(time.time())}")
if total_playback_mode:
    with open("optimal_solution_storage.txt") as f:
        optimal_solution_storage_in_order_of_discovery = eval(f.read())
    with open("combinatorially_obtained_solution_storage.txt") as f:
        combinatorially_obtained_solution_storage = eval(f.read())
    foobar = numpy.load("optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix.npz", allow_pickle=True)
    optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix = foobar["optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix"]
    #
    all_previous_solutions = optimal_solution_storage_in_order_of_discovery+combinatorially_obtained_solution_storage
    print(f"Adding now constraints for all previous solutions. This is going to be {len(all_previous_solutions)} constraints.")
    fractions = [float(pp)/100.0 for pp in range(0,100)]
    checkpoint_numbers = [int(len(all_previous_solutions)*i) for i in fractions]
    for i, solution in enumerate(all_previous_solutions):
        constraint_to_exclude_this_solution = sum(choices[key] for key in solution) <= len(solution)-1
        model.Add(constraint_to_exclude_this_solution)
        reverse_sign_choices = []
        for r, c in solution:
            reverse_sign_choices.append(choices[(r,-c)])
        constraint_to_exclude_this_solution = sum(reverse_sign_choices) <= len(solution)-1
        model.Add(constraint_to_exclude_this_solution)

        if i in checkpoint_numbers:
            print(f"I already added {i} constraints")
        



if not total_playback_mode:
    if os.path.exists("./parking_lot/") :
        shutil.rmtree("./parking_lot/")
    os.mkdir("./parking_lot/")

same_sense_problem_counter = 0
same_sense_problems = []
current_solution_length = 0
current_index_of_optimal_solution = -1
COMBINATORIAL_OFFSET = 1 ## how much space do you want ahead of the current solution length for excluding finding combinatorial solutions 
while True:
    if playback_mode:
        if current_index_of_optimal_solution+1 < len(playback_optimal_solution_storage_in_order_of_discovery):
            if total_playback_mode:
                current_index_of_optimal_solution += 1
                continue
            playback_this_solution = playback_optimal_solution_storage_in_order_of_discovery[current_index_of_optimal_solution+1]
            solver = PlaybackSolver([choices[(r,c)] for (r,c) in playback_this_solution])
            print("...doing playback here...")
        else:
            print("... leaving playback mode.")
            playback_mode = False
            solver = cp_model.CpSolver()
            breakpoint()

    status = solver.Solve(model, solution_printer)

    if status != cp_model.OPTIMAL:
        print("", flush=True)
        print("status:",status,solver.StatusName())
        break

    ## print this solution
    print()
    print(f"Another optimal solution (this is number {len(optimal_solution_storage_in_order_of_discovery)}):")
    selected_variables_for_this_solution = [key for key in choices if solver.Value(choices[key]) != 0]
    print(selected_variables_for_this_solution)

    ## exclude current solution for the next run
    selected_variables_for_this_solution = [key for key in choices if solver.Value(choices[key]) != 0]
    constraint_to_exclude_this_solution = sum(choices[key] for key in selected_variables_for_this_solution) <= len(selected_variables_for_this_solution)-1
    model.Add(constraint_to_exclude_this_solution)

    ## exclude reverse-sense version of this solution 
    selected_variables_for_this_solution = [key for key in choices if solver.Value(choices[key]) != 0]
    reverse_sign_choices = []
    for r, c in selected_variables_for_this_solution:
        reverse_sign_choices.append(choices[(r,-c)])
    constraint_to_exclude_this_solution = sum(reverse_sign_choices) <= len(selected_variables_for_this_solution)-1
    model.Add(constraint_to_exclude_this_solution)


    ## store also as coefficients, and check whether this is actually an optimal solution considering the currently available set of optimal solutions
    coefficient_array = variables_to_coefficient_array(selected_variables_for_this_solution, collected_rxns)
    optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix = numpy.append(optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix, coefficient_array, axis=0)
    #
    was_it_optimal = None
    selected_variables_for_this_solution = [key for key in choices if solver.Value(choices[key]) != 0]
    #print(optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix)
    print("Calculating the rref...")
    if playback_mode:
        rref,pivots=None,[None]*optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix.shape[0]
    else:
        rref,pivots= sympy.Matrix(optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix).T.rref()
    print("...done.")
    #print(rref, pivots,optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix.shape[0])
    if len(pivots)<optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix.shape[0]:
        was_it_optimal = False
        print("this was not an actual optimal solution, but just a linear combination of previous solutions")
        optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix = optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix[:-1,:]
        combinatorially_obtained_solution_storage.append(sorted(selected_variables_for_this_solution, key=sort_variables))
        breakpoint()
        continue
    else:
        was_it_optimal = True
        optimal_solution_storage_in_order_of_discovery.append(sorted(selected_variables_for_this_solution, key=sort_variables))
        current_index_of_optimal_solution +=1
    

    ## find previous optimal solutions that would just extend the cycle
    ## -- if necessary, also extend previous solutions
    selected_variables_for_this_solution = [key for key in choices if solver.Value(choices[key]) != 0]
    need_to_construct_new_combinations = False
    if len(selected_variables_for_this_solution) >= current_solution_length+COMBINATORIAL_OFFSET:
        print("I will construct new combinations for the previous parking_lots.")
        need_to_construct_new_combinations = True
    if len(selected_variables_for_this_solution) > current_solution_length:
        print(f"Setting current_solution_length to {len(selected_variables_for_this_solution)}.")
        current_solution_length = len(selected_variables_for_this_solution)
    #
    with open(f"./parking_lot/optimal_solution_{current_index_of_optimal_solution}.pickle","wb") as f:
            pickle.dump([sorted(selected_variables_for_this_solution, key=sort_variables)], f)
    #
    for dirpath,dirnames,filenames in os.walk("./parking_lot/"):
        filenames = sorted(filenames, key=lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', s)])
        print(dirpath,dirnames,filenames)

        for filename in filenames:
            f_nr = int(filename.split("_")[-1].split(".")[0])
            if not need_to_construct_new_combinations:
                # skip over other parking lots, if we don't need to look after them
                if f_nr != current_index_of_optimal_solution:
                    continue
            with open(dirpath + filename,"rb") as f:
                parking_lot = pickle.load(f)    
            queue = collections.deque(parking_lot)
            parking_lot = []
            print(f"Working on the parking_lot {filename}, which contained {len(queue)} element(s)")
            #
            print("Entering the queueing system...")
            while queue:
                ## we will add the current queue element to the constraints, and look for further extended versions of this elements

                selected_variables_for_this_solution = queue.pop()
                print("--")
                print("remaining_queue_length: ", len(queue))
                print("I'm looking at element: ", selected_variables_for_this_solution)

                ## this is relevant for need_to_construct_new_combinations -- if it's still too long, 
                ## move it to the parking lot
                if need_to_construct_new_combinations:
                    if len(selected_variables_for_this_solution) > current_solution_length+COMBINATORIAL_OFFSET:
                        print("putting the solution without touching on the parking lot") 
                        print("++")
                        parking_lot.append(newly_combined_solution)
                        continue


                ## should probably not happen at all
                if (not selected_variables_for_this_solution in combinatorially_obtained_solution_storage) and (not selected_variables_for_this_solution in optimal_solution_storage_in_order_of_discovery):
                    breakpoint()
                    constraint_to_exclude_this_solution = sum(choices[key] for key in selected_variables_for_this_solution) <= len(selected_variables_for_this_solution)-1
                    model.Add(constraint_to_exclude_this_solution)
                    combinatorially_obtained_solution_storage.append(sorted(selected_variables_for_this_solution, key=sort_variables))

                    ## also exclude reverse-sign extended solution
                    constraint_to_exclude_this_solution = sum(choices[(r,-c)] for r,c in selected_variables_for_this_solution) <= len(selected_variables_for_this_solution)-1
                    model.Add(constraint_to_exclude_this_solution)
                    combinatorially_obtained_solution_storage.append(sorted([(r,-c) for r,c in selected_variables_for_this_solution], key=sort_variables))

                previous_queue_length = len(queue)
                previous_parking_lot_length = len(parking_lot)

                for i, solution in enumerate(optimal_solution_storage_in_order_of_discovery[:f_nr]):
                    reactions_of_this_optimal_solution = [reaction_id for (reaction_id,c) in solution]    

                    same_sense = None
                    matches_in_this_solution = []
                    for (r,c) in selected_variables_for_this_solution:
                        if r in reactions_of_this_optimal_solution:
                            matches_in_this_solution.append((r,c))
                            matched_reaction_variable = [(matched_r,matched_c) for (matched_r,matched_c) in solution if matched_r == r]
                            assert len(matched_reaction_variable) == 1
                            matched_r, matched_c = matched_reaction_variable[0]
                            same_sense_this = True if c == matched_c else False
                            #print(f"same_sense_this={same_sense_this}")
                            ## make sure all reactions are consistently aligned (same or anti sense)
                            if same_sense is None:
                                same_sense = same_sense_this
                            else:
                                if same_sense != same_sense_this:
                                    #print("this solution cannot be used as an extension as a whole because the overlap is complicated")
                                    same_sense_problem_counter +=1
                                    same_sense_problems.append((selected_variables_for_this_solution, solution))
                                    ##TODO: one possibility now would be to split the solution into two extensions and apply them separately 
                                    continue
                    if len(matches_in_this_solution)>0:
                        assert same_sense is not None
                    else:
                        #print("not an overlapping solution...")
                        continue 
                    
                    matches_in_this_solution = sorted(matches_in_this_solution, key=sort_variables)
                    matches_in_this_solution = tuple(matches_in_this_solution)

                    print(f"Found overlap in solution {i} ... ({list(solution)})")
                    print(f"len(matches_in_this_solution): {len(matches_in_this_solution)}")
                    print(matches_in_this_solution)
                    
                    
                    newly_combined_solution = selected_variables_for_this_solution[:]
                    #print("start:")
                    #print(newly_combined_solution)
                    
                    for foo in matches_in_this_solution: newly_combined_solution.remove(foo)
                    #print("after removal:")
                    #print(newly_combined_solution)
                    
                    ## construct extension
                    if same_sense:
                        extension = [(reaction_id,-c) for (reaction_id,c) in solution if not reaction_id in [x for (x,_) in matches_in_this_solution]] 
                    else:
                        extension = [(reaction_id,c) for (reaction_id,c) in solution if not reaction_id in [x for (x,_) in matches_in_this_solution]] 
                    #print("extension:")
                    #print(extension)

                    ## overlap might be longer than "extension"
                    #if len(extension) < len(matches_in_this_solution):
                    #    print("This extension would be a reduction...")
                    ## could also be equal...
                    #if len(extension) == len(matches_in_this_solution):
                    #    print("This extension is not actually an extension...")
                    
                    newly_combined_solution.extend(extension)
                    newly_combined_solution = sorted(newly_combined_solution, key=sort_variables)
                    print("after extension:")
                    print(newly_combined_solution)

                    if newly_combined_solution in combinatorially_obtained_solution_storage:
                        print("discarding this solution from further processing because it was already discovered") 
                        print("+")
                        continue
                    elif newly_combined_solution in queue:
                        print("discarding this solution from further processing because it is already in queue") 
                        print("+")
                        ## occurs when solution was added twice to the queue on independent ways
                        continue

                    ## if it is extended "too much", we rather move it to a parking lot, to get done with the actual search with priority
                    if len(newly_combined_solution) > current_solution_length+COMBINATORIAL_OFFSET:
                        print("putting the extended solution on the parking lot") 
                        print("+")
                        parking_lot.append(newly_combined_solution)
                    else:
                        print("putting the extended solution back in the queue for further processing") 
                        print("+")
                        queue.append(newly_combined_solution)
                    
                    constraint_to_exclude_this_solution = sum(choices[key] for key in newly_combined_solution) <= len(newly_combined_solution)-1
                    model.Add(constraint_to_exclude_this_solution)
                    combinatorially_obtained_solution_storage.append(sorted(newly_combined_solution, key=sort_variables))

                    ## also exclude reverse-sign extended solution
                    constraint_to_exclude_this_solution = sum(choices[(r,-c)] for r,c in newly_combined_solution) <= len(newly_combined_solution)-1
                    model.Add(constraint_to_exclude_this_solution)
                    combinatorially_obtained_solution_storage.append(sorted([(r,-c) for r,c in newly_combined_solution], key=sort_variables))

                if previous_queue_length == len(queue) and previous_parking_lot_length == len(parking_lot):
                    print("I didn't add any extended solution to queue nor parking_lot.")

            if len(parking_lot) == 0:
                os.remove(f"./parking_lot/optimal_solution_{f_nr}.pickle")
            
            if len(parking_lot) > 0:
                print(f"I will put the remaining {len(parking_lot)} solutions on the parking lot file optimal_solution_{f_nr}.pickle, because they are longer than '{current_solution_length+COMBINATORIAL_OFFSET}'...")
                with open(f"./parking_lot/optimal_solution_{f_nr}.pickle","wb") as f:
                    pickle.dump(parking_lot, f)

            ## after processing each parking_lot, also write out the new combinatorially obtained solutions
            with open("combinatorially_obtained_solution_storage.txt", "w") as f:
                f.write(repr(combinatorially_obtained_solution_storage))
            
    
    with open("optimal_solution_storage.txt", "w") as f:
        f.write(repr(optimal_solution_storage_in_order_of_discovery))    
    with open("combinatorially_obtained_solution_storage.txt", "w") as f:
        f.write(repr(combinatorially_obtained_solution_storage))    
    numpy.savez_compressed("optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix.npz", optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix=optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix)



##TODO process the parking lot one final time
##.....

## write-out one final time
with open("optimal_solution_storage.txt", "w") as f:
    f.write(repr(optimal_solution_storage_in_order_of_discovery))    

with open("combinatorially_obtained_solution_storage.txt", "w") as f:
    f.write(repr(combinatorially_obtained_solution_storage))

for solution in combinatorially_obtained_solution_storage:
    #print(variables_to_coefficient_array(solution, collected_rxns))
    pass
print()

print("same sense problems: ", same_sense_problem_counter)
with open("same_sense_problems.txt", "w") as f:
    f.write(repr(same_sense_problems))

numpy.savez_compressed("optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix.npz", optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix=optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix)
print(optimal_solution_storage_in_order_of_discovery_as_coefficient_matrix)

## storage of solver callback containing discovered solutions (but not all combinatorially obtained ones) 
with open("solution_storage.txt", "w") as f:
    f.write(repr(solution_printer.solution_storage_in_order_of_discovery()))

print("done.")
print("number of all solutions:", len(combinatorially_obtained_solution_storage))
print("number of optimal solutions:", len(optimal_solution_storage_in_order_of_discovery))

"""
TIMING
0h05min: playback until solution 127, which is first cycle with 4 members
1h09min: solution 212, which is first cycle with 5 members
1h25min: solution 215, which is a combination
1h31min: solution 218, still 5 members
2h06min: solution 228, which is first cycle with 6 members; created 124,008 combinatorial solutions == 51MB on disk; calculating the parking lots since 20min... 
3h00min: ... working on parking_lot 4/96 ...
4h41min: ... working on parking_lot 9/96 ...
6h18min: currently consuming 12.4% memory (of 8GB)


TIMING:
0h06min: playback until solution 212 with "python ... > log"; 
0h12min: ... working on parking_lot 32/86 ...
0h19min: ... working on parking_lot 49/86 ...
???????: found solution 228, starting processing the 96 parking_lots
0h33min: ... working on parking_lot 10/96 ... 75MB combi storage
0h46min: ...                        16/96
1h04min: ...                        18/96, 109MB,  9.1% mem
4h18min: ...                        22/96, 272MB, 12.8% mem
11h50min: ..                        24/96, 364MB, 763349 combinations, 16,8% mem
"""
