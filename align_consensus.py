#!/usr/bin/python

## Does alignment given a reference transcript and a consensus network.
## Pruning the CN to n-best is not handled here;
## there should be a separate main script for that.
## this handles alignment of a reference to any arbitrary CN.
## likewise the scoring will be handled elsewhere.

import copy

####### CONSENSUS NET DATA STRUCTURES #######

## consensus net: [list of segments]

## segment: a 'sausage link.'  
## segment.start: start time
## segment.end: end time
## segment.edges: [list of edges]
class Segment:
    def __init__(self, start, end, edges): 
        self.start = start
        self.end = end
        self.edges = edges
    def __str__(self):
        ## just print the 1-best for testing purposes.
        return self.edges[0].label

## make sure anything using magical indices is turned
## into the right kind of object.

## edge: one edge in a consensus net.
## edge.label: the hypothesized word.
## edge.prob: that word's posterior.
## edge.metadata: optional; any additional info you want to carry
## around with the edges.  It can be anything at all.
class Edge:
    def __init__(self, label, prob, meta=None):
        self.label = label
        self.prob = prob
        self.metadata = meta

## deault gap segment for constructing the final alignment
hyp_gap_sgmt = Segment(-1, -1, [Edge('<None>', -1)])

## assume that preprocessing has happened: 
## transcript has been stripped of non-words
## and ~SIL, <s>, </s>, and <HES> have been converted to <epsilon>.

## TODO think some more about what to do with non-lexical items....

## refs look the same
## edges are now a list of Segments

####### LEVENSHTEIN DISTANCE + ALIGNMENT SUBROUTINES #######

## cost = the best cost we found for this table entry
## bp_i and bp_j are the backpointer coordinates
class Entry:
    def __init__(self, cost, bp_i, bp_j):
        self.cost = cost
        self.bp_i = bp_i
        self.bp_j = bp_j

## initialize the table 
def init_table(ref, hyp):
    result = []
    for i in range(0,len(ref)+1):
        result.append([])
        for j in range(0,len(hyp)+1):
            result[i].append(Entry(None,None,None))
    # initialize first row and column
    result[0][0] = Entry(0,0,0)
    for i in range(1,len(ref)+1):
        result[i][0] = Entry(i, i-1, 0)        
    for j in range(1,len(hyp)+1):
        result[0][j] = Entry(j, 0, j-1)
    return result

## Fill in the table
def lev(ref, hyp, table):
    # ref is a list of words
    # hyp is a list of tuples, the first item of each is the hyp word
    i = len(ref)
    j = len(hyp)
    # first check if it's in the table - if so, return.
    if table[i][j].cost != None:
        return table[i][j].cost
    # Substitution or Match
    # If ref word is in the hyp list,
    # cost = 1-p(matching edge).
    # epsilon is not a proper substitution, it's just wrong,
    # so not allowing those as a match.
    cost = 0
    match = [x for x in hyp[-1].edges if x.label == ref[-1]]
    if match == []:
        cost = 1 # complete substitution
    else:
        cost = 1 - float(match[0].prob) # exact match.
        # floating point weirdnesses: you can get 1.0001
        # so put that cost back to 0 instead of a negative #.
        if cost < 0:
            cost = 0
    # 3 possibilities for the recursion
    sub_or_equal_cost = lev(ref[0:-1], hyp[0:-1], table) + cost
    # Insertion: cost = 1 if no <epsilon> in hyp
    # else cost = 1-P(epsilon) (we are 'allowing' this insertion)
    ins_penalty = 0
    epsilons = [x for x in hyp[-1].edges if x.label == '<epsilon>']
    if epsilons == []:
        ins_penalty = 1
    else:
        ins_penalty = 1 - float(epsilons[0].prob)
    insertion_cost = lev(ref, hyp[0:-1], table) + ins_penalty
    # Deletion: always costs 1.
    deletion_cost = lev(ref[0:-1], hyp, table) + 1
    best = min(insertion_cost, deletion_cost, sub_or_equal_cost)
    # place the best one in the table with its corresponding backpointer
    # and return the result.
    # prioritize match > insertions > deletions > substitutions.
    if sub_or_equal_cost == best and cost == 0:
        table[i][j] = Entry(best, i-1, j-1)
    elif insertion_cost == best:
        table[i][j] = Entry(best, i, j-1)
    elif deletion_cost == best:
        table[i][j] = Entry(best, i-1, j)
    elif sub_or_equal_cost == best and cost > 0:
        table[i][j] = Entry(best, i-1, j-1)
    else:
        print 'something has gone very wrong!'
        print 'best: ' + str(best)
        print 'sub/eql cost: ' + str(sub_or_equal_cost)
        print 'ins cost: ' + str(insertion_cost)
        print 'deletion cost: ' + str(deletion_cost)
        print 'cost cost: ' + str(cost)
        exit()
    return best

# backtrack + build up the alignment
# the alignment will take the form of:
# ( [reference list of words including <None> for gaps], 
#   [hypothesis list of Segments with hyp_gap_sgmt for gaps] )
# we are carrying around the timestamps with the hypotheses
# because we will need them to determine the +/- labels.
def backtrack(ref, hyp, table):
    ref_copy = copy.deepcopy(ref)
    hyp_copy = copy.deepcopy(hyp) 
    ref_res = []
    hyp_res = []
    last_entry = table[-1][-1]
    backpointer = (last_entry.bp_i, last_entry.bp_j)
    current_cell = (len(table)-1, len(table[0])-1)
    while current_cell != (0, 0):
        # pointing diagonally: same or substitution
        if backpointer == (current_cell[0]-1, current_cell[1]-1):
            ref_res.insert(0, ref_copy.pop())
            hyp_res.insert(0, hyp_copy.pop())
        # pointing to the left:
        # it's an insertion; add a gap to the reference.
        elif backpointer == (current_cell[0], current_cell[1]-1):
            ref_res.insert(0, '<None>')
            hyp_res.insert(0, hyp_copy.pop())
        # pointing up:
        # it's a deletion; add a gap to the hypothesis.
        elif backpointer == (current_cell[0]-1, current_cell[1]): 
            ref_res.insert(0, ref_copy.pop())
            hyp_res.insert(0, copy.deepcopy(hyp_gap_sgmt))
        else:
            print 'something different has gone very wrong!'
            print 'current cell: ' + str(current_cell)
            print 'backpointer: ' + str(backpointer)
            break
        # update cells
        current_cell = backpointer
        bp_cell = table[current_cell[0]][current_cell[1]]
        backpointer = (bp_cell.bp_i, bp_cell.bp_j)
    return (ref_res, hyp_res)
