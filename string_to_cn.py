## turn a string into a consensus net (hyp) (when you want to align plain string to string)
## or a reference structure (ref).

import align_consensus as ac

def string_to_hyp(s):
    hyp = []
    words = s.strip().split()
    for w in words:
        sgmt = Segment(0, 0, [Edge(w)])
        hyp.append(sgmt)
    return hyp

def string_to_ref(s):
    return s.strip().split()
