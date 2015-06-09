# cn-aligner
Aligner + Levenshtein distance for confusion networks with strings.

Aligns a hypothesis consensus network to a reference string.

Can also be used to align a string to a string; the hypothesis string just has to be formulated as a CN.

For example:

def string_to_hyp(s):
    hyp = []
    words = s.strip().split()
    for w in words:
        sgmt = Segment(0, 0, [Edge(w, 1)])
        hyp.append(sgmt)
    return hyp
    
    
