# function find
def find(item, seq):
    i = 0
    while i < len(seq):
        if item == int(seq[i]):
            return 1
        i += 1
    return 0
    ##############
