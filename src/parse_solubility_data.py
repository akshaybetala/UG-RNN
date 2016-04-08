import numpy as np



def load_solubility_data(path = '../data/Delaney_solubility.txt'):
    f = file(path,'r')
    txt = f.read()
    f.close()
    #print txt
    lines = [line.split('"') for line in txt.split('\n')]
    parsed = []
    for line in lines[:-1]:
        if len(line)>1:
            parsed.append([x.replace(' ','') for x in ([line[1]]+line[2].split(',')) if x != ''])
        else:
            parsed.append([x.replace(' ','') for x in line[0].split(',')])
    
    
    if 0:
        print parsed[0]
        print '---------------------------'
        print parsed[1]
        print parsed[2]
        print parsed[372]
        print
    
    prediction_targets = np.array([float(x[1]) for x in parsed[1:]],dtype='float32')
    SMILES = [x[-1] for x in parsed[1:]]
    return SMILES, prediction_targets


if __name__=='__main__':
    SMILES, prediction_targets = load_solubility_data(path = 'Delaney_solubility.txt')
    
    print prediction_targets.min(), prediction_targets.mean(), prediction_targets.max()

    #print SMILES