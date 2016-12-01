# Auhtor: Alessandro Lusci
import os
import sys

import openbabel
import pybel

import tools


def getRingAtomsMulti(molfilePath):
    for m in pybel.readfile("mol", molfilePath):
        # print m.OBMol.GetFormula()
        outString = ""
        rings = m.OBMol.GetSSSR()
        numR = 0;
        for r in rings:
            numR = numR + 1
        outString += "Ring count " + str(numR) + "\n"
        for r in rings:
            outString += "Ring size " + str(r.Size()) + "\n"
            path = r._path
            for p in path:
                outString += str(p - 1) + "\n"

        outString += "Atom count " + str(m.OBMol.NumAtoms()) + "\n"
        outString += "Index HowManyRings RingSize Hybridization Hydro_count Aromaticity AntiClockwise_chiral" + "\n"
        i = 0
        for a in openbabel.OBMolAtomIter(m.OBMol):
            outString += str(i) + " " + str(a.MemberOfRingCount()) + " " + str(a.MemberOfRingSize()) + " " + str(
                a.GetHyb()) + " " + str(a.ImplicitHydrogenCount()) + " "
            if a.IsAromatic():
                outString += "1" + " "
            else:
                outString += "0" + " "
            if a.IsClockwise():
                outString += "1" + "\n"
            else:
                outString += "0" + "\n"
            i = i + 1

    return outString


def generateInput(args):
    reduceGraph = args["reduceGraph"]
    addLogP = 1
    if reduceGraph:
        print "REDUCE GRAPH ON"
        # open files
    try:
        smileFile = open(args["inputFile"])
    except IOError:
        exit("CANNOT OPEN INPUT FILE: " + args["inputFile"])
    try:
        logPFile = open(args["inputFile"] + ".logP")
    except IOError:
        addLogP = 0
        print("NO LOGP FILE FOUND: " + args["inputFile"] + ".logP")
    try:
        inputFile = open(args["inputFile"] + ".mol", "w")
    except IOError, e:
        exit(e)
    ###########
    # count molecules in dataset
    smileList = smileFile.readlines()
    if addLogP == 1:
        logPList = logPFile.readlines()
    totalMolecules = len(smileList)
    if addLogP == 1 and totalMolecules != len(logPList):
        exit("ERROR: smileFile and logPFile contain a different number of molecules")
    inputFile.write(str(totalMolecules) + "\n")
    ###########################
    # main loop
    for l in range(0, len(smileList)):
        print "molecule " + str(l)
        smile = smileList[l].strip("\n\r")
        print smile
        target = "0.0"
        logP = "0.0"
        if addLogP == 1:
            logP = logPList[l].strip("\n")
        molecule = pybel.readstring("can", smile)
        # count atoms number in molecule
        atomsNumber = 0
        for atom in molecule:
            atomsNumber += 1
        if atomsNumber >= 100:
            print "Warning molecule " + str(l) + " contains " + str(atomsNumber) + " atoms"
            exit("EXITING")
        ###############################
        inputFile.write("train" + str(l) + ".can")
        molecule.write("mol", "temp", "overwrite=true")
        if reduceGraph == 0:
            try:
                temp = open("temp", "r")
            except IOError:
                exit("CANNOT OPEN TEMP FILE")
                # delete M RAD and M CHG lines
            tools.deleteLines(os.getcwd())
            # standardize molfile
            tools.standardizeMol(os.getcwd())
            inputFile.write(temp.read() + "$$$$\n" + target + "\n")
            inputFile.write(logP + "\n")
        else:
            # try:
            #	temp = open("temp","r")
            # except IOError:
            #	exit("CANNOT OPEN TEMP FILE")
            ringsInfo = getRingAtomsMulti("temp")
            print ringsInfo
            # standardize molfile
            tools.standardizeMol(os.getcwd())
            # temp.close()
            try:
                temp = open("temp", "a")
            except IOError:
                exit("CANNOT OPEN TEMP FILE")
            temp.write("$$$$\n" + target + "\n" + ringsInfo)
            temp.close()
            # delete M RAD and M CHG lines
            tools.deleteLines(os.getcwd())
            # reduce graph here
            tools.reduceGraph(os.getcwd())
            try:
                tempReduced = open("tempReduced", "r")
            except IOError:
                exit("CANNOT OPEN " + tempReduced + " FILE")

            inputFile.write(tempReduced.read())
            inputFile.write(logP + "\n")
    return len(smileList)


##########
#####

def main(argv):
    generateDataset(argv)


if __name__ == "__main__":
    main(sys.argv)
