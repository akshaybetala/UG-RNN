import string


def reduceGraph(molfilePath):
    try:
        molfile = open(molfilePath + "/temp")
    except IOError:
        exit("CANNOT OPEN " + molfilePath + "/temp")
    try:
        molfileReduced = open(molfilePath + "/tempReduced", "w")
    except IOError:
        exit("CANNOT OPEN " + molfilePath + "/tempReduced")
        # get number of atoms
    line = molfile.readlines()
    info = line[3].split(" ");
    # print info
    numOfAtoms = -1
    numOfBonds = -1
    for i in range(0, len(info)):
        if info[i] != '':
            if numOfAtoms == -1:
                numOfAtoms = int(info[i])
            elif numOfBonds == -1:
                numOfBonds = int(info[i])
                # print numOfAtoms
                # print numOfBonds
    atoms = line[4:(4 + numOfAtoms)]
    # print atoms
    bonds = line[(4 + numOfAtoms):(4 + numOfAtoms + numOfBonds)]
    # print bonds
    print line[7 + numOfAtoms + numOfBonds]
    numOfRings = int(line[7 + numOfAtoms + numOfBonds].split(" ")[2])
    # print numOfRings
    if numOfRings != 0:
        rings = []
        seek = 7 + numOfAtoms + numOfBonds
        for i in range(0, numOfRings):
            rings.append([])
            seek = seek + 1
            atomsInRing = int(line[seek].split(" ")[2])
            for j in range(0, atomsInRing):
                seek = seek + 1
                rings[i].append(int(line[seek].strip("\n")) + 1)
                # print rings
        atomsSymb = []
        for i in range(0, len(atoms)):
            atomsSymb.append(atoms[i][31])
            # print atomsSymb
            # add rings to atoms
            # now just based on numberOfAtoms in ring but we should build a database of rings
            # update atoms
        atomsNew = []
        atomsPos = {}
        connectedRings = []
        # for i in range(0,len(atoms)):
        #	atomsPos.append([])
        for i in range(0, numOfRings):
            atomLine = "    0.0000    0.0000    0.0000 R" + str(len(rings[i])) + "  0  0  0  0  0\n"
            atomsNew.append(atomLine)
        seek = numOfRings - 1
        for i in range(0, len(atoms)):
            inRing = 0
            currRing = -1
            atomInRing = []
            for j in range(0, numOfRings):
                if i + 1 in rings[j]:
                    # print "atom "+str(i+1)+" in ring "+str(j+1)
                    atomsPos[str(i + 1)] = [i + 1, j + 1]
                    inRing += 1
                    currRing = j + 1
                    atomInRing.append(j + 1)
            atomInRing = sorted(atomInRing)
            if inRing >= 2:
                # create couples
                l = 0
                # print atomInRing
                for k in range(l + 1, len(atomInRing)):
                    if [atomInRing[l], atomInRing[k]] not in connectedRings:
                        connectedRings.append([atomInRing[l], atomInRing[k]])
                    l += 1
            if inRing == 0:
                # print "atom "+str(i+1)+" not in ring"
                atomsNew.append(atoms[i])
                seek += 1
                atomsPos[str(i + 1)] = [i + 1, seek + 1]

                # print rings
                # print connectedRings
                # print atomsNew
                # print atomsPos

                # update bonds
        bondsNew = []
        bondsNewTemp = []
        # insert bonds between rings
        # bond between rings is wighted with 0
        # for i in range(0,len(connectedRings)):
        # bondsNew.append("  "+str(connectedRings[i][0])+"  "+str(connectedRings[i][1]) +"  0  0  0  0\n")
        # bondsNewTemp.append([connectedRings[i][0]),connectedRings[i][1]])
        # update bonds between atoms
        # HERE: debug this part but we are almost done
        for i in range(0, len(bonds)):
            bond = bonds[i].split(" ")
            # extract atom info
            atom1 = 0
            atom2 = 0
            bondType = 0
            for j in range(0, len(bond)):
                # print bond
                if bond[j] != '':
                    if atom1 == 0:
                        atom1 = int(bond[j])
                    elif atom2 == 0:
                        atom2 = int(bond[j])
                    elif bondType == 0:
                        bondType = int(bond[j])

                        # check if atoms are in ring
                        # print bond
            if [atomsPos[str(atom1)][1], atomsPos[str(atom2)][1]] not in connectedRings and [atomsPos[str(atom2)][1],
                                                                                             atomsPos[str(atom1)][
                                                                                                 1]] not in connectedRings and \
                            atomsPos[str(atom1)][1] != atomsPos[str(atom2)][1]:
                bondsNewTemp.append([atomsPos[str(atom1)][1], atomsPos[str(atom2)][1], bondType])

                # append connectedRings
        for i in range(0, len(connectedRings)):
            connectedRings[i].append(0)
        for i in range(0, len(connectedRings)):
            bondsNewTemp.append(connectedRings[i])
        print bondsNewTemp
        print connectedRings
        # sort bondsNewTemp
        bondsNewTemp = sorted(bondsNewTemp)
        # print bondsNewTemp

        for i in range(0, len(bondsNewTemp)):
            if bondsNewTemp[i][0] > 9 and bondsNewTemp[i][1] <= 9:
                bondsNew.append(" " + str(bondsNewTemp[i][0]) + "  " + str(bondsNewTemp[i][1]) + "  " + str(
                    bondsNewTemp[i][2]) + "  0  0  0\n")
            if bondsNewTemp[i][0] <= 9 and bondsNewTemp[i][1] > 9:
                bondsNew.append("  " + str(bondsNewTemp[i][0]) + " " + str(bondsNewTemp[i][1]) + "  " + str(
                    bondsNewTemp[i][2]) + "  0  0  0\n")
            if bondsNewTemp[i][0] <= 9 and bondsNewTemp[i][1] <= 9:
                bondsNew.append("  " + str(bondsNewTemp[i][0]) + "  " + str(bondsNewTemp[i][1]) + "  " + str(
                    bondsNewTemp[i][2]) + "  0  0  0\n")
            if bondsNewTemp[i][0] > 9 and bondsNewTemp[i][1] > 9:
                bondsNew.append(" " + str(bondsNewTemp[i][0]) + " " + str(bondsNewTemp[i][1]) + "  " + str(
                    bondsNewTemp[i][2]) + "  0  0  0\n")
    else:
        atomsNew = atoms
        bondsNew = bonds
        # write molfile reduced
        # print atomsNew
        # print bondsNew
    molfileReduced.write(line[0])
    molfileReduced.write(line[1])
    molfileReduced.write(line[2])
    if len(atomsNew) <= 9 and len(bondsNew) <= 9:
        molfileReduced.write(
            "  " + str(len(atomsNew)) + "  " + str(len(bondsNew)) + "  0  0  0  0  0  0  0  0999 V2000\n")
    if len(atomsNew) <= 9 and len(bondsNew) > 9:
        molfileReduced.write(
            "  " + str(len(atomsNew)) + " " + str(len(bondsNew)) + "  0  0  0  0  0  0  0  0999 V2000\n")
    if len(atomsNew) > 9 and len(bondsNew) <= 9:
        molfileReduced.write(
            " " + str(len(atomsNew)) + "  " + str(len(bondsNew)) + "  0  0  0  0  0  0  0  0999 V2000\n")
    if len(atomsNew) > 9 and len(bondsNew) > 9:
        molfileReduced.write(
            " " + str(len(atomsNew)) + " " + str(len(bondsNew)) + "  0  0  0  0  0  0  0  0999 V2000\n")
    for i in range(0, len(atomsNew)):
        molfileReduced.write(atomsNew[i])
    for i in range(0, len(bondsNew)):
        molfileReduced.write(bondsNew[i])
    molfileReduced.write(line[4 + numOfAtoms + numOfBonds])
    molfileReduced.write(line[5 + numOfAtoms + numOfBonds])
    molfileReduced.write(line[6 + numOfAtoms + numOfBonds])

    # close files
    molfile.close()
    molfileReduced.close()


def deleteLines(molfilePath):
    try:
        molfile = open(molfilePath + "/temp")
    except IOError:
        exit("CANNOT OPEN " + molfilePath + "/temp")
    line = molfile.readlines()
    molfile.close()
    # print line
    length = len(line)
    for i in range(0, length):
        if string.find(line[i], 'RAD') != -1:
            # print line[i]
            line.remove(line[i])
            break
    length = len(line)
    for i in range(0, length):
        if string.find(line[i], 'CHG') != -1:
            # print line[i]
            line.remove(line[i])
            break
    try:
        molfile = open(molfilePath + "/temp", "w")
    except IOError:
        exit("CANNOT OPEN " + molfilePath + "/temp")
    molfile.write("".join(line))
    molfile.close()


def standardizeMol(datasetName):
    try:
        tempFile = file(datasetName + "/temp")
    except IOError, e:
        exit(e)
    tempList = tempFile.readlines()
    tempFile.close()
    atomsLine = tempList[3].split(" ")

    atoms = 0
    for el in atomsLine:
        if el != "" and atoms != 0:
            bonds = int(el)
            break
        if el != "":
            atoms = int(el)

    for i in range(1, atoms + 1):
        count = 0
        atomLine = tempList[3 + i].split(" ")
        for j in range(0, len(atomLine)):
            if atomLine[j] != "":
                count += 1
            if count == 9:
                tempList[3 + i] = " ".join(atomLine[0:j + 1])
                break

    for i in range(atoms + 1, atoms + 1 + bonds):
        count = 0
        bondLine = tempList[3 + i].split(" ")
        for j in range(0, len(bondLine)):
            if bondLine[j] != "":
                count += 1
            if count == 6:
                tempList[3 + i] = " ".join(bondLine[0:j + 1])
                break

    for i in range(0, len(tempList)):
        tempList[i] = tempList[i].strip("\n")
    tempString = "\n".join(tempList[0:len(tempList)])

    try:
        tempFile = open(datasetName + "/temp", "w")
        tempFile.write(tempString + "\n")
        tempFile.close()
    except IOError, e:
        exit(e)
