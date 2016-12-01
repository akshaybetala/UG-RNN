from Metrics import *


class TP:
    testID = ""
    models = 0
    folds = 0
    val = 0
    testFixed = 0
    averageOverFold = 0
    wx = 0

    def __init__(self, argv):
        usageString = "USAGE: python ugrnnCheckResults.py -t testID -m numberOfModels -f numberOfFolds -vl -aof";
        for i in range(0, len(argv)):
            if (argv[i] == "-t"):
                self.testID = argv[i + 1]
            elif (argv[i] == "-m"):
                self.models = int(argv[i + 1])
            elif (argv[i] == "-f"):
                self.folds = int(argv[i + 1])
            elif (argv[i] == "-tf"):
                self.testFixed = 1
            elif (argv[i] == "-vl"):
                self.val = 1
            elif (argv[i] == "-aof"):
                self.averageOverFold = 1
            elif (argv[i] == "-wx"):
                self.wx = 1
        if self.testID == "" or self.models == 0 or self.folds == 0:
            exit(usageString)

    def createTPList(self):
        t = []
        p = []
        if self.val != 1:
            flagFileTail = "/flag"
            logTPFileTail = "/logTP"
        else:
            flagFileTail = "/flagVal"
            logTPFileTail = "/logTPVal"
        for f in range(0, self.folds):
            t.append([])
            p.append([])
            for m in range(0, self.models):
                t[f].append([])
                p[f].append([])
                # load tp files
        completedFolds = 0
        for f in range(0, self.folds):
            a = 0
            for m in range(0, self.models):
                try:
                    flagFile = open(self.testID + "/fold" + str(f) + "/model" + str(m) + flagFileTail)
                except IOError:
                    print("CANNOT OPEN FLAG FILE")
                    a = 1
                if a == 1:
                    break
                else:
                    completedFolds += 1
                    try:
                        logTPFile = open(self.testID + "/fold" + str(f) + "/model" + str(m) + logTPFileTail)
                    except IOError:
                        exit("CANNOT OPEN logTP FILE")
                    logTPList = logTPFile.readlines()
                    for i in range(0, len(logTPList)):
                        temp = logTPList[i].split("\t")
                        temp[1] = temp[1].strip("\n")
                        p[f][m].append(float(temp[0]))
                        t[f][m].append(float(temp[1]))
            if a == 1:
                break
        completedFolds = int(float(completedFolds) / float(self.models))
        print "Completed folds: " + str(completedFolds)
        if completedFolds == 0:
            exit("Exiting...completedFolds==0")
        self.completedFolds = completedFolds
        self.t = t
        self.p = p

    def computeAverageValues(self):
        aP = []
        aT = []
        for i in range(0, self.completedFolds):
            aP.append([])
            aT.append([])
        for f in range(0, self.completedFolds):
            aT[f] = self.t[f][0]
            for i in range(0, len(self.p[f][0])):
                averageP = 0.0
                for m in range(0, self.models):
                    averageP += self.p[f][m][i]
                averageP = averageP / float(self.models)
                aP[f].append(averageP)
        self.aP = aP
        self.aT = aT

    def joinLists(self):
        tT = []
        pT = []
        for f in range(0, self.completedFolds):
            for i in range(0, len(self.t[f][0])):
                tT.append(self.t[f][0][i])
            for j in range(0, len(self.aP[f])):
                pT.append(self.aP[f][j])
        tTM = []
        pTM = []
        for m in range(0, self.models):
            tTM.append([])
            pTM.append([])
            for f in range(0, self.completedFolds):
                for i in range(0, len(self.t[f][0])):
                    tTM[m].append(self.t[f][m][i])
                    pTM[m].append(self.p[f][m][i])
        self.tT = tT
        self.pT = pT
        self.tTM = tTM
        self.pTM = pTM

    def computeMetricsTestFixed(self):
        foldLength = len(self.tT) / self.folds
        tTFixed = []
        pTFixed = []
        pTFixedAv = []
        tTFixedAv = []
        for i in range(0, foldLength):
            pTFixedAv.append(0.0)
            tTFixedAv.append(0.0)
        averageRMSE = 0.0
        averageAAE = 0.0
        averagePEARSON = 0.0
        averageStdAAE = 0.0
        for f in range(0, self.folds):
            first = f * foldLength
            last = first + foldLength
            tTFixed.append(self.tT[first:last])
            pTFixed.append(self.pT[first:last])
            for i in range(0, foldLength):
                pTFixedAv[i] = pTFixedAv[i] + pTFixed[f][i]
                tTFixedAv[i] = tTFixedAv[i] + tTFixed[f][i]
            metric = Metrics(tTFixed[f], pTFixed[f])
            metric.computeRMSE()
            metric.computeAAE()
            if len(tTFixed) > 1:
                metric.computePEARSON()
                metric.computeStdAAE()
                print "Fold" + str(f) + "\nRMSE:" + str(metric.RMSE) + " AAE:" + str(metric.AAE) + " PEARSON:" + str(
                    metric.PEARSON) + " STD_R:" + str(metric.stdAAE)
                averageRMSE += metric.RMSE
                averageAAE += metric.AAE
                averagePEARSON += metric.PEARSON
                averageStdAAE += metric.stdAAE
                print "FOLD AVERAGE...\nRMSE:" + str(averageRMSE / float(self.folds)) + " AAE:" + str(
                    averageAAE / float(self.folds)) + " PEARSON:" + str(
                    averagePEARSON / float(self.folds)) + " STD_R:" + str(averageStdAAE / float(self.folds))
            else:
                print "Fold" + str(f) + "\nRMSE:" + str(metric.RMSE) + " AAE:" + str(metric.AAE)

        pTFixedAv = [x / float(self.folds) for x in pTFixedAv]
        tTFixedAv = [x / float(self.folds) for x in tTFixedAv]
        if len(pTFixedAv) > 1:
            metric = Metrics(tTFixedAv, pTFixedAv)
            metric.computeRMSE()
            metric.computeAAE()
            metric.computePEARSON()
            metric.computeStdAAE()
            print "AVERAGE...\nRMSE:" + str(metric.RMSE) + " AAE:" + str(metric.AAE) + " PEARSON:" + str(
                metric.PEARSON) + " STD_R:" + str(metric.stdAAE)

    def computeMetricsAverageOverFold(self):
        RMSE = []
        AAE = []
        PEARSON = []

        for f in range(0, self.completedFolds):
            metric = Metrics(self.aT[f], self.aP[f])
            metric.computeRMSE()
            RMSE.append(metric.RMSE)
            metric.computeAAE()
            AAE.append(metric.AAE)
            if len(self.aT[f]) > 1:
                metric.computePEARSON()
                PEARSON.append(metric.PEARSON)
                print "Fold" + str(f) + "\nRMSE:" + str(RMSE[f]) + " AAE:" + str(AAE[f]) + " PEARSON:" + str(PEARSON[f])
            else:
                print "Fold" + str(f) + "\nRMSE:" + str(RMSE[f]) + " AAE:" + str(AAE[f])

        averageRMSE = 0.0
        averagePEARSON = 0.0
        averageAAE = 0.0
        for f in range(0, self.completedFolds):
            averageRMSE += RMSE[f]
            averageAAE += AAE[f]
            averagePEARSON += PEARSON[f]
        averageRMSE = averageRMSE / float(self.completedFolds)
        averageAAE = averageAAE / float(self.completedFolds)
        averagePEARSON = averagePEARSON / float(self.completedFolds)

        stdRMSE = 0.0
        stdPEARSON = 0.0
        stdAAE = 0.0
        for f in range(0, self.completedFolds):
            stdRMSE += math.pow((RMSE[f] - averageRMSE), 2)
            stdAAE += math.pow((AAE[f] - averageAAE), 2)
            stdPEARSON += math.pow((PEARSON[f] - averagePEARSON), 2)
        stdRMSE = math.sqrt(stdRMSE / float(self.completedFolds))
        stdAAE = math.sqrt(stdAAE / float(self.completedFolds))
        stdPEARSON = math.sqrt(stdPEARSON / float(self.completedFolds))

        print "FOLD AVERAGE...\nRMSE:" + str(averageRMSE) + " " + str(stdRMSE) + " AAE:" + str(averageAAE) + " " + str(
            stdAAE) + " PEARSON:" + str(averagePEARSON) + " " + str(stdPEARSON)

    def computeMetrics(self):
        if self.testFixed and self.completedFolds == self.folds:
            self.computeMetricsTestFixed()
        elif self.averageOverFold:
            self.computeMetricsAverageOverFold()
        else:
            for m in range(0, self.models):
                metric = Metrics(self.tTM[m], self.pTM[m])
                metric.computeRMSE()
                metric.computeAAE()
                if len(self.tTM[m]) > 1:
                    metric.computePEARSON()
                    print "Model" + str(m) + "\nRMSE:" + str(metric.RMSE) + " AAE:" + str(
                        metric.AAE) + " PEARSON:" + str(metric.PEARSON)
                else:
                    print "Model" + str(m) + "\nRMSE:" + str(metric.RMSE) + " AAE:" + str(metric.AAE)
            metric = Metrics(self.tT, self.pT)
            metric.computeRMSE()
            metric.computeAAE()
            if len(self.tT) > 1:
                metric.computePEARSON()
                metric.computeStdAAE()
                metric.computeStd()
                print "AVERAGE...\nRMSE:" + str(metric.RMSE) + " AAE:" + str(metric.AAE) + " PEARSON:" + str(
                    metric.PEARSON)

    def writeXls(self):
        self.folds = self.completedFolds
        if self.wx == 1 and self.completedFolds == self.folds:
            # open file
            xlsFile = open(self.testID + "/results.xls", "w")
            # open smile file
            smileFile = open(self.testID + "/smiles")
            smileList = smileFile.readlines()
            smileList = [x.strip("\n") for x in smileList]
            smileFile.close()
            # get smile positions
            smilePos = []
            for f in range(0, self.folds):
                # read smile positions
                smilePosFile = open(self.testID + "/fold" + str(f) + "/test")
                smilePos.append(smilePosFile.readlines()[1].split(" "))
                smilePosFile.close()
                # write results
            for f in range(0, self.folds):
                xlsFile.write(
                    "smile_fold" + str(f + 1) + "\tprediction_fold" + str(f + 1) + "\ttarget_fold" + str(f + 1) + "\t")
            xlsFile.write("\n")
            for m in range(0, len(self.aT[0])):
                for f in range(0, self.folds):
                    try:
                        xlsFile.write(smileList[int(smilePos[f][m])] + "\t")
                        xlsFile.write(str(self.aP[f][m]) + "\t")
                        xlsFile.write(str(self.aT[f][m]) + "\t")
                    except:
                        xlsFile.write("\t\t\t")
                xlsFile.write("\n")
            xlsFile.close()
