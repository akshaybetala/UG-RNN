# AUTHOR: ALESSANDRO LUSCI
import sys

from pySources.TP import *


def main(argv):
    tp = TP(argv)
    tp.createTPList()
    tp.computeAverageValues()
    tp.joinLists()
    tp.computeMetrics()


if __name__ == "__main__":
    main(sys.argv)
