# file fastQCres.py

# This is a single fastq file's QC results from fastQC

class fastQCres:
    def __init__(self):
        self.fileLocation = ""
        self.fileName = "fastqc_data.txt"
        self.QCtables = {'Basic Statistics': [], 'Per base sequence quality': [], 'Per sequence quality scores': [],
                         'Per base sequence content': [], 'Per base GC content': [], 'Per sequence GC content': [],
                         'Per base N content': [], 'Sequence Length Distribution': [],
                         'Sequence Duplication Levels': [], 'Overrepresented sequences': [], 'Kmer Content': []}

    def lines2tableFormat(self, currentKey, tableLines):
        def baseNumerToStartEnd(baseNumber):
            bNbits = baseNumber.split("-")
            if len(bNbits) == 2:
                baseStart = bNbits[0]
                baseEnd = bNbits[1]
            else:
                baseStart = bNbits[0]
                baseEnd = bNbits[0]
            return ([baseStart, baseEnd])

        def addLines(tableLines):
            for li in range(1, len(tableLines)):
                tmpLine = []
                bits = tableLines[li].split()
                start = 0
                if tableLines[0][:5] == "#Base":
                    basePos = baseNumerToStartEnd(bits[0])
                    tmpLine.append(basePos[0])
                    tmpLine.append(basePos[1])
                    start = 1

                if tableLines[li][0] != "#":
                    for g in range(start, len(bits)):
                        tmpLine.append(bits[g])

                if tableLines[0][:27] == '#Total Duplicate Percentage':
                    tdpBits = tableLines[0].split()
                    tmpLine.append(tdpBits[-1])

                self.QCtables[currentKey].append(tmpLine[:])

        header = []
        if currentKey == "Per base sequence quality":
            header = ["BaseStart", ",BaseEnd", "Mean", "Median", "Lower_Quartile", "Upper_Quartile", "10th_Percentile",
                      "90th_Percentile"]
            self.QCtables[currentKey].append(header)
            addLines(tableLines)

        elif currentKey == "Per sequence GC content":
            header = ["GC_Content", "Count"]
            self.QCtables[currentKey].append(header)
            addLines(tableLines)

        elif currentKey == "Per base N content":
            header = ["BaseStart", ",BaseEnd", "N_Count"]
            self.QCtables[currentKey].append(header)
            addLines(tableLines)

        elif currentKey == "Sequence Duplication Levels":
            header = ["Duplication Level", "Relative count", "Total Duplicate Percentage"]
            self.QCtables[currentKey].append(header)
            addLines(tableLines)

        else:
            header = []
            if len(tableLines) > 0:
                he = tableLines[0][1:].split()
                if he[0] == "Base":
                    header.append("BaseStart")
                    header.append("BaseEnd")
                else:
                    header.append(he[0])

                for h in range(1, len(he)):
                    header.append(he[h])
            self.QCtables[currentKey].append(header)
            addLines(tableLines)

    def readInQCresult(self):
        import os

        fastQCfile = os.path.join(self.fileLocation, self.fileName)
        FQCF = open(fastQCfile, "r")
        lines = FQCF.readlines()
        FQCF.close()
        currentKey = ""
        tableLines = []
        for li in range(0, len(lines)):
            if lines[li][:2] == ">>" and lines[li][:12] != ">>END_MODULE":
                bits = lines[li][2:-1].split()
                currentKey = ""
                for g in range(0, len(bits) - 1):
                    currentKey = currentKey + bits[g] + " "
                currentKey = currentKey[:-1]
            elif lines[li][:12] == ">>END_MODULE":
                self.lines2tableFormat(currentKey, tableLines)
                tableLines = []
            else:
                tableLines.append(lines[li][:-1])
