import os
import sys


class sampleConfig:
    name = ""
    reads1 = ""
    reads2 = ""
    minQualMean = 0
    trimLeft = 0
    trimQualRight = 0
    trimQualWindow = 0
    trimQualStep = 0
    minLen = 0
    removeDup = ""
    performSubset = ""


confFile = sys.argv[1]
installationDirectory = sys.argv[2]

configFile = open(confFile)
projectName = ((configFile.readline().rstrip()).split('\t'))[1]




#Iterating the config file
while True:
    #Reading a sample config block
    sample = sampleConfig()
    header = configFile.readline()
    if not header:
        break

    sample.name = ((configFile.readline().rstrip()).split('\t'))[1]

    sample.reads1 = ((configFile.readline().rstrip()).split('\t'))[1]
    sample.reads2 = ((configFile.readline().rstrip()).split('\t'))[1]

    header = configFile.readline()

    sample.minQualMean = int(((configFile.readline().rstrip()).split('\t'))[1])

    sample.trimLeft = ((configFile.readline().rstrip()).split('\t'))[1]
    if not sample.trimLeft == "no":
        sample.trimLeft = int(sample.trimLeft)

    sample.trimQualRight = ((configFile.readline().rstrip()).split('\t'))[1]
    if not sample.trimQualRight == "no":
        sample.trimQualRight = int(sample.trimQualRight)

    sample.trimQualWindow = ((configFile.readline().rstrip()).split('\t'))[1]
    if not sample.trimQualWindow == "no":
        sample.trimQualWindow = int(sample.trimQualWindow)

    sample.trimQualStep = ((configFile.readline().rstrip()).split('\t'))[1]
    if not sample.trimQualStep == "no":
        sample.trimQualStep = int(sample.trimQualStep)

    sample.minLen = ((configFile.readline().rstrip()).split('\t'))[1]
    if not sample.minLen == "no":
        sample.minLen = int(sample.minLen)


    header = configFile.readline()
    

    #Qualiti filetering fastq 1
    print("Quality filtering for file ",sample.reads1)
    comand = installationDirectory+"src/conda/bin/prinseq-lite.pl -fastq " + sample.reads1 + " -out_format 3 -min_qual_mean " + str(sample.minQualMean)
    if not sample.trimLeft == "no":
        comand += " -trim_left " + str(sample.trimLeft)

    if not sample.trimQualRight == "no":
        comand += " -trim_qual_right " + str(sample.trimQualRight) + " -trim_qual_type mean "

    if not sample.trimQualWindow == "no":
        comand += " -trim_qual_window " + str(sample.trimQualWindow)

    if not sample.trimQualStep == "no":
        comand += " -trim_qual_step " + str(sample.trimQualStep)

    if not sample.minLen == "no":
        comand += " -min_len " + str(sample.minLen)

    comand += " -out_good " + sample.name+"_hq_1 "
    comand += " -out_bad badReads1 >" + sample.name+"_hq_1_filterStats.txt"


    os.system(comand)




    #Qualiti filetering fastq 2
    print("Quality filtering for file ",sample.reads2)
    comand = installationDirectory+"src/conda/bin/prinseq-lite.pl -fastq " + sample.reads2 + " -out_format 3 -min_qual_mean " + str(sample.minQualMean)
    if not sample.trimLeft == "no":
        comand += " -trim_left " + str(sample.trimLeft)

    if not sample.trimQualRight == "no":
        comand += " -trim_qual_right " + str(sample.trimQualRight) + " -trim_qual_type mean "

    if not sample.trimQualWindow == "no":
        comand += " -trim_qual_window " + str(sample.trimQualWindow)

    if not sample.trimQualStep == "no":
        comand += " -trim_qual_step " + str(sample.trimQualStep)

    if not sample.minLen == "no":
        comand += " -min_len " + str(sample.minLen)
    

    comand += " -out_good " + sample.name+"_hq_2 "
    comand += " -out_bad badReads2 > "+ sample.name+"_hq_2_filterStats.txt "

    os.system(comand)

    #Sort paired reads in files
    comand = installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts//assembly/utils/mergeQualFilteredMates.py " + sample.name+"_hq_1.fastq " + sample.name+"_hq_2.fastq"
    os.system(comand)
    







