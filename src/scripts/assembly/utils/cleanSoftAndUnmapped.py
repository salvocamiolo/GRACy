import sys

filename = sys.argv[1]

infile = open(filename)
outfile = open(filename+"_cleaned.sam","w")


line1 = infile.readline().rstrip()
if not line1:
    print("cleanSoftAndUnmapped failed. Now exit")
    exit()
outfile.write(line1+"\n")


while line1[0]=="@":
    line1 = infile.readline().rstrip()
    outfile.write(line1+"\n")
    if not line1:
        break


while True:
    line2 = infile.readline().rstrip()
    if not line2:
        break
    
    fields1 = line1.split("\t")
    fields2 = line2.split("\t")

    if not fields1[0] == fields2[0]:
        print("Not a pair!!")
        exit()
    else:
        if not "S" in fields1[5] and not "S" in fields2[5] and  not "*" in fields1[5] and not "*" in fields2[5]:
            outfile.write(line1+"\n"+line2+"\n")
    line1 = infile.readline().rstrip()


outfile.close()






