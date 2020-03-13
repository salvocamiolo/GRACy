def reverseComplement(sequence):
    outSequence = ""
    accepted = ["A","C","T","G","a","c","t","g","N","n"]
    for a in range(len(sequence)-1,-1,-1):
        if sequence[a]=="A" or sequence[a]=="a":
            outSequence+="T"
        elif sequence[a]=="T" or sequence[a]=="t":
            outSequence+="A"
        elif sequence[a]=="G" or sequence[a]=="g":
            outSequence+="C"
        elif sequence[a]=="C" or sequence[a]=="c":
            outSequence+="G"
        else:
            outSequence+="N"
        

    return outSequence