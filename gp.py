'''
Created on Apr 5, 2011

@author: Ross Eldridge

Python 2.x.x
syntax: gp.py [test sequence FASTA file] [training sequence FASTA file]

The test sequence FASTA file should be a single entry, but the training
FASTA file may be a multiple entry FASTA file.

Default Consensus Sequence probabilities obtained from:
C.B.Harley, R.P.Reynolds, Nucleic Acid Research, Vol. 15, No. 5, 1987, pp. 2343-2361

Ideas to implement for the future:
- UP element at -50 (occasional replacement for -35 sequence)
- Possibility of multiple promoter sequences in the promoter region up-regulating
- Alternative start codons (GTG, TTG)
- Alternative promoters for other sigma factors (all here are for s70 sigma factor)
- Implement circular DNA support
- Adjust training scheme - instead of training based on a default pribnow/-35 sequence,
    attempt to find unique consensus sequences using multiple sequence alignment on
    large coding regions
'''

import sys
import getopt
import math
from time import clock

def expandDictionary(s,deb):
    """ Expands an ACTG scoring dictionary to
        include values for all ambiguous nucleotides
        as well.
        
        deb = debugging flag
    """
    if s["A"]+s["C"]+s["G"]+s["T"] != 1 and deb:
        print "Warning! Scoring dictionary unbalanced.\n"
        print s["A"],s["C"],s["G"],s["T"],"\n"
    s["R"] = (s["A"]+s["G"])/2
    s["Y"] = (s["C"]+s["T"])/2
    s["S"] = (s["G"]+s["C"])/2
    s["W"] = (s["A"]+s["T"])/2
    s["K"] = (s["G"]+s["T"])/2
    s["M"] = (s["A"]+s["C"])/2
    s["B"] = (s["C"]+s["G"]+s["T"])/3
    s["D"] = (s["A"]+s["G"]+s["T"])/3
    s["H"] = (s["A"]+s["C"]+s["T"])/3
    s["V"] = (s["A"]+s["C"]+s["G"])/3
    s["N"] = (s["A"]+s["C"]+s["G"]+s["T"])/4
    return s

def readFasta(filename):
    """ Reads a fasta file by lines, removing
        all headers and leaving only the sequence.
        Currently set up such that given a multi-
        sequence fasta file, the sequences will be appended
        together to form a single sequence.
    """
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    for k in lines:
        if k[0] == '>':
            lines.remove(k)
    line = ''.join(lines).translate(None,'\n')
    return line.upper()

def findbyORF(seq,str):
    """ Returns a list, each of 6 entries being a list of indices of 
        string str in (very large) string seq.  These indices
        are sorted into the 6 open reading frames (3 forward and
        three backward), making up the return list.
        
        NOTE: The indices for the backwards ORF return the "end index,"
        which would theoretically be the beginning index if the strand
        were reversed.
        
        For example, searching for:
            ATG   in    CCATGCCCGTACC
                        -->>>---<<<--
        would return an index of 2 in ORFs[2] (forward #3), and an index
        of 10 (NOT 8) in ORFs[5] (backward #3).
    """
    l = len(str)
    str2 = str[::-1]
    f = []
    r = []
    i = 0
    temp = seq.find(str,i)
    while temp != -1:
        f.append(temp)
        i = temp+1
        temp = seq.find(str,i)
    i = 0
    temp = seq.find(str2,i)
    while temp != -1:
        r.append(temp+(l-1))
        i = temp+1
        temp = seq.find(str2,i)
    ORFs = []
    for i in range(6):
        ORFs.append([])
    for k in f:
        if k%3 == 0:
            ORFs[0].append(k)
        elif k%3 == 1:
            ORFs[1].append(k)
        elif k%3 == 2:
            ORFs[2].append(k)
    for k in r:
        if (l-k-1)%3 == 0:
            ORFs[3].append(k)
        elif (l-k-1)%3 == 1:
            ORFs[4].append(k)
        elif (l-k-1)%3 == 2:
            ORFs[5].append(k)
    return ORFs

def alignSeqs(s1,s2):
    """ Aligns two sequences with no gaps.  Dynamic programming
        not required, and a simple loop is sufficient.  Function
        returns a list of indexes relative to s2 of maximum scores.
        The shorter sequence should always be s1.
        
        Example: given s1=ACCTGG, s2 = TTAGACCGTGCCTG
              Alignment           Function would return maxindex = [4]
            ----ACCTGG----
            TTAGACCGTGCCTG
    """
    scores = []
    for i in range(len(s2)-len(s1)):
        score = 0
        for j in range(len(s1)):
            score += match(s1[j],s2[i+j])
        scores.append(score)
    maxscore = max(scores)
    maxindex = []
    for i,j in enumerate(scores):
        if j == maxscore:
            maxindex.append(i)
    return maxindex

def match(c1,c2):
    """ Returns the value of matching two nucleotides,
        taking into account ambiguous nucleotides.  For example,
        matching an 'A' to an 'N' would be a 1 in 4 chance, but
        so would matching an 'R' ('A' or 'G') to an 'S' ('G' or 'C').
        Value returned is based proportionally on the chance of a match
        and a perfect match score of 12.
    """
    actg = {"A","C","T","G"}
    if(c1 == c2):
        return 12 # perfect match = score of 12
    elif(c1 in actg and c2 in actg):
        return 0 # Simple check to exit the function early
    matchups = {"A":["R","W","M"],"C":["Y","S","M"],
                 "T":["Y","W","K"],"G":["R","S","K"],
                 "R":[],"Y":[],"S":[],"W":[],"K":[],"M":[],
                 "B":[],"D":[],"H":[],"V":[],"N":[]}
    if(c2 in matchups[c1] or c1 in matchups[c2]):
        return 6 # 50% chance of match = score of 6
    matchups = {"A":["D","H","V"],"C":["B","H","V"],
                 "T":["B","D","H"],"G":["B","D","V"],
                 "R":[],"Y":[],"S":[],"W":[],"K":[],"M":[],
                 "B":[],"D":[],"H":[],"V":[],"N":[]}
    if(c2 in matchups[c1] or c1 in matchups[c2]):
        return 4 # 33.3% chance of match = score of 4
    matchups = {"A":["N"],"C":["N"],"T":["N"],"G":["N"],
                 "R":["S","W","K","M"],"Y":["S","W","K","M"],
                 "S":["R","Y","K","M"],"W":["R","Y","K","M"],
                 "K":["R","Y","S","W"],"M":["R","Y","S","W"],
                 "B":[],"D":[],"H":[],"V":[],"N":[]}
    if(c2 in matchups[c1] or c1 in matchups[c2]):
        return 3 # 25% chance of match = score of 3
    return 0 # Anything below 25% chance of match

def outputORFs(pORF,seq,outg,outt):
    """ Outputs the putative sequences to a gene prediction file,
        as well as the translated proteins to a separate file.
    """
    if outg == "":
        outg = "GenePrediction.txt"
    if outt == "":
        outt = "PredictedProteins.fasta"
    with open(outg,'w') as f:
        print "Outputting",outg
        with open(outt,'w') as g:
            print "Outputting",outt
            f.write("Start".rjust(13))
            f.write("Stop".rjust(13))
            f.write("LenNuc".rjust(8))
            f.write("LenAA".rjust(8))
            f.write("+-".rjust(3))
            f.write("PribBox".rjust(11))
            f.write("PribScore".rjust(10))
            f.write("-35Seq".rjust(11))
            f.write("-35Score".rjust(10))
            f.write("\n\n")
            for i in pORF:
                # GenePrediction.txt
                #f.write("%s %s %s %s \n" % (i.start,i.stop,i.len,i.dir))
                f.write(str(i.start).rjust(13))
                f.write(str(i.stop).rjust(13))
                f.write(str(i.len).rjust(8))
                f.write(str(i.len/3).rjust(8))
                f.write(i.dir.rjust(3))
                f.write(i.pribnow.rjust(11))
                f.write("{0:.5f}".format(i.pribscore).rjust(10))
                f.write(i.seq35.rjust(11))
                f.write("{0:.5f}".format(i.seq35score).rjust(10))
                f.write("\n")
                # PredictedProteins.txt
                g.write(">CR|%d|%d|%dbp|%daa\n" % \
                        (i.start,i.stop,i.len,(i.len)/3))
                if i.dir == "+":
                    g.write(translate(seq[i.start:i.stop]))
                elif i.dir == "-":
                    g.write(translate(seq[i.start:i.stop:-1]))
                else:
                    g.write("error - invalid direction")
                g.write("\n")
    #f.close()
    #g.close()

def translate(s):
    """ Translates a DNA sequence into a protein sequence.
        Any ambiguous codons will be translated as "X".
    """
    trans = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L",
             "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
             "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
             "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
             "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
             "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
             "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
             "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
             "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
             "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
             "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
             "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
             "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
             "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
             "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
             "GGT":"G","GGC":"G","GGA":"G","GGG":"G"}
    p = ""
    if(len(s) % 3 != 0):
        print "Error: sequence to be translated not in same ORF\n"
        return "eRR"
    iter = 0
    while iter < len(s):
        if (s[iter:iter+3] in trans):
            p += trans[s[iter:iter+3]]
        else:
            p+="X"
        iter += 3
        if iter % 210 == 0 and iter < len(s):
            p += '\n'
    return p

class ORF:
    """Open Reading Frame class
        Keeps track of information for a single ORF.  Does *NOT* keep
        actual sequence in memory - only absolute start/stop locations.
    """
    def __init__(self, start, stop, prib, pribscore, s35, s35score):
        self.start = start
        self.stop = stop
        self.len = abs(stop - start)
        if start < stop:
            self.dir = '+'
        else:
            self.dir = '-'
        self.pribnow = prib
        self.pribscore = pribscore
        self.seq35 = s35
        self.seq35score = s35score

def train(trainfile,pribnow,seq35):
    """
    Uses a tighter set of constraints to attempt to train the Pribnow table, -35 table, and Codon counts from another dataset.
    """
    # settings for training set are more rigid
    settings = {"min_size":250,"max_size":6000,"pcut":-2.75,"tcut":-2.75}
    seq = readFasta(trainfile)
    startORF = findbyORF(seq,'ATG')
    amberORF = findbyORF(seq,'TAG') #Different STOP codons
    ochreORF = findbyORF(seq,'TAA')
    umberORF = findbyORF(seq,'TGA')
    stopORF = []
    for i in range(6):
        stopORF.append(amberORF[i] + ochreORF[i] + umberORF[i])
        stopORF[i].sort()
    pORF = findCRs(seq,startORF,stopORF,settings,pribnow,seq35,"")
    print "Script being trained on",len(pORF),"high probability sequences."
    # pribnow table
    p = [{"A":0.0,"C":0.0,"G":0.0,"T":0.0},
            {"A":0.0,"C":0.0,"G":0.0,"T":0.0},
            {"A":0.0,"C":0.0,"G":0.0,"T":0.0},
            {"A":0.0,"C":0.0,"G":0.0,"T":0.0},
            {"A":0.0,"C":0.0,"G":0.0,"T":0.0},
            {"A":0.0,"C":0.0,"G":0.0,"T":0.0}]
    # -35 sequence table
    s = [{"A":0.0,"C":0.0,"G":0.0,"T":0.0},
            {"A":0.0,"C":0.0,"G":0.0,"T":0.0},
            {"A":0.0,"C":0.0,"G":0.0,"T":0.0},
            {"A":0.0,"C":0.0,"G":0.0,"T":0.0},
            {"A":0.0,"C":0.0,"G":0.0,"T":0.0},
            {"A":0.0,"C":0.0,"G":0.0,"T":0.0}]
    # Codon counts
    cc = {"TTT":0.0,"TTC":0.0,"TTA":0.0,"TTG":0.0,
         "TCT":0.0,"TCC":0.0,"TCA":0.0,"TCG":0.0,
         "TAT":0.0,"TAC":0.0,"TAA":0.0,"TAG":0.0,
         "TGT":0.0,"TGC":0.0,"TGA":0.0,"TGG":0.0,
         "CTT":0.0,"CTC":0.0,"CTA":0.0,"CTG":0.0,
         "CCT":0.0,"CCC":0.0,"CCA":0.0,"CCG":0.0,
         "CAT":0.0,"CAC":0.0,"CAA":0.0,"CAG":0.0,
         "CGT":0.0,"CGC":0.0,"CGA":0.0,"CGG":0.0,
         "ATT":0.0,"ATC":0.0,"ATA":0.0,"ATG":0.0,
         "ACT":0.0,"ACC":0.0,"ACA":0.0,"ACG":0.0,
         "AAT":0.0,"AAC":0.0,"AAA":0.0,"AAG":0.0,
         "AGT":0.0,"AGC":0.0,"AGA":0.0,"AGG":0.0,
         "GTT":0.0,"GTC":0.0,"GTA":0.0,"GTG":0.0,
         "GCT":0.0,"GCC":0.0,"GCA":0.0,"GCG":0.0,
         "GAT":0.0,"GAC":0.0,"GAA":0.0,"GAG":0.0,
         "GGT":0.0,"GGC":0.0,"GGA":0.0,"GGG":0.0}
    for i in p:
        i = expandDictionary(i,0)
    for i in s:
        i = expandDictionary(i,0)
    for i in pORF:
        for j in range(6):
            p[j][i.pribnow[j]] += 1
            s[j][i.seq35[j]] += 1
        codon = ""
        # Codon bias analysis
        if(i.start < i.stop):
            orf = seq[i.start:i.stop+2]
        else:
            orf = seq[i.start:i.stop-2:-1]
        for k in orf:
            codon += k
            if len(codon) == 3:
                if codon in cc:
                    cc[codon] += 1
                codon = ""
    cc = codonTable(cc)
    for i in p:
        sum = i["A"]+i["C"]+i["G"]+i["T"]
        for j in i:
            i[j] /= sum
        i = expandDictionary(i,1)
    for i in s:
        sum = i["A"]+i["C"]+i["G"]+i["T"]
        for j in i:
            i[j] /= sum
        i = expandDictionary(i,1)
    return (p,s,cc)

def codonTable(cc):
    """
    Takes a table of codon sums and divides by the sum per
    amino acid to make a probability table of codon bias.
    
    Stop codons included as well.
    """
    Phe = max(cc["TTT"] + cc["TTC"],1)
    Leu = max(cc["TTA"] + cc["TTG"] + cc["CTT"] + cc["CTC"] + cc["CTA"] + cc["CTG"],1)
    Ser = max(cc["TCT"] + cc["TCC"] + cc["TCA"] + cc["TCG"] + cc["AGT"] + cc["AGC"],1)
    Tyr = max(cc["TAT"] + cc["TAC"],1)
    Stop = max(cc["TAA"] + cc["TAG"] + cc["TGA"],1)
    Cys = max(cc["TGT"] + cc["TGC"],1)
    Trp = max(cc["TGG"],1)
    Pro = max(cc["CCT"] + cc["CCC"] + cc["CCA"] + cc["CCG"],1)
    His = max(cc["CAT"] + cc["CAC"],1)
    Gln = max(cc["CAA"] + cc["CAG"],1)
    Arg = max(cc["CGT"] + cc["CGC"] + cc["CGA"] + cc["CGG"] + cc["AGA"] + cc["AGG"],1)
    Ile = max(cc["ATT"] + cc["ATC"] + cc["ATA"],1)
    Met = max(cc["ATG"],1)
    Thr = max(cc["ACT"] + cc["ACC"] + cc["ACA"] + cc["ACG"],1)
    Asn = max(cc["AAT"] + cc["AAC"],1)
    Lys = max(cc["AAA"] + cc["AAG"],1)
    Val = max(cc["GTT"] + cc["GTC"] + cc["GTA"] + cc["GTG"],1)
    Ala = max(cc["GCT"] + cc["GCC"] + cc["GCA"] + cc["GCG"],1)
    Asp = max(cc["GAT"] + cc["GAC"],1)
    Glu = max(cc["GAA"] + cc["GAG"],1)
    Gly = max(cc["GGT"] + cc["GGC"] + cc["GGA"] + cc["GGG"],1)
    cc["TTT"] /= Phe
    cc["TTC"] /= Phe
    cc["TTA"] /= Leu
    cc["TTG"] /= Leu
    cc["TCT"] /= Ser
    cc["TCC"] /= Ser
    cc["TCA"] /= Ser
    cc["TCG"] /= Ser
    cc["TAT"] /= Tyr
    cc["TAC"] /= Tyr
    cc["TAA"] /= Stop
    cc["TAG"] /= Stop
    cc["TGT"] /= Cys
    cc["TGC"] /= Cys
    cc["TGA"] /= Stop
    cc["TGG"] /= Trp
    cc["CTT"] /= Leu
    cc["CTC"] /= Leu
    cc["CTA"] /= Leu
    cc["CTG"] /= Leu
    cc["CCT"] /= Pro
    cc["CCC"] /= Pro
    cc["CCA"] /= Pro
    cc["CCG"] /= Pro
    cc["CAT"] /= His
    cc["CAC"] /= His
    cc["CAA"] /= Gln
    cc["CAG"] /= Gln
    cc["CGT"] /= Arg
    cc["CGC"] /= Arg
    cc["CGA"] /= Arg
    cc["CGG"] /= Arg
    cc["ATT"] /= Ile
    cc["ATC"] /= Ile
    cc["ATA"] /= Ile
    cc["ATG"] /= Met
    cc["ACT"] /= Thr
    cc["ACC"] /= Thr
    cc["ACA"] /= Thr
    cc["ACG"] /= Thr
    cc["AAT"] /= Asn
    cc["AAC"] /= Asn
    cc["AAA"] /= Lys
    cc["AAG"] /= Lys
    cc["AGT"] /= Ser
    cc["AGC"] /= Ser
    cc["AGA"] /= Arg
    cc["AGG"] /= Arg
    cc["GTT"] /= Val
    cc["GTC"] /= Val
    cc["GTA"] /= Val
    cc["GTG"] /= Val
    cc["GCT"] /= Ala
    cc["GCC"] /= Ala
    cc["GCA"] /= Ala
    cc["GCG"] /= Ala
    cc["GAT"] /= Asp
    cc["GAC"] /= Asp
    cc["GAA"] /= Glu
    cc["GGG"] /= Glu
    cc["GGT"] /= Gly
    cc["GGC"] /= Gly
    cc["GGA"] /= Gly
    cc["GGG"] /= Gly
    return cc

def findCRs(seq,startORF,stopORF,settings,pribnow,seq35,codonbias):
    """
    Finds plausible coding regions in FASTA sequence "seq" using a list of start and stop codon locations and tables of biases.
    """
    pribtemplate = consensus(pribnow)
    seq35template = consensus(seq35)
    pORF = []
    for i in (0,1,2): #forward ORFs
        curstart = 0
        curstop = 0
        prevstop = 0
        while(curstart < len(startORF[i]) and curstop < len(stopORF[i])):
            j = startORF[i][curstart]
            k = stopORF[i][curstop]
            if j < prevstop:
                curstart += 1
            else:
                # Start codon is too close to or even past stop codon
                if k - j < settings["min_size"]:
                    curstop += 1
                    prevstop = k
                # Start codon is too far away from stop codon
                elif k - j > settings["max_size"]:
                    curstart += 1
                else:
                    pscore = 1
                    tscore = 1
                    if j > 15:
                        ind = alignSeqs(pribtemplate,seq[j-15:j-5])
                        prib = seq[j-15+ind[0]]+seq[j-14+ind[0]]+\
                            seq[j-13+ind[0]]+seq[j-12+ind[0]]+\
                            seq[j-11+ind[0]]+seq[j-10+ind[0]]
                    else:
                        prib = ""
                        pscore = 0.000001
                    if j > 40:
                        ind = alignSeqs(seq35template,seq[j-40:j-30])
                        s35 = seq[j-40+ind[0]]+seq[j-39+ind[0]]+\
                            seq[j-38+ind[0]]+seq[j-37+ind[0]]+\
                            seq[j-36+ind[0]]+seq[j-35+ind[0]]
                    else:
                        s35 = ""
                        tscore = 0.000001
                    for m in range(len(prib)):
                        pscore *= pribnow[m][prib[m]]
                    for m in range(len(s35)):
                        tscore *= seq35[m][s35[m]]
                    pscore = math.log10(pscore)
                    tscore = math.log10(tscore)
                    cscore = codonScore(seq[j:k+2],codonbias)
                    if pscore + tscore + cscore > settings["pcut"] + settings["tcut"]:
                        pORF.append(ORF(j,k,prib,pscore,s35,tscore))
                        """print "start",j,"stop,",k,"numgene",len(pORF),"FoR","+"+str(i+1),\
                            "Pribnow",pORF[len(pORF)-1].pribnow,"PribScore",\
                            pORF[len(pORF)-1].pribscore,prib,seq[j-15:j-5],"CodonScore",cscore
                            debugging"""
                    curstart += 1
    for i in (3,4,5): #backward ORFs
        curstart = len(startORF[i])-1
        curstop = len(stopORF[i])-1
        prevstop = len(seq)
        while(curstart >= 0 and curstop >= 0):
            j = startORF[i][curstart]
            k = stopORF[i][curstop]
            if j > prevstop:
                curstart -= 1
            else:
                # Start codon is too close to or even past stop codon
                if j - k < settings["min_size"]:
                    curstop -= 1
                    prevstop = k
                elif j - k > settings["max_size"]:
                    curstart -= 1
                else:
                    pscore = 1
                    tscore = 1
                    if j < len(seq)-15:
                        ind = alignSeqs(pribtemplate[::-1],seq[j+15:j+5:-1])
                        prib = seq[j+15-ind[0]]+seq[j+14-ind[0]]+seq[j+13-ind[0]]+\
                        seq[j+12-ind[0]]+seq[j+11-ind[0]]+seq[j+10-ind[0]]
                    else:
                        prib = ""
                        pscore = 0.000001
                    if j < len(seq)-40:
                        ind = alignSeqs(pribtemplate[::-1],seq[j+40:j+30:-1])
                        s35 = seq[j+40-ind[0]]+seq[j+39-ind[0]]+seq[j+38-ind[0]]+\
                        seq[j+37-ind[0]]+seq[j+36-ind[0]]+seq[j+35-ind[0]]
                    else:
                        s35 = ""
                        tscore = 0.000001
                    for m in range(len(prib)):
                        pscore *= pribnow[m][prib[m]]
                    for m in range(len(s35)):
                        tscore *= seq35[m][s35[m]]
                    pscore = math.log10(pscore)
                    tscore = math.log10(tscore)
                    cscore = codonScore(seq[j:k-2:-1],codonbias)
                    if pscore + tscore + cscore > settings["pcut"] + settings["tcut"]:
                        pORF.append(ORF(j,k,prib,pscore,s35,tscore))
                        """print "start",j,"stop,",k,"numgene",len(pORF),"FoR","-"+str(i-2),\
                            "Pribnow",pORF[len(pORF)-1].pribnow,"PribScore",\
                            pORF[len(pORF)-1].pribscore,prib,seq[j+15:j+5:-1],"CodonScore",cscore
                            debugging"""
                    curstart -= 1
    return pORF

def codonScore(seq,ct):
    if ct == "":
        return 0
    if len(seq) < 3:
        return -10
    codon = ""
    score = 1.0
    for i in seq:
        codon += i
        if len(codon) == 3:
            if codon in ct:
                score += math.log(max(ct[codon],0.01))
            codon = ""
    return (score / (len(seq)/3)) + 0.5

def consensus(table):
    str = ""
    for i in table:
        l = [i["A"],i["C"],i["G"],i["T"]]
        m = max(l)
        d = {l[0]:"A",l[1]:"C",l[2]:"G",l[3]:"T"}
        str += d[m]
    return str

def usage():
    print "Usage: gp.py -h -i [in file] -t [in file] -g [out file] -p [out file]"
    print "   Informational:"
    print "      -h, --help : outputs this usage information"
    print "   Required:"
    print "      -i, --input : filename of test set file in FASTA format"
    print "   Optional:"
    print "      -t, --training : filename of training set file in FASTA format"
    print "                default = script has default settings if not present"
    print "      -g, --genes : designated filename for gene prediction output"
    print "                default = 'GenePrediction.txt'"
    print "      -p, --prot : designed filename for translated output"
    print "                default = 'PredictedProteins.fasta'"
    
def main(argv):
    clock()
    testfile = ""
    trainfile= ""
    outfileg = ""
    outfilet = ""
    try:
        opts, args = getopt.getopt(argv, "hi:t:g:p:", \
                                   ["help", "input=", "training=","genes=","prot="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--input"):
            testfile = arg
        elif opt in ("-t", "--training"):
            trainfile = arg
        elif opt in ("-g", "--genes"):
            outfileg = arg
        elif opt in ("-p", "--prot"):
            outfilet = arg
    if len(args) > 0:
        print "Unused arguments:",args,"\n"
    settings = {"min_size":150,"max_size":6000,"pcut":-3.3,"tcut":-3.3}
    # Default Pribnow sequence has an AT bias.
    pribnow = [{"A":0.0300,"C":0.0800,"G":0.0700,"T":0.8200},
                {"A":0.8900,"C":0.0300,"G":0.0100,"T":0.0700},
                {"A":0.2600,"C":0.1000,"G":0.1200,"T":0.5200},
                {"A":0.5900,"C":0.1200,"G":0.1500,"T":0.1400},
                {"A":0.4900,"C":0.2100,"G":0.1100,"T":0.1900},
                {"A":0.0325,"C":0.0525,"G":0.0225,"T":0.8925}]
    # Default -35 sequence
    seq35 = [{"A":0.0300,"C":0.0900,"G":0.1000,"T":0.7800},
                {"A":0.1000,"C":0.0300,"G":0.0500,"T":0.8200},
                {"A":0.0300,"C":0.1400,"G":0.6800,"T":0.1500},
                {"A":0.5775,"C":0.1275,"G":0.0975,"T":0.1975},
                {"A":0.3175,"C":0.5175,"G":0.0675,"T":0.0975},
                {"A":0.5400,"C":0.0500,"G":0.1700,"T":0.2400}]
    codonbias = ""
    if trainfile != "":
        pribnow, seq35, codonbias = train(trainfile,pribnow,seq35)
    for i in pribnow:
        i = expandDictionary(i,1)
    for i in seq35:
        i = expandDictionary(i,1)
    testseq = readFasta(testfile)
    startORF = findbyORF(testseq,'ATG')
    amberORF = findbyORF(testseq,'TAG')
    ochreORF = findbyORF(testseq,'TAA')
    umberORF = findbyORF(testseq,'TGA')
    stopORF = []
    for i in range(6):
        stopORF.append(amberORF[i] + ochreORF[i] + umberORF[i])
        stopORF[i].sort()
    pORF = findCRs(testseq,startORF,stopORF,settings,pribnow,seq35,codonbias)
    outputORFs(pORF,testseq,outfileg,outfilet)
    end_scr = clock()
    print "Total of",len(pORF),"putative coding regions found."
    print "Script completed in",end_scr,"seconds.\n"

if __name__ == '__main__':
    main(sys.argv[1:])