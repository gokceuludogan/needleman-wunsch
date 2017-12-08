#!/usr/bin/envpython
import os,sys,re

penalty = 0
reward = 1
def getSequences(testSeq):
    textSeq = open(os.getcwd()+"/"+testSeq).read()
    sequences = re.findall(r'\w+', textSeq)
    print "The sequences are\n"
    print sequences , "\n"
    return sequences
def max(matchCost, delCost, inCost, mismatchCost):
    Max = matchCost
    index = 0 
    if(mismatchCost > Max):
        Max = mismatchCost
        index = 3
    if(delCost > Max):
        Max = delCost
        index = 1
    if(inCost > Max):
        Max = inCost
        index = 2
   
    return index     

def calculatePairwaiseAlignments(sequences):
    maxScore = 0
   
    print "Pairwise Alignments:"
    for i in range(0,len(sequences)):
        for j in range(i+1, len(sequences)):
            print "\n"
            (score, alignSeq1, alignSeq2) = globalAlignment(sequences[i], sequences[j])
            print "Alignments of sequence " , i+1 ," and sequence ", j+1 
            print alignSeq1
            print alignSeq2
            print "Score: " , score
            if score > maxScore:
                profileSeq1 = alignSeq1
                profileSeq2 = alignSeq2
                maxScore = score
                index1 = i
                index2 = j
    (profile, freqMatrix) = createProfileAndFreqMatrix(profileSeq1,profileSeq2)
  
  
    sequence = sequences[findSequence(index1,index2)]
    profileAlignment(profile, sequence, freqMatrix)

def globalAlignment(seq1, seq2):
   
    lenSeq1 = len(seq1)
    lenSeq2 = len(seq2)

    score = [[0 for i in range(lenSeq2+1)] for j in range(lenSeq1+1)] # represent score matrix
    
    backtrace = [[' ' for x in range(lenSeq2+1)] for y in range(lenSeq1+1)]  # contains information about which operation is chosen

    score[0][0] = 0
    #backtrace[0][0] = 'N'

    for i in range(lenSeq1+1):
    
        score[i][0] = penalty
        backtrace[i][0] = 'D'       # 'D' is used for delete operation
    
    for i in range(lenSeq2+1):
    
        score[0][i] = penalty
        backtrace[0][i] = 'I'       # 'I' represents insertion

    for i in range(1,lenSeq1+1):
        for j in range(1,lenSeq2+1):
            matchScore = 0
            if seq1[i-1] == seq2[j-1]:  
                matchScore = score[i-1][j-1] + reward
                backtrace[i][j] = 'M'  # 'M' represents match  
        
            delScore = score[i-1][j] + penalty
            inScore = score[i][j-1] + penalty
            mismatchScore = score[i-1][j-1] + penalty
            maxScore = max(matchScore, delScore, inScore, mismatchScore)
            if( maxScore == 1):
                score[i][j] = delScore
                backtrace[i][j] = 'D'
            elif( maxScore == 2):
                score[i][j] = inScore
                backtrace[i][j] = 'I'
            elif( maxScore == 3):
                score[i][j] = mismatchScore
                backtrace[i][j] = 'N'
            elif maxScore == 0:
                score[i][j] = matchScore
                backtrace[i][j] = 'M'
            else:
                print "Houston we have a problem about scoring!"      


    row = lenSeq1
    col = lenSeq2
    alSeq1 = ""
    alSeq2 = ""
    while row > -1 and col > -1:
        opCode = backtrace[row][col]
        if opCode == 'I':
            if col > 0: 
                alSeq2 = alSeq2 + seq2[col-1]
            else: 
                alSeq2 = alSeq2 + "-"
      
            alSeq1 = alSeq1 + "-"    
            col = col - 1
        elif opCode == 'D':

            if row > 0:
                alSeq1 = alSeq1 + seq1[row-1]
            else: 
                alSeq1 = alSeq1 + "-"
            alSeq2 = alSeq2 + "-"
            row = row - 1
        elif opCode == 'N':
            alSeq1 = alSeq1 + seq1[row-1]
            alSeq2 = alSeq2 + seq2[col-1]
            row = row - 1
            col = col - 1
        elif opCode == 'M':
            if row > 0:
                alSeq1 = alSeq1 + seq1[row-1]
            else: 
                alSeq1 = alSeq1 + "-"
            if col > 0:    
                alSeq2 = alSeq2 + seq2[col-1]
            else:
                alSeq2 = alSeq2 + "-"
            row = row - 1
            col = col - 1            
        else:
            print "Houston we have a problem about backtracing!"    
    
  
    alSeq1 = alSeq1[::-1]
    alSeq1 = alSeq1[1:]
    alSeq2 = alSeq2[::-1]
    alSeq2 = alSeq2[1:]
   
    return (score[lenSeq1][lenSeq2], alSeq1, alSeq2)

def profileAlignment(profile, sequence, freqMatrix):
  
    lenSeq = len(sequence)
    lenProfile = len(profile)
    score = [[0 for i in range(lenProfile+1)] for j in range(lenSeq+1)]    
    backtrace = [[' ' for x in range(lenProfile+1)] for y in range(lenSeq+1)]  # contains information about which operation is chosen

    score[0][0] = 0
    #backtrace[0][0] = 'N'
    penalty = 0
    for i in range(lenSeq+1):
    
        score[i][0] = penalty
        backtrace[i][0] = 'D'       # 'D' is used for delete operation
    
    for i in range(lenProfile+1):
    
        score[0][i] = penalty
        backtrace[0][i] = 'I'       # 'I' represents insertion
    for i in range(1,lenSeq+1):
        for j in range(1,lenProfile+1):

            matchScore = 0
            if sequence[i-1] == "A":
                matchScore = score[i-1][j-1] + freqMatrix[0][j-1]                 
            elif sequence[i-1] == "T":
                matchScore = score[i-1][j-1] + freqMatrix[1][j-1] 
            elif sequence[i-1] == "C":
                matchScore = score[i-1][j-1] + freqMatrix[2][j-1] 
            elif sequence[i-1] == "G":
                matchScore = score[i-1][j-1] + freqMatrix[3][j-1]
            else:
                matchScore = score[i-1][j-1] + freqMatrix[4][j-1]
        
            delScore = score[i-1][j] + penalty
            inScore = score[i][j-1] + penalty
            mismatchScore = score[i-1][j-1] + penalty
            maxScore = max(matchScore, delScore, inScore, mismatchScore)
            if maxScore == 0:
                score[i][j] = matchScore
                backtrace[i][j] = "M"
            elif( maxScore == 1):
                score[i][j] = delScore
                backtrace[i][j] = 'D'
            elif( maxScore == 2):
                score[i][j] = inScore
                backtrace[i][j] = 'I'
            elif( maxScore == 3):
                score[i][j] = mismatchScore
                backtrace[i][j] = 'N'
            else:
                print "Houston we have a problem about scoring in profile alignment!"  
  
  
    
    (alProfile, alSeq) = alignProfileWithSeq(sequence, profile, backtrace)
    createFinalProfile(alProfile, alSeq, score[lenSeq-1][lenProfile])
def createFinalProfile(alProfile, alSeq,score):
    profile = alProfile.split()
    alSeq1 = ""
    alSeq2 = ""

    print "\nThe alignment of profile with the sequence"
    print alProfile.strip()
    for index in range(0,len(profile)):
        if len(profile[index]) != 1:
            print "",alSeq[index],"",
        else:
            print alSeq[index],    
    for i in range(0, len(profile)):
        if len(profile[i]) == 1:
            alSeq1 = alSeq1 + profile[i]
            alSeq2 = alSeq2 + profile[i]

        else:
            alSeq1 = alSeq1 + profile[i][0].upper()
            alSeq2 = alSeq2 + profile[i][2].upper()



    print "\n"    
    print "Multiple Sequence Alignment\n"
    print "The max alignment score of the final multiple sequence:", score
    print "\nThe final multiple alignment:\n"
    print alSeq1
    print alSeq2
    print alSeq




def alignProfileWithSeq(sequence, profile, backtrace):
    row = len(sequence)
    col = len(profile)
    alSeq = ""
    alProfile = ""
    while row > -1 and col > -1:
        opCode = backtrace[row][col]
        if opCode  == "M":
            if row > 0:
                alSeq = sequence[row-1] + alSeq 

            else:
                alSeq = "-" + alSeq

            if col > 0:    
                alProfile = profile[col-1] + " " + alProfile  

            else:
                alProfile = "-" + " " + alProfile 
            row = row - 1
            col = col - 1   

        elif opCode == "I":
            if col > 0: 
                alProfile = profile[col-1] + " " + alProfile  
            else: 
                alProfile = "-" + " " + alProfile  
            alSeq = "-" +alSeq
            col = col - 1

        elif opCode == "D":
            if row > 0:
                alSeq = sequence[row-1] + alSeq 
            else: 
                alSeq = "-" + alSeq 
            alProfile = "-" + " " + alProfile 
            row = row - 1

        elif opCode == "N":
            alSeq = sequence[row-1] + alSeq 
            alProfile = profile[col-1] + " " + alProfile  
            row = row - 1
            col = col - 1

        else:
            print "Houston we have problem in creating final profile!"   

    alSeq = alSeq[1:]
    alProfile = alProfile[2:]
    return (alProfile, alSeq)
def findSequence(index1, index2):
    for i in range(3):
        if i != index1 and i != index2:
            return i
def createProfileAndFreqMatrix(seq1, seq2):
    freq = [[0 for i in range(len(seq1))] for j in range(5)] # represents the frequencies of nucleotides and gaps in the profile 
    # i.e. for AT-Ct/g
    #   A T - C t/g 
    # A 1 0 0 0 0
    # T 0 1 0 0 1/2
    # C 0 0 0 1 0
    # G 0 0 0 0 1/2
    # - 0 0 1 0 0

    profile = [' ' for i in range(len(seq1))]
    for i in range(0, len(seq1)):
        if seq1[i] == seq2[i]:

            profile[i] = seq1[i]
            if seq1[i] == "A":

                freq[0][i] = 1
            elif seq1[i] == "T":
                freq[1][i] = 1
            elif seq1[i] == "C":
                freq[2][i] = 1
            elif seq1[i] == "G":
                freq[3][i] = 1
            elif seq1[i] == "-":
                freq[4][i] = 1
            else:
                #print seq1[i]
                print "Houston the sequence contains some weird nucleotides!"            
        elif seq1[i] == "-":
            profile[i] =  "-/" +seq2[i].lower() 
            freq[4][i] = 1.0/2
            if seq2[i] == "A":
                freq[0][i] = 1.0/2
            elif seq2[i] == "T":
                freq[1][i] = 1.0/2
            elif seq2[i] == "C":
                freq[2][i] = 1.0/2
            elif seq2[i] == "G":
                freq[3][i] = 1.0/2
            else:
                
                print "Houston the sequence contains some weird nucleotides!"         

        elif seq2[i] == "-":
            profile[i] = seq1[i].lower() + "/-"
            freq[4][i] = 1.0/2
            if seq1[i] == "A":
                freq[0][i] = 1.0/2
            elif seq1[i] == "T":
                freq[1][i] = 1.0/2
            elif seq1[i] == "C":
                freq[2][i] = 1.0/2
            elif seq1[i] == "G":
                freq[3][i] = 1.0/2
            else:
                
                print "Houston the sequence contains some weird nucleotides!" 
        else:
            profile[i] = seq1[i].lower() +"/" + seq2[i].lower()
            if seq1[i] == "A":
                freq[0][i] = 1.0/2
            elif seq1[i] == "T":
                freq[1][i] = 1.0/2
            elif seq1[i] == "C":
                freq[2][i] = 1.0/2
            elif seq1[i] == "G":
                freq[3][i] = 1.0/2
            if seq2[i] == "A":
                freq[0][i] = 1.0/2
            elif seq2[i] == "T":
                freq[1][i] = 1.0/2
            elif seq2[i] == "C":
                freq[2][i] = 1.0/2
            elif seq2[i] == "G":
                freq[3][i] = 1.0/2  

               
    return (profile, freq)

sequences = []

if len(sys.argv) == 2:    
    testSeq = sys.argv[1]
    sequences = getSequences(testSeq)
    calculatePairwaiseAlignments(sequences)
elif len(sys.argv) == 4:
    seq1 = sys.argv[1]
    seq2 = sys.argv[2]
    seq3 = sys.argv[3]
    sequences.append(seq1)
    sequences.append(seq2)
    sequences.append(seq3)
    calculatePairwaiseAlignments(sequences)

else:
    print "usage: \n python aligment.py <path-of-input-file> \n or \n python alignment.py <sequence1> <sequence2> <sequence3>"  



