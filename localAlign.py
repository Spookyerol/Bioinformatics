# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 13:38:50 2019

@author: Erdal Guclu
"""
import numpy as np


def initializeMatrices(lenA, lenB):
    matrix = [[None for i in range(lenB+2)] for j in range(lenA+2)]
    backtr = [[None for i in range(lenB+2)] for j in range(lenA+2)]
    
    #initialise empty unused cells
    matrix[0][0] = " "
    matrix[0][1] = " "
    matrix[1][0] = " "
    matrix[1][1] = 0
    
    backtr[0][0]= ""
    backtr[0][1]= "  "
    backtr[1][0]= ""
    
    for i in range(2, lenA+2): #initialise indices of seqA
        matrix[i][0] = i-2
        matrix[i][1] = 0
        backtr[i][0] = i-2
        backtr[i][1] = "S"
    for j in range(2, lenB+2): #initialise indices of seqB
        matrix[0][j] = j-2
        matrix[1][j] = 0
        backtr[0][j] = j-2
        backtr[1][j] = "S"
    
    backtr[1][1] = "FN"
        
    return [matrix, backtr]

def checkInput(alphabet, scoringMatrix, A, B):
    for symbol in A:
        if(symbol not in alphabet):
            print("Error: Alphabet does not cover symbols appearing in seqA")
            return -1
    for symbol in B:
        if(symbol not in alphabet):
            print("Error: Alphabet does not cover symbols appearing in seqB")
            return -1
    if(len(scoringMatrix[0]) != len(alphabet)+1):
        print("Error: Invalid scoring matrix for this alphabet")
        return -1
    
    return 1


def getMatch(alphabet, scoringMatrix, A, B, i, j):
    symbolA = ""
    symbolB = ""
    for k in range(0, len(alphabet)):
        if(A[i] == alphabet[k]):
            symbolA = k
        if(B[j] == alphabet[k]):
            symbolB = k
            
    return scoringMatrix[symbolA][symbolB]
    


def findIndices(backtr, lenA, lenB, i, j):
    m = i + 2
    n = j + 2
    matchIndexA = []
    matchIndexB = []
    while(backtr[m][n] != "FN"):
        if backtr[m][n] == "L":
            n -= 1
        elif backtr[m][n] == "U":
            m -= 1
        elif backtr[m][n] == "D":
            matchIndexA.insert(0,m-2)
            matchIndexB.insert(0,n-2)
            m -= 1
            n -= 1
        elif(backtr[m][n] == "S"): #end of local alignment
            break
        
    return matchIndexA, matchIndexB



def dynprog(alphabet, scoringMatrix, seqA, seqB):
    if(checkInput(alphabet, scoringMatrix, seqA, seqB) == -1):
        return -1
    lenA = len(seqA)
    lenB = len(seqB)
    matrices = initializeMatrices(lenA, lenB)
    bestScore = [0, (0,0)] #stores best local score and the endpoint symbols
    for i in range(2, lenA+2):
        for j in range(2, lenB+2):
            #diagonal, up, left and fresh start respectively
            matchScore = getMatch(alphabet, scoringMatrix, seqA, seqB, i-2, j-2)
            score = max(matchScore + matrices[0][i-1][j-1], matrices[0][i-1][j] + scoringMatrix[0][len(scoringMatrix[0])-1], matrices[0][i][j-1] + scoringMatrix[len(scoringMatrix[0])-1][0], 0)

            if(score == 0):
                matrices[1][i][j] = "S"
            elif(score == matchScore + matrices[0][i-1][j-1]):
                matrices[1][i][j] = "D"
            elif(score == matrices[0][i-1][j] + scoringMatrix[0][len(scoringMatrix[0])-1]):
                matrices[1][i][j] = "U"
            elif(score == matrices[0][i][j-1] + scoringMatrix[len(scoringMatrix[0])-1][0]):
                matrices[1][i][j] = "L"
            
            matrices[0][i][j] = score
            if(score > bestScore[0]):
                bestScore[0] = score
                bestScore[1] = (i-2,j-2)
    
    indices = findIndices(matrices[1], lenA, lenB, bestScore[1][0], bestScore[1][1])
    print(np.matrix(matrices[0]))
    #print(np.matrix(matrices[1]))
    print("Best Local Score is:  ", bestScore[0])
    print("Resulting indices:  ", indices[0],indices[1])
    return bestScore[0],indices[0],indices[1]

#dynprog("ACTG", [[1,-1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-1,-2],[-1,-1,-1,1,-2],[-2,-2,-2,-2,-2]], "AAAC", "AGC")
dynprog("ACTG", [[1,-1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-1,-2],[-1,-1,-1,1,-2],[-2,-2,-2,-2,-2]], "TAATA", "TACTAA")
dynprog("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAACA", "CBACCCBA")


def initializeMatricesLin(lenA, lenB):
    if(lenA <= lenB):
        matrix = [[None for i in range(2)] for j in range(lenA+1)]
        backtr = [[None for i in range(2)] for j in range(lenA+1)]
    else:
        matrix = [[None for i in range(lenB+1)] for j in range(2)]
        backtr = [[None for i in range(lenB+1)] for j in range(2)]
        
    for i in range(0, len(matrix[0])): #initialise indices of seqA
        matrix[0][i] = 0
        backtr[0][i] = "S"
    for j in range(0, len(matrix)): #initialise indices of seqB
        matrix[j][0] = 0
        backtr[j][0] = "S"
        
    matrix[0][0] = 0
    backtr[0][0] = "FN"
        
    #print(np.matrix(matrix))
    #print(np.matrix(backtr))
    return [matrix, backtr]

#initializeMatricesLin(6,9)

def align(alphabet, scoringMatrix, seqA, seqB):
    lenA = len(seqA)
    lenB = len(seqB)
    
    matrices = initializeMatricesLin(lenA, lenB)
    bestScore = [0, (0,0)] #stores best local score and the endpoint symbols
    if(lenA <= lenB):
        Ashort = True
    else:
        Ashort = False
    
    short = min(lenA, lenB)
    count = 0
    while(count < short):
        for i in range(1, len(matrices[0])):
            for j in range(1, len(matrices[0][0])):
                #diagonal, up, left and fresh start respectively
                #print(seqA[i-1], seqB[j-1+count])
                if(Ashort):
                    matchScore = getMatch(alphabet, scoringMatrix, seqA, seqB, i-1, j-1+count)
                else:
                    matchScore = getMatch(alphabet, scoringMatrix, seqA, seqB, i-1+count, j-1)
                    
                score = max(matchScore + matrices[0][i-1][j-1], matrices[0][i-1][j] + scoringMatrix[0][len(scoringMatrix[0])-1], matrices[0][i][j-1] + scoringMatrix[len(scoringMatrix[0])-1][0], 0)
    
                if(score == 0):
                    matrices[1][i][j] = "S"
                elif(score == matchScore + matrices[0][i-1][j-1]):
                    matrices[1][i][j] = "D"
                elif(score == matrices[0][i-1][j] + scoringMatrix[0][len(scoringMatrix[0])-1]):
                    matrices[1][i][j] = "U"
                elif(score == matrices[0][i][j-1] + scoringMatrix[len(scoringMatrix[0])-1][0]):
                    matrices[1][i][j] = "L"
                
                matrices[0][i][j] = score
                
                if(score > bestScore[0]):
                    bestScore[0] = score
                    if(Ashort): 
                        bestScore[1] = (i-1,j-1+count)
                    else:
                        bestScore[1] = (i-1+count,j-1)
                    
                if(Ashort):                                 #"forget" the previous column
                    matrices[0][i-1][j-1] = matrices[0][i-1][j] 
                    matrices[1][i-1][j-1] = matrices[1][i-1][j] 
                else:                                       #"forget" the previous row
                    matrices[0][i-1][j-1] = matrices[0][i][j-1] 
                    matrices[1][i-1][j-1] = matrices[1][i][j-1]
        if(Ashort):                                 #"forget" the last remaining colunm cell
            matrices[0][len(matrices[0])-1][len(matrices[0][0])-2] = matrices[0][len(matrices[0])-1][len(matrices[0][0])-1] 
            matrices[1][len(matrices[0])-1][len(matrices[0][0])-2] = matrices[1][len(matrices[0])-1][len(matrices[0][0])-1]
        else:                                       #"forget" the last remaining row cell 
            matrices[0][len(matrices[0])-2][len(matrices[0][0])-1] = matrices[0][len(matrices[0])-1][len(matrices[0][0])-1] 
            matrices[1][len(matrices[0])-2][len(matrices[0][0])-1] = matrices[1][len(matrices[0])-1][len(matrices[0][0])-1]
        #print(np.matrix(matrices[1]))
        count += 1
    #print(np.matrix(matrices[0]))
    return bestScore, matrices

def findMidpoint(alphabet, scoringMatrix, seqA, seqB, matrices, lenA, lenB):
    if(lenB == 1):
        return getMatch(alphabet, scoringMatrix, seqA, seqB[0], i-1, j-1+count)
    #else:
        
def findIndicesLin(alphabet, scoringMatrix, seqA, seqB):
    lenA = len(seqA)
    lenB = len(seqB)
    
    global matchIndexA
    matchIndexA = []
    global matchIndexB
    matchIndexB = []
    
    #print("lenA ", lenA)
    #print("lenB ", lenB)
    if(lenA == 1 or lenB == 1):
        alignment = align(alphabet, scoringMatrix, seqA, seqB)
        backtr = alignment[1][1]
        start = alignment[0]
        
        m = start[1][0] + 1
        n = start[1][1] + 1
        #print(np.matrix(backtr))
        while(backtr[m][n] != "FN"):
            if backtr[m][n] == "L":
                n -= 1
            elif backtr[m][n] == "U":
                m -= 1
            elif backtr[m][n] == "D":
                matchIndexA.insert(0,m-1)
                matchIndexB.insert(0,n-1)
                print(matchIndexA, matchIndexB)
                m -= 1
                n -= 1
            elif(backtr[m][n] == "S"): #end of local alignment
                break
    else:
        findIndicesLin(alphabet, scoringMatrix, seqA[0:lenA//2], seqB[0:lenB//2])
        findIndicesLin(alphabet, scoringMatrix, seqA[lenA//2:lenA], seqB[lenB//2:lenB])
    print(matchIndexA, matchIndexB)

def dynproglin(alphabet, scoringMatrix, seqA, seqB):
    if(checkInput(alphabet, scoringMatrix, seqA, seqB) == -1):
        return -1
    
    bestScore = align(alphabet, scoringMatrix, seqA, seqB)[0]
    
    seqA = seqA[0:bestScore[1][0]]
    seqB = seqB[0:bestScore[1][1]]
    indices = findIndicesLin(alphabet, scoringMatrix, seqA, seqB)
    
    #print(np.matrix(matrices[0]))
    #print(np.matrix(matrices[1]))
    print("Best Local Score is:  ", bestScore[0])
    #print("Resulting indices:  ", indices[0],indices[1])

#dynproglin("ACTG", [[1,-1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-1,-2],[-1,-1,-1,1,-2],[-2,-2,-2,-2,-2]], "AAAC", "AGC")
dynproglin("ACTG", [[1,-1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-1,-2],[-1,-1,-1,1,-2],[-2,-2,-2,-2,-2]], "TAATA", "TACTAA")
#dynproglin("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAACA", "CBACCCBA")
         
"""
lenA = len(seqA)
lenB = len(seqB)

matrices = initializeMatricesLin(lenA, lenB)
bestScore = [0, (0,0)] #stores best local score and the endpoint symbols
if(lenA <= lenB):
    Ashort = True
else:
    Ashort = False

short = min(lenA, lenB)
count = 0
while(count < short):
    for i in range(1, len(matrices[0])):
        for j in range(1, len(matrices[0][0])):
            #diagonal, up, left and fresh start respectively
            print(seqA[i-1], seqB[j-1+count])
            if(Ashort):
                matchScore = getMatch(alphabet, scoringMatrix, seqA, seqB, i-1, j-1+count)
            else:
                matchScore = getMatch(alphabet, scoringMatrix, seqA, seqB, i-1+count, j-1)
                
            score = max(matchScore + matrices[0][i-1][j-1], matrices[0][i-1][j] + scoringMatrix[0][len(scoringMatrix[0])-1], matrices[0][i][j-1] + scoringMatrix[len(scoringMatrix[0])-1][0], 0)

            if(score == 0):
                matrices[1][i][j] = "S"
            elif(score == matchScore + matrices[0][i-1][j-1]):
                matrices[1][i][j] = "D"
            elif(score == matrices[0][i-1][j] + scoringMatrix[0][len(scoringMatrix[0])-1]):
                matrices[1][i][j] = "U"
            elif(score == matrices[0][i][j-1] + scoringMatrix[len(scoringMatrix[0])-1][0]):
                matrices[1][i][j] = "L"
            
            matrices[0][i][j] = score
            
            if(score > bestScore[0]):
                bestScore[0] = score
                if(Ashort): 
                    bestScore[1] = (i-1,j-1+count)
                else:
                    bestScore[1] = (i-1+count,j-1)
                
            if(Ashort):                                 #"forget" the previous column
                matrices[0][i-1][j-1] = matrices[0][i-1][j] 
                matrices[1][i-1][j-1] = matrices[1][i-1][j] 
            else:                                       #"forget" the previous row
                matrices[0][i-1][j-1] = matrices[0][i][j-1] 
                matrices[1][i-1][j-1] = matrices[1][i][j-1]
    if(Ashort):                                 #"forget" the last remaining colunm cell
        matrices[0][len(matrices[0])-1][len(matrices[0][0])-2] = matrices[0][len(matrices[0])-1][len(matrices[0][0])-1] 
        matrices[1][len(matrices[0])-1][len(matrices[0][0])-2] = matrices[1][len(matrices[0])-1][len(matrices[0][0])-1]
    else:                                       #"forget" the last remaining row cell 
        matrices[0][len(matrices[0])-2][len(matrices[0][0])-1] = matrices[0][len(matrices[0])-1][len(matrices[0][0])-1] 
        matrices[1][len(matrices[0])-2][len(matrices[0][0])-1] = matrices[1][len(matrices[0])-1][len(matrices[0][0])-1]
    #print(np.matrix(matrices[1]))
    count += 1

#print(np.matrix(matrices[0]))
#print(np.matrix(matrices[1]))
print("Best Local Score is:  ", bestScore[0])
"""
            
            
            
            
            
            
            
            
            