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
    #print(np.matrix(matrices[0]))
    #print(np.matrix(matrices[1]))
    print("Best Local Score is:  ", bestScore[0])
    print("Resulting indices:  ", indices[0],indices[1])
    return bestScore[0],indices[0],indices[1]

dynprog("ACTG", [[1,-1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-1,-2],[-1,-1,-1,1,-2],[-2,-2,-2,-2,-2]], "AAAC", "AGC")
dynprog("ACTG", [[1,-1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-1,-2],[-1,-1,-1,1,-2],[-2,-2,-2,-2,-2]], "TAATA", "TACTAA")
dynprog("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAACA", "CBACCCBA")

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            