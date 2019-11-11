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
            
            m = 0
            for symbol in alphabet:
                if(seqA[i-2] == symbol):
                    break
                m += 1
            
            n = 0
            for symbol in alphabet:
                if(seqB[j-2] == symbol):
                    break
                n += 1
            
            score = max(matchScore + matrices[0][i-1][j-1], matrices[0][i-1][j] + scoringMatrix[m][len(scoringMatrix[m])-1], matrices[0][i][j-1] + scoringMatrix[len(scoringMatrix[n])-1][n], 0)

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

dynprog("ACTG", [[1,-1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-1,-2],[-1,-1,-1,1,-2],[-2,-2,-2,-2,-2]], "AAAC", "AGC")
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
#dynproglin("ACTG", [[1,-1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-1,-2],[-1,-1,-1,1,-2],[-2,-2,-2,-2,-2]], "TAATA", "TACTAA")
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

def findDiagonals(seqA, seqB, ktup):
    matchPairs = []
    
    lenA = len(seqA)
    lenB = len(seqB)
    
    for i in range(lenA-ktup+1):
        for j in range(lenB-ktup+1):
            if(seqA[i:i+ktup] == seqB[j:j+ktup]):
                matchPairs.insert(0,(i,j))
    
    #print(matchPairs)   

    diagonals = {}
    
    for pair in matchPairs:
        diff = pair[0] - pair[1]
        if(diff not in diagonals):
            diagonals[diff] = [(pair[0],pair[1])]
        else:
            #print((pair[0],pair[1]))
            diagonals[diff].append((pair[0],pair[1]))
    
    #print(diagonals)
     
    return diagonals

def extendDiagonals(alphabet, scoringMatrix, seqA, seqB, diagonals):
    for diagonal in diagonals:
        extendedDiagonal = diagonals[diagonal].copy()
        for pair in diagonals[diagonal]:
            canForward = True
            canBackward = True
            currentForward = pair
            currentBackward = pair
            #print(pair)
            while(canForward or canBackward):
                if(currentForward[0]+1 >= len(seqA) or currentForward[1]+1 >= len(seqB)):
                    #print(currentForward[0]+1,currentForward[1]+1)
                    canForward = False
                if(currentBackward[0]-1 <= -1 or currentBackward[1]-1 <= -1):
                    canBackward = False
                
                if(canForward):
                    #print("forward",canForward)
                    #print(currentForward[0]+1,currentForward[1]+1)
                    if(getMatch(alphabet, scoringMatrix, seqA, seqB, currentForward[0]+1,currentForward[1]+1) < 0):
                        canForward = False
                    else:
                        currentForward = (currentForward[0]+1,currentForward[1]+1)
                        if(currentForward not in extendedDiagonal):
                            extendedDiagonal.append((currentForward[0],currentForward[1]))
                        
                if(canBackward):
                    #print("backward",canBackward)
                    if(getMatch(alphabet, scoringMatrix, seqA, seqB, currentBackward[0]-1,currentBackward[1]-1) < 0):
                        canBackward = False
                    else:
                        currentBackward = (currentBackward[0]-1,currentBackward[1]-1)
                        if(currentBackward not in extendedDiagonal):
                            extendedDiagonal.insert(0,(currentBackward[0],currentBackward[1]))
        
        extendedDiagonal.sort()      
        diagonals[diagonal] = extendedDiagonal
                
    print(diagonals)
    return diagonals

def heuralign(alphabet, scoringMatrix, seqA, seqB):
    if(checkInput(alphabet, scoringMatrix, seqA, seqB) == -1):
        return -1
    
    ktup = 2
    
    diagonals = findDiagonals(seqA, seqB, ktup)
    #print(diagonals)
    diagonals = extendDiagonals(alphabet, scoringMatrix, seqA, seqB, diagonals)
    
    width = 7
    
    bestAlignment = (None, [], [])
    
    for diagonal in diagonals:
        
        usedWidthUp = 0
        usedWidthLeft = 0
        
        pairs = diagonals[diagonal]
        
        indicesA = []
        indicesB = []
        
        currentScore = 0
        searchScore = 0
        
        for pair in pairs: 
            currentScore += getMatch(alphabet, scoringMatrix, seqA, seqB, pair[0], pair[1])
            indicesA.append(pair[0])
            indicesB.append(pair[1])
            
            nonMatchExtend = (pair[0]-1,pair[1]-1)
            while(nonMatchExtend not in pairs and nonMatchExtend[0] >= 0 and nonMatchExtend[1] >= 0):
                if((nonMatchExtend[0]+1,nonMatchExtend[1]+1) == pairs[0]):
                    break
                currentScore += getMatch(alphabet, scoringMatrix, seqA, seqB, nonMatchExtend[0], nonMatchExtend[1])
                #print("n",nonMatchExtend)
                nonMatchExtend = (nonMatchExtend[0]-1,nonMatchExtend[1]-1)
            
        
        searchPos = (pairs[0][0]-1, pairs[0][1]-1)
        print(searchPos, currentScore)
        while((usedWidthUp < width//2 and searchPos[0] >= 0) and (usedWidthLeft < width//2 and searchPos[1] >= 0)):
            #print(searchPos, searchPos[0] >= 0)
            stepScoreDiagonal = getMatch(alphabet, scoringMatrix, seqA, seqB, searchPos[0], searchPos[1])
            #print("sc",stepScoreDiagonal,seqA[searchPos[0]],seqB[searchPos[1]])
            if(currentScore + searchScore + stepScoreDiagonal >= currentScore):
                print("f",searchScore,stepScoreDiagonal)
                currentScore += stepScoreDiagonal
                currentScore += searchScore
                indicesA.insert(0,searchPos[0])
                indicesB.insert(0,searchPos[1])
                searchPos = (searchPos[0]-1, searchPos[1]-1)
                searchScore = 0
            else:
                stepScoreUp = -90000
                stepScoreLeft = -90000
                #find the correct symbol to match gap against
                if(usedWidthUp < width//2 and searchPos[0] >= 0):
                    i = 0
                    for symbol in alphabet:
                        if(seqA[searchPos[0]] == symbol):
                            break
                        i += 1
                    
                    stepScoreUp = scoringMatrix[i][len(scoringMatrix[i])-1]
                    print(stepScoreUp, searchScore)
                
                if(usedWidthLeft < width//2 and searchPos[1] >= 0):
                    j = 0
                    for symbol in alphabet:
                        if(seqB[searchPos[1]] == symbol):
                            break
                        j += 1
                    
                    stepScoreLeft = scoringMatrix[j][len(scoringMatrix[j])-1]
                
                if(stepScoreUp >= stepScoreLeft):
                    searchScore += stepScoreUp
                    searchPos = (searchPos[0]-1, searchPos[1])
                    usedWidthUp += 1
                    usedWidthLeft -= 1
                else:
                    searchScore += stepScoreLeft
                    searchPos = (searchPos[0], searchPos[1]-1)
                    usedWidthUp -= 1
                    usedWidthLeft += 1
        #print("p",pair[0],pair[1])
        stepScoreDiagonal = getMatch(alphabet, scoringMatrix, seqA, seqB, searchPos[0], searchPos[1])
        if(currentScore + searchScore + stepScoreDiagonal >= currentScore and searchPos[1] >= 0 and searchPos[0] >= 0):
            print(searchPos, currentScore, stepScoreDiagonal)
            currentScore += stepScoreDiagonal
            currentScore += searchScore
            indicesA.insert(0,searchPos[0])
            indicesB.insert(0,searchPos[1])
            searchPos = (searchPos[0]-1, searchPos[1]-1)
            searchScore = 0
        
        if(bestAlignment[0] == None):
            bestAlignment = (currentScore, indicesA, indicesB)
        elif(bestAlignment[0] <= currentScore):
            bestAlignment = (currentScore, indicesA, indicesB)
        else:
            pass
    
    print("Best Local Score found is:  ", bestAlignment[0])
    print("Resulting indices:  ", bestAlignment[1], bestAlignment[2])
    return bestAlignment[0], bestAlignment[1], bestAlignment[2]


"""
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
"""






#heuralign("ACTG", [[1,-1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-1,-2],[-1,-1,-1,1,-2],[-2,-2,-2,-2,-2]], "AAAC", "AGC")
#heuralign("ACTG", [[1,-1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-1,-2],[-1,-1,-1,1,-2],[-2,-2,-2,-2,-2]], "TAATA", "TACTAA")
#heuralign("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAACA", "CBACCCBA")

"TAATA"
"TACTAA  "

"AABBAACA"
"  CB-ACCCBA"

"[0,1,2][3,4,5]"


            
            
            
            
            