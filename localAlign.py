# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 13:38:50 2019

@author: Erdal Guclu
"""
import numpy as np

#initialise the score and backtrack matrices for local alignments
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
    
    backtr[1][1] = "S"
        
    return [matrix, backtr]

#makes sure the input is valid
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

#matches and finds the resulting score between two characters in a scoring matrix
def getMatch(alphabet, scoringMatrix, A, B, i, j):
    symbolA = ""
    symbolB = ""
    for k in range(0, len(alphabet)):
        if(A[i] == alphabet[k]):
            symbolA = k
        if(B[j] == alphabet[k]):
            symbolB = k
            
    return scoringMatrix[symbolA][symbolB]
    
#generates the indices for the best local alignment starting at (i,j)
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

#finds the optimal local alignment and returns its score
def dynprog(alphabet, scoringMatrix, seqA, seqB):
    if(checkInput(alphabet, scoringMatrix, seqA, seqB) == -1):
        return -1
    lenA = len(seqA)
    lenB = len(seqB)
    if(lenA == 0 or lenB == 0):
        print("Error: Cannot align empty string")
        return 0,[],[]
    matrices = initializeMatrices(lenA, lenB)
    bestScore = [0, (0,0)] #stores best local score and the endpoint symbols
    
    for i in range(2, lenA+2):
        for j in range(2, lenB+2):
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
            
            #diagonal, up, left and fresh start respectively
            score = max(matchScore + matrices[0][i-1][j-1], matrices[0][i-1][j] + scoringMatrix[m][len(scoringMatrix[m])-1], matrices[0][i][j-1] + scoringMatrix[n][len(scoringMatrix[n])-1], 0)
            
            #update backtrack with relevant entry
            if(score == 0):
                matrices[1][i][j] = "S"
            elif(score == matchScore + matrices[0][i-1][j-1]):
                matrices[1][i][j] = "D"
            elif(score == matrices[0][i-1][j] + scoringMatrix[m][len(scoringMatrix[m])-1]):
                matrices[1][i][j] = "U"
            elif(score == matrices[0][i][j-1] + scoringMatrix[n][len(scoringMatrix[n])-1]):
                matrices[1][i][j] = "L"
            
            matrices[0][i][j] = score
            if(score > bestScore[0]):
                bestScore[0] = score
                bestScore[1] = (i-2,j-2)
    
    indices = findIndices(matrices[1], lenA, lenB, bestScore[1][0], bestScore[1][1])
    print(bestScore[0],indices[0],indices[1])
    return bestScore[0],indices[0],indices[1]

#dynprog("ACTG", [[1,-1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-1,-2],[-1,-1,-1,1,-2],[-2,-2,-2,-2,-2]], "AAAC", "AGC")
#dynprog("ACTG", [[1,-1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-1,-2],[-1,-1,-1,1,-2],[-2,-2,-2,-2,-2]], "TAATA", "TACTAA")
#dynprog("ACTG", [[2,-1,-1,-1,-2],[-1,2,-1,-2],[-1,-1,2,-1,-2],[-1,-1,-1,2,-2],[-2,-2,-2,-2,0]], "AGTACGCA", "TATGC")
#dynprog("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAACA", "CBACCCBA")     
#dynprog("ABCD", [[ 1,-5,-5,-5,-1],[-5, 1,-5,-5,-1],[-5,-5, 5,-5,-4],[-5,-5,-5, 6,-4],[-1,-1,-4,-4,-9]],"AAAAACCDDCCDDAAAAACC","CCAAADDAAAACCAAADDCCAAAA")
#dynprog("ABCD", [[ 1,-5,-5,-5,-1],[-5, 1,-5,-5,-1],[-5,-5, 5,-5,-4],[-5,-5,-5, 6,-4],[-1,-1,-4,-4,-9]],"AACAAADAAAACAADAADAAA","CDCDDD")
#dynprog("ABCD", [[ 1,-5,-5,-5,-1],[-5, 1,-5,-5,-1],[-5,-5, 5,-5,-4],[-5,-5,-5, 6,-4],[-1,-1,-4,-4,-9]],"DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD","DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD")

#generates the indices for the global alignment
def findIndicesGlobal(backtr, lenA, lenB):
    m = len(backtr) - 1
    n = len(backtr[0]) - 1
    matchIndexA = []
    matchIndexB = []
    
    while(backtr[m][n] != "FN"):
        if(backtr[m][n] == "L"):
            n -= 1
        elif(backtr[m][n] == "U"):
            m -= 1
        elif(backtr[m][n] == "D"):
            matchIndexA.insert(0,m-2)
            matchIndexB.insert(0,n-2)
            m -= 1
            n -= 1
        elif(backtr[m][n] == "FN"): #end of alignment
            break
        
    return matchIndexA, matchIndexB

#initalises the score and backtrack matrix for global alignments
def initializeMatricesGlobal(alphabet, scoringMatrix, seqA, seqB):
    lenA = len(seqA)
    lenB = len(seqB)
    
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
        
        m = 0
        for symbol in alphabet:
            if(seqA[i-2] == symbol):
                break
            m += 1
        matrix[i][1] = matrix[i-1][1] + scoringMatrix[m][len(scoringMatrix[m])-1]
        
        backtr[i][0] = i-2
        backtr[i][1] = "U"
    for j in range(2, lenB+2): #initialise indices of seqB
        matrix[0][j] = j-2
        
        n = 0
        for symbol in alphabet:
            if(seqB[j-2] == symbol):
                break
            n += 1
        matrix[1][j] = matrix[1][j-1] + scoringMatrix[n][len(scoringMatrix[n])-1]
        
        backtr[0][j] = j-2
        backtr[1][j] = "L"
    
    backtr[1][1] = "FN" #end marker for the backtrack
        
    return [matrix, backtr]

#finds the best global alignment and returns both the alignment and its score
def dynprogGlobal(alphabet, scoringMatrix, seqA, seqB):
    if(checkInput(alphabet, scoringMatrix, seqA, seqB) == -1):
        return -1
    
    lenA = len(seqA)
    lenB = len(seqB)
    if(lenA == 0 or lenB == 0):
        print("Error: Cannot align empty string")
        return 0,[],[]
    
    matrices = initializeMatricesGlobal(alphabet, scoringMatrix, seqA, seqB)
    
    for i in range(2, lenA+2):
        for j in range(2, lenB+2):
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
            
            #diagonal, up and left respectively
            score = max(matchScore + matrices[0][i-1][j-1], matrices[0][i-1][j] + scoringMatrix[m][len(scoringMatrix[m])-1], matrices[0][i][j-1] + scoringMatrix[n][len(scoringMatrix[n])-1])
            
            #add relevant entry to backtrack matrix
            if(score == matchScore + matrices[0][i-1][j-1]):
                matrices[1][i][j] = "D"
            elif(score == matrices[0][i-1][j] + scoringMatrix[m][len(scoringMatrix[m])-1]):
                matrices[1][i][j] = "U"
            elif(score == matrices[0][i][j-1] + scoringMatrix[n][len(scoringMatrix[n])-1]):
                matrices[1][i][j] = "L"
            
            matrices[0][i][j] = score
    
    indices = findIndicesGlobal(matrices[1], lenA, lenB)
    return matrices[0][len(matrices[0]) - 1][len(matrices[0][0]) - 1],indices[0],indices[1]

#finds the best local alignment score in linear space
def getScoreLocal(alphabet, scoringMatrix, seqA, seqB):
    #initialise matrix
    matrix = [[None for i in range(len(seqB)+1)] for j in range(2)]
    matrix[0][0] = 0
    for i in range(1, len(seqB)+1): 
        n = 0
        for symbol in alphabet:
            if(seqB[i-1] == symbol):
                break
            n += 1
            
        matrix[0][i] = matrix[0][i-1] + scoringMatrix[n][len(scoringMatrix[n])-1]
    
    highestScore = [0, (0,0)] #best score and location of best score
    
    for i in range(0, len(seqA)):
        for j in range(0, len(seqB)+1):
            n = 0
            if(j != 0):
                #wont match anything on leftmost column
                matchScore = getMatch(alphabet, scoringMatrix, seqA, seqB, i, j-1)
                for symbol in alphabet:
                    if(seqB[j-1] == symbol):
                        break
                    n += 1
            m = 0
            for symbol in alphabet:
                if(seqA[i] == symbol):
                    break
                m += 1
            
            if(j == 0):
                #leftmost column
                matrix[1][j] = 0
            else:
                #diagonal, up, left and fresh start 
                matrix[1][j] = max(matchScore + matrix[0][j-1], matrix[0][j] + scoringMatrix[m][len(scoringMatrix[m])-1], matrix[1][j-1] + scoringMatrix[n][len(scoringMatrix[n])-1], 0)
            
            if(matrix[1][j] >= highestScore[0]):
                highestScore[0] = matrix[1][j]
                highestScore[1] = (i+1,j) #return position in terms of matrix indices
                
        matrix[0][:] = matrix[1][:]
    #print(highestScore)
    return matrix[1], highestScore

#computes the global score matrix in linear space and returns last line as well as position of best score in the matrix
def getLineGlobal(alphabet, scoringMatrix, seqA, seqB):
    #initialise matrix
    matrix = [[None for i in range(len(seqB)+1)] for j in range(2)]
    matrix[0][0] = 0
    
    for i in range(1, len(seqB)+1): 
        n = 0
        for symbol in alphabet:
            if(seqB[i-1] == symbol):
                break
            n += 1
            
        matrix[0][i] = matrix[0][i-1] + scoringMatrix[n][len(scoringMatrix[n])-1]
    
    highestScore = [0, (0,0)] #the best score and associated position
    
    for i in range(0, len(seqA)):
        for j in range(0, len(seqB)+1):
            n = 0
            if(j != 0):
                #wont be matching anything when on leftmost column
                matchScore = getMatch(alphabet, scoringMatrix, seqA, seqB, i, j-1)
                for symbol in alphabet:
                    if(seqB[j-1] == symbol):
                        break
                    n += 1
            m = 0
            for symbol in alphabet:
                if(seqA[i] == symbol):
                    break
                m += 1
                
            
            if(j == 0):
                #leftmost column
                matrix[1][j] = matrix[0][j] + scoringMatrix[m][len(scoringMatrix[m])-1]
            else:
                #diagonal, up, left 
                matrix[1][j] = max(matchScore + matrix[0][j-1], matrix[0][j] + scoringMatrix[m][len(scoringMatrix[m])-1], matrix[1][j-1] + scoringMatrix[n][len(scoringMatrix[n])-1])

            if(matrix[1][j] >= highestScore[0]):
                highestScore[0] = matrix[1][j]
                highestScore[1] = (i,j-1)
        
        matrix[0][:] = matrix[1][:]
    return matrix[1], highestScore

#recursive linear space algorithm for local alignments with wrapper function
def dynproglin(alphabet, scoringMatrix, seqA, seqB):
    if(checkInput(alphabet, scoringMatrix, seqA, seqB) == -1):
        return -1
    
    #find the endpoint of the local alignment
    localEnd = getLineGlobal(alphabet, scoringMatrix, seqA, seqB)
    seqA = seqA[0:localEnd[1][1][0]+1]
    seqB = seqB[0:localEnd[1][1][1]+1]
    
    #find the score of best local alignment
    topScore = getScoreLocal(alphabet, scoringMatrix, seqA, seqB)[1][0]
    
    seqA = seqA[::-1]
    seqB = seqB[::-1]
    
    lenA = len(seqA)
    lenB = len(seqB)
    
    #find the startpoint of the local alignment using the reversed sequences
    localStart = getLineGlobal(alphabet, scoringMatrix, seqA, seqB)
    seqA = seqA[0:localStart[1][1][0]+1]
    seqB = seqB[0:localStart[1][1][1]+1]
    
    startIndex = localStart[1][1]
    
    seqA = seqA[::-1]
    seqB = seqB[::-1]
    
    #tracks the position in the sequences through the recursion tree
    offsetA = lenA - startIndex[0] - 1
    offsetB = lenB - startIndex[1] - 1
    
    #recursive private function
    def dynproglin(alphabet, scoringMatrix, seqA, seqB, offsetA, offsetB):
        
        if(len(seqA) == 0 or len(seqB) == 0):
            return [[], []]
        elif(len(seqA) == 1 or len(seqB) == 1):
            #use global algorithm on the base case
            bestScore = dynprogGlobal(alphabet, scoringMatrix, seqA, seqB)
            indices = [[x + offsetA for x in bestScore[1]],[y + offsetB for y in bestScore[2]]]
            
            return indices
        else:
            lenA = len(seqA)
            lenB = len(seqB)
        
            Amid = lenA // 2
            
            #find the last lines of score matrix before and after midpoint
            L = getLineGlobal(alphabet, scoringMatrix, seqA[0:Amid], seqB)
            ScoreL = L[0]
            partA = seqA[Amid:lenA]
            R = getLineGlobal(alphabet, scoringMatrix, partA[::-1], seqB[::-1])
            ScoreR = R[0]
            
            #compute Bmid to get location of midpoint
            argMax = [sum(x) for x in zip(ScoreL, ScoreR[::-1])]
            Bmid = argMax.index(max(argMax))
            
            #recurse with respect to the found midpoint
            bestScoreL = dynproglin(alphabet, scoringMatrix, seqA[0:Amid], seqB[0:Bmid], offsetA, offsetB)
            bestScoreR = dynproglin(alphabet, scoringMatrix, seqA[Amid:lenA], seqB[Bmid:lenB], offsetA+Amid, offsetB+Bmid)
            
            indices = [[], []]
            indices[0] = bestScoreL[0] + bestScoreR[0]
            indices[1] = bestScoreL[1] + bestScoreR[1]
            
            return indices
        
    indices = dynproglin(alphabet, scoringMatrix, seqA, seqB, offsetA, offsetB)
    print(topScore, indices[0], indices[1])
    return topScore, indices[0], indices[1]

#dynproglin("ACTG", [[1,-1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-1,-2],[-1,-1,-1,1,-2],[-2,-2,-2,-2,-2]], "AAAC", "AGC")
#dynproglin("ACTG", [[1,-1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-1,-2],[-1,-1,-1,1,-2],[-2,-2,-2,-2,-2]], "TAATA", "TACTAA")
#dynproglin("ACTG", [[2,-1,-1,-1,-2],[-1,2,-1,-2],[-1,-1,2,-1,-2],[-1,-1,-1,2,-2],[-2,-2,-2,-2,0]], "AGTACGCA", "TATGC")
#dynproglin("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAACA", "CBACCCBA")     
#dynproglin("ABCD", [[ 1,-5,-5,-5,-1],[-5, 1,-5,-5,-1],[-5,-5, 5,-5,-4],[-5,-5,-5, 6,-4],[-1,-1,-4,-4,-9]],"AAAAACCDDCCDDAAAAACC","CCAAADDAAAACCAAADDCCAAAA")
#dynproglin("ABCD", [[ 1,-5,-5,-5,-1],[-5, 1,-5,-5,-1],[-5,-5, 5,-5,-4],[-5,-5,-5, 6,-4],[-1,-1,-4,-4,-9]],"AACAAADAAAACAADAADAAA","CDCDDD")
#dynproglin("ABCD", [[ 1,-5,-5,-5,-1],[-5, 1,-5,-5,-1],[-5,-5, 5,-5,-4],[-5,-5,-5, 6,-4],[-1,-1,-4,-4,-9]],"DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD","DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD")

#uses the two sequences and a ktup value to find all the diagonals 
def findDiagonals(seqA, seqB, ktup):
    matchPairs = []
    lenA = len(seqA)
    lenB = len(seqB)
    
    #find all the matching subsequences of length ktup
    while(matchPairs == []):
        for i in range(lenA-ktup+1):
            for j in range(lenB-ktup+1):
                if(seqA[i:i+ktup] == seqB[j:j+ktup]):
                    matchPairs.insert(0,(i,j))
        ktup -= 1
        if(ktup == 0):
            break
        
    #place all the match pairs into their corresponding diagonals
    diagonals = {}
    for pair in matchPairs:
        diff = pair[0] - pair[1]
        if(diff not in diagonals):
            diagonals[diff] = [(pair[0],pair[1])]
        else:
            diagonals[diff].append((pair[0],pair[1]))
     
    return diagonals

#extends all found diagonals so that they are connected
def extendDiagonals(alphabet, scoringMatrix, seqA, seqB, diagonals):
    for diagonal in diagonals:
        extendedDiagonal = diagonals[diagonal].copy()
        for pair in diagonals[diagonal]:
            canForward = True  #extend the diagonal diagonally upwards on the matrix
            canBackward = True #extend the diagonal diagonally downwards on the matrix
            currentForward = pair
            currentBackward = pair
            
            while(canForward or canBackward):
                #stop extending in appropriate direction if one or both sequences are exhausted
                if(currentForward[0]+1 >= len(seqA) or currentForward[1]+1 >= len(seqB)):
                    canForward = False
                if(currentBackward[0]-1 <= -1 or currentBackward[1]-1 <= -1):
                    canBackward = False
                if(canForward):
                    if(getMatch(alphabet, scoringMatrix, seqA, seqB, currentForward[0]+1,currentForward[1]+1) < 0):
                        canForward = False
                        for pair in diagonals[diagonal]: #if not a positive match it could still be a mismatch in the alignment in the middle of the diagonal
                            
                            if((currentForward[0]+1,currentForward[1]+1) < pair): 
                                canForward = True
                                currentForward = (currentForward[0]+1,currentForward[1]+1)
                                if(currentForward not in extendedDiagonal):
                                    extendedDiagonal.append((currentForward[0],currentForward[1]))
                    else: #is a matching symbol
                        currentForward = (currentForward[0]+1,currentForward[1]+1)
                        if(currentForward not in extendedDiagonal):
                            extendedDiagonal.insert(0,(currentForward[0],currentForward[1]))
                if(canBackward):
                    if(getMatch(alphabet, scoringMatrix, seqA, seqB, currentBackward[0]-1,currentBackward[1]-1) < 0): #not a match
                        canBackward = False
                        for pair in diagonals[diagonal]: #if not a positive match it could still be a mismatch in the alignment in the middle of the diagonal
                            
                            if((currentBackward[0]-1,currentBackward[1]-1) > pair):
                                canBackward = True
                                currentBackward = (currentBackward[0]-1,currentBackward[1]-1)
                                if(currentBackward not in extendedDiagonal):
                                    extendedDiagonal.append((currentBackward[0],currentBackward[1]))
                    else: #is a matching symbol
                        currentBackward = (currentBackward[0]-1,currentBackward[1]-1)
                        if(currentBackward not in extendedDiagonal):
                            extendedDiagonal.append((currentBackward[0],currentBackward[1]))
        
        #sort and update the new diagonal
        extendedDiagonal.sort()      
        diagonals[diagonal] = extendedDiagonal
    return diagonals

def heuralign(alphabet, scoringMatrix, seqA, seqB):
    if(checkInput(alphabet, scoringMatrix, seqA, seqB) == -1):
        return -1
    
    ktup = 6 #length of substrings to search matches for
    
    diagonals = findDiagonals(seqA, seqB, ktup)
    diagonals = extendDiagonals(alphabet, scoringMatrix, seqA, seqB, diagonals)
    
    width = 7 #width-1//2 is the number of rows upward and leftwards banded DP can search
    
    bestAlignment = (None, [], [])
    
    for diagonal in diagonals:
        #keeps track of where the search is happening
        usedWidthUp = 0
        usedWidthLeft = 0
        
        pairs = diagonals[diagonal]
        
        indicesA = []
        indicesB = []
        
        currentScore = 0 #the score that the current best alignment has
        searchScore = 0 #the score accumulated by the search path still not on the best alignment
        
        #find the score along the entire diagonal before starting the search
        for pair in pairs:
            #the members of the diagonal will be on the current alignment
            currentScore += getMatch(alphabet, scoringMatrix, seqA, seqB, pair[0], pair[1])
            indicesA.append(pair[0])
            indicesB.append(pair[1])
            
        searchPos = (pairs[0][0]-1, pairs[0][1]-1) #current location of search
        
        #loop as long as we are within bounds of the bandedDP
        while((usedWidthUp < width//2 and searchPos[0] >= 0) and (usedWidthLeft < width//2 and searchPos[1] >= 0)):
            stepScoreDiagonal = getMatch(alphabet, scoringMatrix, seqA, seqB, searchPos[0], searchPos[1]) #search in the diagonal
            if(currentScore + searchScore + stepScoreDiagonal >= currentScore): #diagonal results in a better alignment than the one we had
                currentScore += stepScoreDiagonal
                currentScore += searchScore
                indicesA.insert(0,searchPos[0])
                indicesB.insert(0,searchPos[1])
                searchPos = (searchPos[0]-1, searchPos[1]-1)
                searchScore = 0
            else: #if the diagonal is not better then we extend to other diagonals
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
                if(usedWidthLeft < width//2 and searchPos[1] >= 0):
                    j = 0
                    for symbol in alphabet:
                        if(seqB[searchPos[1]] == symbol):
                            break
                        j += 1
                        
                    stepScoreLeft = scoringMatrix[j][len(scoringMatrix[j])-1]
                
                #decide whether or not its better to search left or up
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
                    
        stepScoreDiagonal = getMatch(alphabet, scoringMatrix, seqA, seqB, searchPos[0], searchPos[1]) #attempt to find a match along new diagonal
        if(currentScore + searchScore + stepScoreDiagonal >= currentScore and searchPos[1] >= 0 and searchPos[0] >= 0): #diagonal results in a better alignment than the one we had
            currentScore += stepScoreDiagonal
            currentScore += searchScore
            indicesA.insert(0,searchPos[0])
            indicesB.insert(0,searchPos[1])
            searchPos = (searchPos[0]-1, searchPos[1]-1)
            searchScore = 0
        
        #only return the best alignment found accross all diagonals
        if(bestAlignment[0] == None):
            bestAlignment = (currentScore, indicesA, indicesB)
        elif(bestAlignment[0] <= currentScore):
            bestAlignment = (currentScore, indicesA, indicesB)
        else:
            pass
        
    print(bestAlignment[0], bestAlignment[1], bestAlignment[2])
    return bestAlignment[0], bestAlignment[1], bestAlignment[2]
    
#heuralign("ACTG", [[1,-1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-1,-2],[-1,-1,-1,1,-2],[-2,-2,-2,-2,-2]], "AAAC", "AGC")
#heuralign("ACTG", [[1,-1,-1,-1,-2],[-1,1,-1,-2],[-1,-1,1,-1,-2],[-1,-1,-1,1,-2],[-2,-2,-2,-2,-2]], "TAATA", "TACTAA")
#heuralign("ACTG", [[2,-1,-1,-1,-2],[-1,2,-1,-2],[-1,-1,2,-1,-2],[-1,-1,-1,2,-2],[-2,-2,-2,-2,0]], "AGTACGCA", "TATGC")
#heuralign("ABC", [[1,-1,-2,-1],[-1,2,-4,-1],[-2,-4,3,-2],[-1,-1,-2,0]], "AABBAACA", "CBACCCBA") 
#heuralign("ABCD", [[ 1,-5,-5,-5,-1],[-5, 1,-5,-5,-1],[-5,-5, 5,-5,-4],[-5,-5,-5, 6,-4],[-1,-1,-4,-4,-9]],"AAAAACCDDCCDDAAAAACC","CCAAADDAAAACCAAADDCCAAAA")
#heuralign("ABCD", [[ 1,-5,-5,-5,-1],[-5, 1,-5,-5,-1],[-5,-5, 5,-5,-4],[-5,-5,-5, 6,-4],[-1,-1,-4,-4,-9]],"AACAAADAAAACAADAADAAA","CDCDDD")
#heuralign("ABCD", [[ 1,-5,-5,-5,-1],[-5, 1,-5,-5,-1],[-5,-5, 5,-5,-4],[-5,-5,-5, 6,-4],[-1,-1,-4,-4,-9]],"DDCDDCCCDCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCDDDCDADCDCDCDCD","DDCDDCCCDCBCCCCDDDCDBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBDCDCDCDCD")

            
            
            
            
            
