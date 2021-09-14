from math import sqrt
from pprint import pprint

def integerRound(x):
    from math import floor, ceil
    a = floor(x)
    b = ceil(x)
    if x-a<=.1:
        result = a
    elif b-x<=.1:
        result = b
    else:
        print("Please apply the function to a smaller difference")
    return result

def dec2Round(x):
    y = 100*x
    remainder = y%1
    y = y - remainder
    if remainder > 0.1:
        modified_digit = 1
    else:
        modified_digit = 0
    y = y + modified_digit
    y = y/100
    return y

def displayMatrix(matrix): #display given matrix
    for row in matrix:
        print(row)

def transpose(A):
    n = len(A)
    A_T = [[0 for i in range(n)] for j in range(n)]
    for i in range(len(A)):
        for j in range(len(A[0])):
            A_T[j][i] = A[i][j]
    return A_T               
        
def applyGJ(A,B): #function to apply GJ elim
    for k in range(len(B)):
        #rows are pivotted
        if abs(A[k][k]) < 1.0e-6:
            #defining upper limit for element = 0
            for i in range(k+1, len(B)): 
                if abs(A[i][k]) > abs(A[k][k]):
                    #swap check
                    for j in range(k, len(B)):
                        #obtaining swapped rows
                        A[k][j], A[i][j] = A[i][j], A[k][j] 
                    B[k], B[i] = B[i], B[k] 
                    break
        A_kk = A[k][k]
        if A_kk == 0:        
            print("No distinct solution was found for this system of equations")
            return
        for j in range(k, len(B)): #column index of row is pivotted
            A[k][j] /= A_kk         #pivot row division for row ech form
        B[k] /= A_kk
        for i in range(len(B)):    #changed rows assigned new indices
            if i == k or A[i][k] == 0: continue
            factor = A[i][k]
            for j in range(k, len(B)): 
                #columns for subtraction assigned indices
                A[i][j] -= factor*A[k][j]
            B[i] -= factor * B[k]
    return B

def productMatrix(B,A): #finds product of two matrices
    try:               
        if  len(A) != len(B[0]): 
            print("Multiplication is undefined for given matrices!") 
        else:
            C = [[0 for i in range(len(A[0]))] for j in range(len(B))]
            for i in range(len(B)):       #rows of the matrix product.
                for j in range(len(A[0])):#columns of the matrix product.
                    for k in range(len(A)):
                        C[i][j] += dec2Round(B[i][k]*A[k][j])
            return C
    except TypeError:
        print("Invalid entries found in matrix!")

def findDeterminant(A):     #determinant function for use in applyGJ()
    if len(A) != len(A[0]): #square matrix check
        print("Determinant is undefined for square matrices")
    else:
        count = 0
        for i in range(len(A) - 1): 
            if abs(A[i][i]) < 1.0e-6:
                for j in range(i+1 , len(A)): 
                    if  abs(A[i][i]) < abs(A[j][i]): 
                        for k in range(i , len(A)):
                            #swapping the rows of the matrix
                            A[i][k], A[j][k] = A[j][k], A[i][k]
                            count += 1
            for j in range(i+1 , len(A)):
                if A[j][i] == 0: continue 
                result = A[j][i]/A[i][i] #dividing the rows.
                for k in range(i , len(A)):
                    A[j][k] -= result*A[i][k]
        initial = 1
        for j in range(len(A)):
            initial *= A[j][j] #product of diagonal elements of matrix
        initial *= (-1)**count
        print(initial) 

def inverseMatrix(A):
    #function for inverse matrix (3x3)
    M = [[ 0.00 for i in range(len(A))] for j in range(len(A))] 
    for i in range(3):
        for j in range(3):
            M[j][j] = 1.00
    for i in range(len(A)):
        A[i].extend(M[i])
    for k in range(len(A)):
        if abs(A[k][k]) < 1.0e-6:
            #GJ elim segment
            for i in range(k+1, len(A)):
                if abs(A[i][k]) > abs(A[k][k]):
                    for j in range(k,2*len(A)):
                        #swap rows
                        A[k][j], A[i][j] = A[i][j], A[k][j]
                    break
        count = A[k][k] #element is pivotted
        if count == 0:  #checking if pivot = 0
            print("The matrix does not have a defined inverse")
            return
        else:
            for j in range(k, 2*len(A)):
                #pivotted row columns
                A[k][j] /= count
            for i in range(len(A)):
                #substracted rows indiced
                if i == k or A[i][k] == 0: continue
                result = A[i][k] 
                for j in range(k, 2*len(A)): 
                    #columns for subtraction indiced
                    A[i][j] -= result*A[k][j]
                    
    A_inv = []
    
    for i in range(len(A)):
        blank_row = []
        for j in range(len(A),len(A[0])):
            blank_row.append(dec2Round(A[i][j]))
        A_inv.append(blank_row)
    return A_inv

#LU Decomposition function set
def swapRows(M, row_old, row_new, columns):
    #to swap rows, wherever needed
    tmp = []
    for i in range (0, int(columns)):
        tmp.append(M[int(row_old)][i])
        M[int(row_old)][i] = M[int(row_new)][i]
        M[int(row_new)][i] = tmp[i]
        
def unSwap(A,B,A_orig):
    #to undo swapping from inversion
    N = len(A)
    fix_list = []
    B_fixed = []
    
    for i in range(N):
        for j in range(N):
            if A[i]==A_orig[j]:
                fix_list.append(j)
    
    for elem in fix_list:
        B_fixed.append(B[elem])
        
    B = B_fixed
        
def crout(A):
    n = len(A)
    L = [[0.0 for i in range(n)] for j in range(n)]
    U = [[0.0 for i in range(n)] for j in range(n)]
    
    for i in range(n):
        U[i][i] = 1
        for j in range(i, n):
            tmp1 = A[j][i] 
            for k in range(i):
                tmp1 -= L[j][k]*U[k][i]
            L[j][i] = tmp1
        for j in range(i+1, n):
            tmp2 = A[i][j]
            for k in range(i):
                tmp2 -= L[i][k]*U[k][j]
            U[i][j] = tmp2/L[i][i]
    return (L, U)

def doolittle(A):
    n = len(A)
    L = [[0.0 for i in range(n)] for j in range(n)]
    U = [[0.0 for i in range(n)] for j in range(n)]
    
    for z in range(n):
        L[z][z] = 1
        f1 = 0
        f2 = 0
        f3 = 0
        for p in range(z):
            f1 = f1 + L[z][p]*U[p][z]
        U[z][z] = (A[z][z] - f1)
        for i in range(z+1, n):
            for p in range(z):
                f2 = f2 + L[z][p]*U[p][i]
            U[z][i] = (A[z][i] - f2)
        for k in range(z+1, n):
            for p in range(z):
                f3 = f3 + L[k][p]*U[p][z]
            L[k][z] = (A[k][z] - f3)/U[z][z]
    return (L, U)

def forwardSolve(L, b): #solves L.y = b
    y = [0 for i in range(len(b))]
    for i in range(len(b)):
        sumj = 0
        for j in range(i):
            sumj = sumj + L[i][j]*y[j]
        y[i] = (b[i] -sumj)/L[i][i]
    return y

def backwardSolve(U, y): #solves U.x = y
    n = len(y)
    x = [0 for i in range(len(y))]
    for i in range(n-1, -1, -1):
        sumj = 0
        for j in range(i+1, n):
            sumj = sumj + U[i][j] * x[j]
        x[i] = (y[i] - sumj)/U[i][i]
    return x 

def partPivot(M, m, r, c):
    #to get partially pivotted matrix for LU decomp
    global n, swaps
    n = 0
    swaps = 0
    pivot = M[int(m)][int(m)]
    for i in range (int(r)):         
        if pivot < M[int(i)][int(m)]:
            #same-column elements are checked
            pivot = M[int(i)][int(m)]
            n += 1
            swaps = i
    if swaps != 0: #swap if allowed
        swapRows(M, m, swaps, c)
    if int(pivot) == 0:
        print ("No unique solution")
        return None

def solveLU(A, b, method):
    if method==1:
        L, U = crout(A)
    elif method==2:
        L, U = doolittle(A)
    #print("L = " + str(L) + "\n")
    #print("U = " + str(U) + "\n")
    y = forwardSolve(L, b)
    x = backwardSolve(U, y)
    return x

def LU_inverse(M):
    #to find inverse of matrix using LU decomposition
    M_orig = [[M[j][i] for i in range(len(M))] for j in range(len(M))]
    
    I = [[0.00 for i in range(len(M))] for j in range(len(M))] 
    for i in range(len(M)):
        for j in range(len(M)):
            #creates an identity matrix
            I[j][j] = 1.00

    if M[1][1] == 0 and M[0][1] != 0:
        #swaps rows of first submatrices to prevent det=0
        swapRows(M, 0, 1, 4)

    partPivot(M, 0, len(M), len(M[0]))
    L_i, U_i = doolittle(M)
    
    y = [[0 for c in range(len(I))] for r in range(len(I))]
    for i in range(len(I)):
        for k in range (len(I[0])):
            y[i][k] = I[i][k]
            for j in range(i):
                y[i][k]=y[i][k]-(L_i[i][j]*y[j][k])
            y[i][k] = y[i][k]/L_i[i][i]
            
    n = len(y)
    x = [[0,0,0,0] for r in range(len(I))]
    if U_i[n-1][n-1] == 0: 
        #check for diagonal elements = 0
        raise ValueError
    
    for i in range(n-1, -1, -1):
        for k in range (len(I[0])): 
            x[i][k] = y[i][k]
            for j in range(i+1,n):
                x[i][k] = x[i][k] -(U_i[i][j]*x[j][k])
            x[i][k] = x[i][k]/U_i[i][i]
    
    x = transpose(x)
    unSwap(M, x, M_orig)
    x = transpose(x)
    
    return(x)

#Cholesky decomposition function set 
def choleskydecomp(A):
    #finds L using Cholesky's algorithm
    n = len(A) #null matrix for L
    L = [[0.0] * n for i in range(n)]
    for i in range(n):
        for k in range(i+1):
            tmp_sum = sum(L[i][j] * L[k][j] for j in range(k))
            if (i == k): # Diagonal elements
                L[i][k] = sqrt(A[i][i] - tmp_sum)
            else:
                L[i][k] = (1.0 / L[k][k] * (A[i][k] - tmp_sum))
    return L
    
def solveCholesky(L, U, b):
    n = len(L)
    y = [0 for i in range(n)]
    x = [0 for i in range(n)]
    for i in range(n):
        sumj = 0
        for j in range(i):
            sumj += L[i][j]*y[j]
        y[i] = (b[i]-sumj)/L[i][i]
    for i in range(n-1, -1, -1):
        sumj = 0
        for j in range(i+1, n):
            sumj += U[i][j]*x[j]
        x[i] = dec2Round((y[i]-sumj)/U[i][i])
    return x