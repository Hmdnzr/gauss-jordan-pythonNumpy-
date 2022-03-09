
from fractions import Fraction as f
import numpy as np

def Upper_Triangle_Solver(U, b, n):
    x = np.zeros((n,1), dtype=int)
    x = x + f()
    for i in range(n-1, -1, -1):
        rhs = f(b[i], U[i][i])
        for j in range(n-1, i, -1):
            rhs -= U[i][j]*x[j]/U[i][i]
        x[i] = rhs
    return x

def Gauss_Elimination_Method(A, b, n, Jordan = 0):
    print('\nGiven matrix A:\n', A, '\nb:\n', b, '\n')
    print('---------------------------------------------------------\n')
    #Note: If any diagonal element is zero, algorithm won't work. Rearrange columns of matrix A and corresponding terms in x and b vectors.
    #^Potential Improvement. 
    
    if(Jordan != 0):      
        x = np.zeros((n,1), dtype=int)
        x = x + f()
        invA = np.zeros((n,n), dtype=int)
        invA = invA + f()
        for i in range(n):
            invA[i][i] = f(1,1)
        for k in range(n-1):
            for i in range(k+1, n):
                coeff = f(A[i][k], A[k][k])
                b[i] = b[i] - b[k]*coeff
                print('R', i+1, '-> R', i+1, '- R', k+1, '* (', A[i][k], '/', A[k][k], ')')
                for j in range(k, n):
                    A[i][j] = A[i][j] - A[k][j]*coeff 
                    invA[i][j] = invA[i][j] - invA[k][j]*coeff
                
                print('\nA:\n', A, '\nb:\n', b, '\n')
                print('Inverted Matrix steps:\n', invA, '\n')
                print('---------------------------------------------------------\n')
        for k in range(n-1, 0, -1):
            for i in range(k-1, -1, -1):
                coeff = f(A[i][k], A[k][k])
                b[i] = b[i] - b[k]*coeff
                print('R', i+1, '-> R', i+1, '- R', k+1, '* (', A[i][k], '/', A[k][k], ')')
                for j in range(k, n):
                    A[i][j] = A[i][j] - A[k][j]*coeff
                    invA[i][j] = invA[i][j] - invA[k][j]*coeff
                print('\nA:\n', A, '\nb:\n', b, '\n')
                print('Inverted Matrix steps:\n', invA, '\n')
                print('---------------------------------------------------------\n')
        for k in range(n):
            b[k] = b[k] / A[k][k]
            invA[k][k] = invA[k][k] / A[k][k]
            print('R', k+1, '->', 'R', k+1, '/ (', A[k][k], ')')
            A[k][k] = f(1,1)
            print('\nA:\n', A, '\nb:\n', b, '\n')
            print('Inverted Matrix steps:\n', invA, '\n')
            print('---------------------------------------------------------\n')
        x = b

    else:
        for k in range(n-1):
            for i in range(k+1, n):
                coeff = f(A[i][k], A[k][k])
                b[i] = b[i] - b[k]*coeff
                print('R', i+1, '-> R', i+1, '- R', k+1, '* (', A[i][k], '/', A[k][k], ')')
                for j in range(k, n):
                    A[i][j] = A[i][j] - A[k][j]*coeff
                print('\nA:\n', A, '\nb:\n', b, '\n')
                print('---------------------------------------------------------\n')     
    print('Final Solution Vector X:\n', Upper_Triangle_Solver(A, b, n),'\n')

def LU_Decomposition(A, b, n, type = 'Crotch'):
    stepNumber = 0
    for i in range(n):
        stepNumber += 1
        print('\nStep', stepNumber, ':\n')
        for j in range(i, n):
            sum = 0
            print('l', j+1, i+1, '=', A[j][i], end='')
            for k in range(i):
                sum += A[k][i]*A[j][k]
                print(' - (', A[k][i], '*', A[j][k], end=')')
            A[j][i] = A[j][i] - sum
            print(' =', A[j][i])
        stepNumber += 1
        print('\nStep', stepNumber, ':\n')
        for j in range(i+1, n):
            sum = 0
            print('u', i+1, j+1, '= (', A[i][j], end='')
            for k in range(i):
                sum += A[k][j]*A[i][k]
                print(' - (', A[k][j], '*', A[i][k], end=')')
            A[i][j] = (A[i][j] - sum) / A[i][i]
            print(') /', A[i][i], ' =', A[i][j])
        print('y', i+1, '= (', b[i], end='' )
        sumRHS = 0
        for k in range(i):
            sumRHS += b[k]*A[i][k]
            print(' - (', b[k], '*', A[i][k], end=')')
        b[i] = (b[i] - sumRHS) / A[i][i]
        print(') /', A[i][i], ' =', b[i])
    print('\n---------------------------------------------------------')
    print('\nFinal D Matrix: \n', A)
    print('\nFinal Y Matrix: \n', b)
    for i in range(n):
        A[i][i] = 1
    print('\n---------------------------------------------------------\n')
    print('Final Solution Vector X:\n', Upper_Triangle_Solver(A, b, n),'\n')
    
def norm(y, x, n):
    squareSum = 0
    for i in range(n):
        squareSum += ((y[i] - x[i])**2)
    return (squareSum**0.5)

def Gauss_Jacobi_Iteration(A, b, n, error = 0.01):
    x = np.zeros((n,1))
    y = np.zeros((n,1))
    iterationNumber = 0
    for i in range(15):
        print('\nIteration:', iterationNumber+1)
        iterationNumber += 1
        for j in range(n):
            sum = 0
            for k in range(n):
                if(k==j):
                    continue
                sum += x[k]*A[j][k]
            y[j] = (b[j] - sum) / A[j][j]
            print('x', j+1, '=', y[j])
        currentError = norm(y, x, n) 
        if(currentError < error):
            print('\nConverged as norm is less than', error)
            break
        x = y.copy()
        
def Gauss_Seidal_Iteration(A, b, n, error = 0.01):
    x = np.zeros((n,1))
    y = np.zeros((n,1))
    t = np.zeros((n,1))
    iterationNumber = 0
    for i in range(15):
        print('\nIteration:', iterationNumber+1)
        iterationNumber += 1
        for j in range(n):
            sum = 0
            for k in range(n):
                if(k==j):
                    continue
                sum += x[k]*A[j][k]
            y[j] = (b[j] - sum) / A[j][j]
            x[j] = y[j]
            print('x', j+1, '=', y[j])
        currentError = norm(y, t, n) 
        if(currentError < error):
            print('\nConverged as norm is less than', error)
            break
        t = x.copy()
        
def Eigen_Estimation(A, n):
    v = np.zeros((n,1))
    v = v + 1.0
    m = [None, None]
    print('A * v = \n')
    print(A, '*', v, '=\n')
    v = np.dot(A, v)
    m[0] = max(abs(v))
    print(v, '= Y 1\n')
    print('m 1 =', m[0])
    v = v / m[0]
    print('v 1 = Y 1 / m 1')
    print('\n---------------------------------------------------------\n')
    print('A * v 1 = \n')
    print(A, '*', v, '=\n')
    v = np.dot(A, v)
    m[1] = max(abs(v))
    print(v, '= Y 2\n')
    print('m 2 =', m[1])
    v = v / m[1]
    print('v 2 = Y 2 / m 2')
    print('\n---------------------------------------------------------\n')
    i = 2
    while(abs(abs(m[1])-abs(m[0])) > 0.01):
        v = np.dot(A, v)
        print('A * v', i, ' = \n')
        print(A, '*', v, '=\n')
        m[0] = m[1]
        m[1] = max(abs(v))
        print(v, '= Y', i, '\n')
        print('m', i, '=', m[1])
        v = v / m[1]
        print('v', i, '= Y', i, '/ m', i)
        print('\n---------------------------------------------------------\n')
        i+=1

a = np.array([[f(1,1), f(2,1)], [f(3,1), f(4,1)]])

b = np.array([f(1,3), f(3,5)])

C = np.array([[f(3), f(2), f(4)],
              [f(-1), f(4), f(10)], 
              [f(1), f(3), f(-1)]])

d = np.array([f(-5,1), f(-12,1), f(4,1)])

D = np.array([[f(1,1), f(1,1), f(1,1), f(1,1)],
              [f(2,1), f(3,1), f(4,1), f(5,1)], 
              [f(3,1), f(-1,1), f(1,1), f(1,1)], 
              [f(1,1), f(-1,1), f(3,1), f(5,1)]])

e = np.array([f(4,1), f(14,1), f(4,1), f(8,1)])

#Gauss_Elimination_Method(a, b, 2)
Gauss_Elimination_Method(D, e, 4, 0)
#LU_Decomposition(D, e, 4)
#Gauss_Jacobi_Iteration(E, g, 3)
#Gauss_Seidal_Iteration(E, g, 3)
#Eigen_Estimation(C, 3)