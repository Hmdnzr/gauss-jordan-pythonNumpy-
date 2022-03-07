import numpy as np
from sympy import Matrix, Rational, mod_inverse, pprint
import pprint
import time


class P_matrix:
    """ A class containing function that perform operations on a matrix in field of p.
    The class takes the field and the matrix as an input at the point of instantiation"""
    def __init__(self,p,matrix) :
        self._p = p
        self._matrix = matrix
        self._inverse = list(range(self._p))
        self.calc_inverse()

    def get_p(self):
        return self._p

    def get_matrix(self):
        return self._matrix

    def add_row(self, row_1, row_2, m = 1):
        

        matrix_A = self.get_matrix() # matrix is a list of a list
        first_row = matrix_A[row_1]
        second_row = matrix_A[row_2]
        for i in range(len(first_row)):
            second_row[i] = (m * first_row[i] + second_row[i]) % self.get_p()      
        return 

    def get_inverse(self):
        return self._inverse

    def calc_inverse(self):
        p_val = self.get_p()
        for i in range(2,p_val):
            for j in range(2, p_val):
                if i*j % p_val == 1 :
                    self._inverse[i] = j
        return


    def switch_row(self, row_1, row_2):
        matrix_A = self.get_matrix()
        matrix_A[row_1], matrix_A[row_2] = matrix_A[row_2], matrix_A[row_1]
        return 

    def multiply_row(self, row, m): # takes row index and the scaler multiplier
        matrix_A = self.get_matrix()
        for i in range(len(matrix_A[row])):
            matrix_A[row][i] = (m * matrix_A[row][i]) % self.get_p() 
        return


def activation_matrix(m,n):
    """  mario definisi dari mar=trix dibawah ini """ 

    matrix_A = np.identity(m*n) # generate an m*n identity matrix
    for p in range(n):
        for q in range(m):
            i = n*q + p   # p controls movement along puzzle column and q controls movement 
            # along puzzle row. i is the matrix row number 

            if q > 0:   # To fill out the identity matrix below the diagonal
              matrix_A[i][i-n] = 1

            if q < m-1: # To fill out the identity matrix above the diagonal
              matrix_A[i][i+n] = 1

            if p > 0:   # To cover 1's right below the diagonal
              matrix_A[i][i-1] = 1

            if p < n-1: # To cover 1's right above the diagonal
              matrix_A[i][i+1] = 1

    matrix_A = matrix_A.astype(int) # converting numpy float array to numpy int array
    matrix_A = matrix_A.tolist() # converting numpy array to list

    return matrix_A


def build_matrix_for_rref(m,n):
    """ This function takes inputs of m-number of rows and n-number of columns and forms a 
    matrix for computing the Gaus-Jordan elimination of a lights out puzzle whose lights 
    start all off and to be all switched on
    
    It returns the matrix as a list of columns with an element of the solution vector at the end"""

    soln_state_vec = [[1]]*m*n # Create an m*n long vector of all ones

    mat_for_rref = np.append(activation_matrix(m,n), soln_state_vec, axis=1) # Attach the solution state vector to the end of the matrix

    return mat_for_rref.tolist()
    

def REF(matrix, m,n):
    """ A function that computes the Row Echelon Form of an m*n matrix.
    It take in an input of the matrix whose last column are the elements of the solution state vector
    and returns the matrix in an REF as a list of columns"""

    pivot_columns=[] # specifies the columns of the pivots(the leading variable of a row).
    # NB: the length of the pivot columns equals the number of rows with pivots

    for j in range((m*n)): # for each column in the matrix
        for i in range(len(pivot_columns), m*n): # for each row beginning from the row whose pivot is to be found
            if matrix.get_matrix()[i][j] != 0: # if the element is a non zero,
                # i is the row in column j with the leading variable.

                matrix.switch_row(i,len(pivot_columns)) # switch matrix[i] to the row whose leading variable is 
                # to be determined(number of pivot columns recorded thus far)

                matrix.multiply_row(len(pivot_columns),matrix.get_inverse()[matrix.get_matrix()[len(pivot_columns)][j]]) # multiply the row
                # with the leading variable by the inverse of the leading varible to make it one.

                pivot_columns.append(j) # store pivot in pivot column

                break # move to the next for loop

        for i in range(len(pivot_columns),m*n): # for all row below that of the leading variable
            if matrix.get_matrix()[i][j] != 0: # if you find a non-zero
                
                matrix.add_row(len(pivot_columns)-1, i, -matrix.get_matrix()[i][j]) # increase leading row by a factor
                # the non-zero variable and add leading row to current row in mod p.
    return np.array(matrix.get_matrix()).astype(int).tolist() , pivot_columns



def RREF(tuple_ref_pivot_column,p):
    piv_col = tuple_ref_pivot_column[1]
    n_i = len(piv_col)
    mat = tuple_ref_pivot_column[0]
    obj_matrix = P_matrix(p,mat)
    for j in reversed(range(len(piv_col))):
        for i in reversed(range(n_i-1)): 
            if obj_matrix.get_matrix()[i][j] != 0: # if you find a non-zero 
                obj_matrix.add_row(n_i-1, i, -obj_matrix.get_matrix()[i][j]) # increase leading row by a factor
                # the non-zero variable and add leading row to current row in mod p.
        n_i -= 1
    return np.array(obj_matrix.get_matrix()).astype(int).tolist() , piv_col


def sol_no_sol(matrix, m, n):
    """
    returns the |sollution space| if there is a solution and 
    "NO SOLUTION" if there is not
    """
    piv_column = matrix[1]
    rref_matrix = matrix[0]
    len_piv_columns = len(piv_column)
    # print("length of pivot column: ", len_piv_columns)
    col_dim = m*n
    nullity = col_dim - len_piv_columns    
    for i in range(len_piv_columns, nullity+len_piv_columns+1):
        if i == m*n:
            return nullity
        elif rref_matrix[i][-1] != 0:
            return 'NO SOLUTIONS'
        else:
            return nullity
        
def bool_sol_no_sol(matrix, m, n):
    if sol_no_sol(matrix, m, n) == 'NO SOLUTIONS':
        return 0
    else:
        return 1
        

def nullity(matrix, m, n):
    piv_column = matrix[1]
    rref_matrix = matrix[0]
    len_piv_columns = len(piv_column)
    col_dim = m*n
    nullity = col_dim - len_piv_columns   
    return nullity




#helper function
if __name__=="__main__":
    a = 2
    b = 4
    p = 3
    my_matrix = build_matrix_for_rref(a,b)
    pp = pprint.PrettyPrinter(depth=6)
    A = P_matrix(p, my_matrix)
    #pp.pprint(A.get_matrix())
    print()
    # A.add_row(0, 1)
    #pp.pprint(A.get_matrix())
    # print()
    # A.switch_row(2,3)
    # pp.pprint(A.get_matrix())
    # print()
    # A.multiply_row(0,3)
    #pp.pprint(A.get_matrix())
    start_time = time.time()
    ref=REF(A,a,b)
    print("It took ", time.time() - start_time, "seconds.")
    #pp.pprint(ref)
    print()
    rref = RREF(ref,p)
    pp.pprint(rref)
    print()
    pp.pprint(nullity(rref, a, b))