import matplotlib.pyplot as plt
from typing import Dict
import numpy as np


def q8_1_1_dense(spring_consts: list[float]):
    '''dense matrix construction'''
    n = len(spring_consts)
    # make a zero array of shape (n-1,n-1)
    A=np.zeros((n-1,n-1)) #both sets of parentheses are required here!
    
    # set the zero'th row
    A[0,0]=spring_consts[0] + spring_consts[1]
    A[0,1]=-1*spring_consts[1]
    
    # set rows 1 to n-3
    for i in range(1,n-2): # the n-2 is not a typo here...
        A[i,i-1]=-1*spring_consts[i]
        A[i,i]=spring_consts[i] + spring_consts[i+1]
        A[i,i+1]=-1*spring_consts[i+1]
        
    # set row n-1
    A[-1,-2]=-1*spring_consts[n-2]
    A[-1,-1]=spring_consts[n-2] + spring_consts[n-1]
    
    return A

def q8_1_dok(spring_consts: list[float]):
    '''construct matrix in DOK format'''
    n = len(spring_consts)

    # for dictionary of keys
    A: Dict[(int, int), float] = {}
    A[(0,0)] = spring_consts[0] + spring_consts[1]
    A[(0,1)] = -1*spring_consts[1]


    # set rows 1 to n-3
    for i in range(1,n-2): # the n-2 is not a typo here...
        A[(i,i-1)]=-1*spring_consts[i]
        A[(i,i)]=spring_consts[i] + spring_consts[i+1]
        A[(i,i+1)]=-1*spring_consts[i+1]

    # set row n-1
    A[(n-1,n-2)]=-1*spring_consts[n-2]
    A[(n-1,n-1)]=spring_consts[n-2] + spring_consts[n-1]
    
    return A

def q8_1_csr(spring_consts: list[float]):
    '''return matrix in CSR format'''
    n = len(spring_consts)

    NC = n-1 # number of columns is the number of springs
    Rows = []
    Vals = []
    Cols = []

    # construct vals array

    Vals.append(spring_consts[0] + spring_consts[1])
    Cols.append(0)

    Vals.append(-1*spring_consts[1])
    Cols.append(1)

    Rows.append(0)

    # set rows 1 to n-3
    for i in range(1,n-2): # the n-2 is not a typo here...
        Vals.append(-1*spring_consts[i])
        Cols.append(i-1)

        Vals.append(spring_consts[i] + spring_consts[i+1])
        Cols.append(i)

        Vals.append(-1*spring_consts[i+1])
        Cols.append(i+1)

        Rows.append(len(Vals) - 3)
    
    Vals.append(-1*spring_consts[n-2])
    Cols.append(NC-2)

    Vals.append(spring_consts[n-2] + spring_consts[n-1])
    Cols.append(NC-1)

    Rows.append(len(Vals) - 2)

    Rows.append(len(Vals))

    return Rows, Cols, Vals, NC

def test_spring_matrix_construction():
    '''test csr construction'''

    # expect to get 
    # R = [0, 2, 5, 7]
    # C = [0, 1, 0, 1, 2, 1, 2]
    # V = [30.0, -20.0, -20.0, 50.0, -30.0, -30.0, 70.0]
    # NC = 3
    spring_consts_1 = [10.0, 20.0, 30.0, 40.0]

    Rows, Cols, Vals, NC = q8_1_csr(spring_consts_1)

    assert Rows == [0, 2, 5, 7], f"Rows incorrect: {Rows}"
    assert Cols == [0, 1, 0, 1, 2, 1, 2], f"Cols incorrect: {Cols}"
    assert Vals == [30.0, -20.0, -20.0, 50.0, -30.0, -30.0, 70.0], f"Vals incorrect: {Vals}"
    assert NC == 3, f"NC incorrect: {NC}"

    # expect to get
    # R = [0,2,4]
    # C = [0,1,0,1]
    # V = [2, -1, -1, 2]
    # NC = 2
    spring_consts_2 = [1, 1, 1]
    Rows, Cols, Vals, NC = q8_1_csr(spring_consts_2)
    assert Rows == [0,2,4], f"Rows incorrect: {Rows}"
    assert Cols == [0,1,0,1], f"Cols incorrect: {Cols}"
    assert Vals == [2, -1, -1, 2], f"Vals incorrect: {Vals}"
    assert NC == 2, f"NC incorrect: {NC}"

    # expect to get
    # R = [0,2,5,8,10]
    # C = [0, 1, 0, 1, 2, 1, 2, 3, 2, 3]
    # V = [9, -4, -4, 7, -3, -3, 5, -2, -2, 3]
    # NC = 4
    spring_consts_3 = [5, 4, 3, 2, 1]
    Rows, Cols, Vals, NC = q8_1_csr(spring_consts_3)
    assert Rows == [0,2,5,8,10], f"Rows incorrect: {Rows}"
    assert Cols == [0, 1, 0, 1, 2, 1, 2, 3, 2, 3], f"Cols incorrect: {Cols}"
    assert Vals == [9, -4, -4, 7, -3, -3, 5, -2, -2, 3], f"Vals incorrect: {Vals}"
    assert NC == 4, f"NC incorrect: {NC}"

    print("Spring matrix construction test passed.")

def q8_2(F: list[float], spr_consts: list[float], lengths: list[float]):
    '''construct vector B'''
    # ensure lengths are correct
    len_F = len(F)
    len_spr = len(spr_consts)
    len_lengths = len(lengths)

    if (len_F != len_spr-1):
        print("Length of F is not equal to length of spring constants minus 1")
        return

    if (len_F != len_lengths-1):
        print("Length of F is not equal to length of lengths minus 1")
        return
    
    if (len_spr != len_lengths):
        print("Length of spring constants is not equal to length of lengths")
        return

    B = [0 for _ in range(len_F)]

    for i in range(len_F):
        # F[i] gives external force with subscript i+1
        # spr_const[i] gives spring constant with subscript i
        # lengths[i] gives length with subscript i
        B[i] = F[i] + spr_consts[i]*lengths[i] - spr_consts[i+1]*lengths[i+1]
    
    return B

def test_vector_b_construction():
    # expect B = [-5, 5, 5, -90]
    B = q8_2([0, 0, 0, -100], [10, 20, 30, 40, 50], [1.5, 1, 0.5, 0.25, 0])
    assert B == [-5.0, 5.0, 5.0, -90.0], f"Vector B incorrect: {B}, should be [-5.0, 5.0, 5.0, -90.0]"

    # expect B = [-5, -5]
    B = q8_2([10, -10], [5, 15, 25], [3, 2, 1])
    assert B == [-5.0, -5.0], f"Vector B incorrect: {B}, should be [-5.0, -5.0]"

    # expect B = [-15, -22, -29]
    B = q8_2([1, 2, 3], [8, 12, 16, 20], [1, 2, 3, 4])
    assert B == [-15.0, -22.0, -29.0], f"Vector B incorrect: {B}, should be [-15.0, -22.0, -29.0]"

    print("Vector B construction test passed.")

def q8_3(A:list[list[float]]):
    # assume A is stored dense

    # compute infinity norm of A, ends up being max row sum of abs vals
    rows = []
    for row in A:
        rows.append(np.sum(np.abs(row)))

    return np.max(rows)

def q8_4():
    '''plot log condition number vs spring constant k4 = k12'''
    n = 16
    rs = [1/16.0 for _ in range(n)] # not sure if we need these if we're only constructing the A matrix?
    ks = [1.0 for _ in range(n)]

    condition_numbers = []
    k4_k12_values = []
    
    for _exp in range(1, 10):
        exp_val = 2**_exp
        ks[4] = exp_val
        ks[12] = exp_val
        
        # Instead of the toDense function, I directly constructed a dense matrix like in the examples from class
        A_dense = q8_1_1_dense(ks)
        
        # Compute condition number (infinity norm)
        A_inv = np.linalg.inv(A_dense)
        norm_A = np.linalg.norm(A_dense, ord=np.inf)
        norm_A_inv = np.linalg.norm(A_inv, ord=np.inf)
        cond_num = norm_A * norm_A_inv
        
        condition_numbers.append(cond_num)
        k4_k12_values.append(exp_val)
        

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(k4_k12_values, np.log(condition_numbers), 'bo-', linewidth=2, markersize=8)
    plt.xlabel('$k_4 = k_{12} = 2^{-k}$', fontsize=12)
    plt.ylabel('$\\log(\\kappa_\\infty(A))$', fontsize=12)
    plt.title('Logarithm of Condition Number vs. Spring Constant', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()



def q9(NC, V, R, C, vector):
    NR = len(R) - 1 # for a square matrix, NC = NR
    Mv = np.zeros(NR)
    
    for i in range(NR):
        start = R[i] # start of row i
        end = R[i + 1] # end of row i
        for j in range(start, end):
            Mv[i] += V[j] * vector[C[j]]
    
    return Mv

def test_q9():
    # matrix [[10, 0, 20],
    #         [0, 0, 30],
    #         [40, 50, 0]]
    # expect 70, 90, 140
    NC = 3
    V = [10, 20, 30, 40, 50]
    R = [0, 2, 3, 5]
    C = [0, 2, 2, 0, 1]

    vector = [1, 2, 3]

    result = q9(NC, V, R, C, vector)
    assert np.array_equal(result, np.array([70, 90, 140])), f"Result incorrect: {result}, should be [70, 90, 140]"

    # matrix [[1, 2],
    #         [0, 0],
    #         [3, 4]]
    # expect [17, 0, 39]
    NC = 2
    V = [1, 2, 3, 4]
    R = [0, 2, 2, 4]
    C = [0, 1, 0, 1]

    vector = [5, 6]
    result = q9(NC, V, R, C, vector)
    assert np.array_equal(result, np.array([17, 0, 39])), f"Result incorrect: {result}, should be [17, 0, 39]"

    # matrix [[7]]
    # expect [21]
    NC = 1
    V = [7]
    R = [0, 1]
    C = [0]

    vector = [4]
    result = q9(NC, V, R, C, vector)
    assert np.array_equal(result, np.array([28])), f"Result incorrect: {result}, should be [28]"

    #matrix [[0, 0, 0],
    #        [0, 5, 0],
    # expect [0, 15, 0]
    NC = 3
    V = [5]
    R = [0, 0, 1, 1]
    C = [1]
    vector = [1, 3, 4]
    result = q9(NC, V, R, C, vector)
    assert np.array_equal(result, np.array([0, 15, 0])), f"Result incorrect: {result}, should be [0, 15, 0]"

    print("Sparse matrix-vector multiplication test passed.")
    

if __name__ == "__main__":
    test_spring_matrix_construction()
    test_vector_b_construction()
    q8_4() #plot
    test_q9()
