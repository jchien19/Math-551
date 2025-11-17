import numpy as np
from scipy.linalg import lu, solve_triangular
from matplotlib import pyplot as plt

# taken from given jupyter notebook
def spring_matrix(k,n):
    # make a zero array of shape (n-1,n-1)
    A=np.zeros((n-1,n-1)) #both sets of parentheses are required here!
    
    # set the zero'th row
    A[0,0]=2
    A[0,1]=-1
    
    # set rows 1 to n-3
    for i in range(1,n-2): # the n-2 is not a typo here...
        A[i,i-1]=-1
        A[i,i]=2
        A[i,i+1]=-1
        
    # set row n-1
    A[-1,-2]=-1
    A[-1,-1]=2
    
    return k*A

"""
Iteration has the form x_k = (I - PA)x_{k-1} + Pb

Can decompose A = L + D + U

Jacobi: let P = D^{-1}

Gauss-Seidel: let P = (L + D)^{-1} 
"""

def q6_1(A, b):
    '''
    question 6 part 1, PLU decompose matrix A, some system for x_tilde and find infinity norm backwards error
    '''
    P, L, U = lu(A) # decompose A

    y = solve_triangular(L, P.T@b, lower=True)

    x_tilde = solve_triangular(U, y)

    numerator_norm = np.linalg.norm(A@x_tilde - b, ord=np.inf)
    denominator_norm = np.linalg.norm(b, ord=np.inf)

    print("Infinity Norm backwards error for x_tilde: ", numerator_norm/denominator_norm)
    # the calculated value was: 1.1102665177675712e-15

    return x_tilde


def q6_iterative(A, b, x0, max_iter=64000):
    '''
    Perform jacobi and gauss-seidel iterative methods

    A - matrix A in Ax=b
    b - vector b in Ax=b
    n - dimensions of vectors and matrices
    max_iter = maximum iterations for both methods, set to 64000
    '''

    # get D inverse
    D_inv = np.diag(1.0 / np.diag(A))

    def jacobi(A, b, x0, max_iter):
        x_current = x0.copy()

        iterations = [] # store iterations

        for _ in range(max_iter):
            x_next = x_current + D_inv@b - D_inv@A@x_current

            iterations.append(np.copy(x_next))

            x_current = x_next

        return x_current, iterations
    
    def gauss_seidel(A, b, x0=None, max_iter=64000, tol=1e-10):
        '''
        (formula is in the book as well)
        Gauss Seidel is: x_{k+1} = (I - (L+D)^{-1}A)x_k + (L+D)^{-1}b

        = (I - (L+D)^{-1}(L + D + U))x_k + (L+D)^{-1}b
        = (I - ((L+D)^{-1}(L + D) + (L+D)^{-1}U)))x_k + (L+D)^{-1}b
        = (I - (I + (L+D)^{-1}U))x_k + (L+D)^{-1}b
        = -(L+D)^{-1}Ux_k + (L+D)^{-1}b
        = (L+D)^{-1}(b - Ux_k)

        => (L+D)x_{k=1} = b - Ux_k

        Where L + D is a lower triangular matrix
        '''

        # Store all iterations
        iterations = []
        
        # Decompose A = L + D + U
        L = np.tril(A, -1)  # strictly lower triangular
        D = np.diag(np.diag(A))  # Diagonal matrix
        U = np.triu(A, 1)  # strictly upper triangular

        
        # For (L+D)^{-1} in the iteration
        L_plus_D = L + D
        x_current = x0.copy()
        
        for _ in range(max_iter):
            rhs = -U@x_current + b

            x_next = solve_triangular(L_plus_D, rhs, lower=True)
            
            iterations.append(np.copy(x_next))
            
            x_current = x_next
        
        return x_current, iterations

    jacobi_res, jacobi_history = jacobi(A, b, x0, max_iter)
    print("\nFinished Jacobi Method\n")

    gauss_seidel_res, gauss_seidel_history = gauss_seidel(A, b, x0, max_iter)
    print("\nFinished Gauss Seidel Method\n")

    return jacobi_history, gauss_seidel_history

if __name__ == "__main__":
    iteration_nums = [_ for _ in range(1, 64001)]
    b = np.array([0.5*np.cos(i) for i in range(1,64)]) # initialize b vector
    A = spring_matrix(1, 64) # initialize matrix A
    x0 = np.arange(1, A.shape[0] + 1) / 64.0 # initialize starting vector
    
    # get x_tilde
    x_tilde = q6_1(A, b)

    # get list of vectors for jacobi and gauss-seidel
    jacobi_history, gauss_seidel_history = q6_iterative(A, b, x0)

    # plot
    jacobi_errors = [np.linalg.norm(x_j - x_tilde, ord=np.inf) for x_j in jacobi_history]
    gauss_seidel_errors = [np.linalg.norm(x_g - x_tilde, ord=np.inf) for x_g in gauss_seidel_history]

    # for part 3, estimate rate of convergence using last error points

    # jacobi
    print("Jacobi error: ", jacobi_errors[-1] / jacobi_errors[-2])
    # gauss_seidel
    print("Gauss-Seidel error: ", gauss_seidel_errors[-1] / gauss_seidel_errors[-2])

    # Take log10 of errs
    jacobi_log_errors = [np.log10(err) for err in jacobi_errors]
    gauss_seidel_log_errors = [np.log10(err) for err in gauss_seidel_errors]
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(iteration_nums, jacobi_log_errors, label='Jacobi', linewidth=2)
    plt.plot(iteration_nums, gauss_seidel_log_errors, label='Gauss-Seidel', linewidth=2)
    
    plt.xlabel('Iteration Number (m)', fontsize=12)
    plt.ylabel('log₁₀(||xₘ - x̃||∞)', fontsize=12)
    plt.title('Convergence of Jacobi and Gauss-Seidel Methods', fontsize=14)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    plt.show()
