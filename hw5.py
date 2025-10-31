import numpy as np
from scipy.linalg import solve_triangular
from matplotlib import pyplot as plt

## Question 1

def q1_2(n, k):
    ls = [_ for _ in range(1, n)]

    lamba_l = [2*k*(1-np.cos(l*np.pi/n)) for l in ls]

    # plot 
    ls = np.array(ls)

    plt.figure()
    plt.plot(ls, lamba_l, marker='o', linestyle='-')
    plt.xlabel('l')
    plt.ylabel('lambda_l')
    plt.title(f'lambda_l vs l (n={n}, k={k})')
    plt.grid(True)
    plt.show()

def inf_norm(A):
    '''infinity norm function from last hw'''
    # assume A is stored dense

    # compute infinity norm of A, ends up being max row sum of abs vals
    rows = []
    for row in A:
        rows.append(np.sum(np.abs(row)))

    return np.max(rows)

def fake_regression():
    '''fake regression problem from homework statement'''
    N=100
    k=14
    t=np.linspace(0,1,N,endpoint=True)
    X=np.zeros((N,k+1))
    for i in range(k+1):
        X[:,i]=t**i
        y=np.exp(np.sin(4*t))
        y=y/2006.787453080206
    return X,y

## Question 2

def q2():
    '''do qr factorization with fake regression function and calculate beta_hat value'''

    X, y = fake_regression()
    Q, R = np.linalg.qr(X, mode='reduced')
    Beta_hat = solve_triangular(R, Q.T@y, lower=False)
    # 1.0000001625366146 is the calculated value of beta hat_15
    return Beta_hat

def q3():
    '''question 3'''
    X, y = fake_regression()
    M = X.T@X
    Q, R = np.linalg.qr(M, mode='reduced')
    Beta_hat = solve_triangular(R, Q.T@X.T@y, lower=False)
    # -0.9832027398277714 is the calculated value for beta hat_15, which is very far from 1
    # considering machine precision

    # condition number with infinity norm
    M_inv = np.linalg.inv(M)
    norm_M = inf_norm(M)
    norm_M_inv = inf_norm(M_inv)
    cond_num = norm_M * norm_M_inv

    print("condition number: ", cond_num)
    # I calculated 5.085812225290763e+18 for my condition number

    # calculate actual error using results from the question 2, use infinity norm
    Beta_hat_actual = q2()

    norm_actual_minus_approx = np.linalg.norm(Beta_hat_actual - Beta_hat, ord=np.inf)
    norm_actual = np.linalg.norm(Beta_hat_actual, ord=np.inf)

    print("Beta hat relative error: ", norm_actual_minus_approx/norm_actual)

    ''' 
    The calculated relative error is 1.9979037115805112

    machine error for doubles is about 2.22e-16 while the calculated condition number is about 5.09e+18
    If we multiply machine error with the calculated condition number, we get about 1.13e+3.

    So our calculated relative error is much lower than the estimate of the error.
    '''

    return Beta_hat


if __name__ == "__main__":
    print("q2 Beta hat 15: ", q2()[14]) # print 15th value
    print("q3 Beta hat 15: ", q3()[14])

    ## plot
    q1_2(10, 50)
