import matplotlib.pyplot as plt
import numpy as np

import matplotlib.cm as cm

def q2_1():

    def r(x: np.array):
        A = np.array([[1, 0], [0, 16]])
        return 0.5 * x.T@A@x
    
    y1 = np.linspace(-3, 3, 400)
    y2 = np.linspace(-3, 3, 400)
    Y1, Y2 = np.meshgrid(y1, y2)
    YY = np.stack([Y1.ravel(), Y2.ravel()], axis=1)
    Rvals = np.array([r(y) for y in YY]).reshape(Y1.shape)

    fig1, ax1 = plt.subplots(figsize=(6, 6))
    levels = [0.5, 1.0, 2.0, 4.0, 8.0]
    cs1 = ax1.contour(Y1, Y2, Rvals, levels=levels)
    ax1.set_aspect('equal', 'box')
    ax1.set_title("r(x) = (1/2)x.T@A@x")
    ax1.set_xlabel('x coord')
    ax1.set_ylabel('y coord')
    ax1.clabel(cs1, inline=True, fontsize=8)
    ax1.grid(True)

    plt.show()

def q2_4():

    # --- Parameters ---
    theta = np.pi / 6   # rotation angle = 30 degrees
    Q = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta),  np.cos(theta)]
    ])
    L = np.diag([1.0, 16.0])   # eigenvalues (positive definite)
    A = Q @ L @ Q.T            # A = Q L Q^T (symmetric positive definite)

    b = np.array([2.0, 1.0])   # you can change this
    x_star = np.linalg.solve(A, b)  # minimizer A^{-1} b

    # --- Define the functions ---
    def r_of_y(y):
        return 0.5 * y.T @ L @ y

    def q_of_x(x):
        return 0.5 * x.T @ A @ x - b.T @ x

    # --- Plot r(y): axis-aligned ellipses in y-plane ---
    y1 = np.linspace(-3, 3, 400)
    y2 = np.linspace(-3, 3, 400)
    Y1, Y2 = np.meshgrid(y1, y2)
    YY = np.stack([Y1.ravel(), Y2.ravel()], axis=1)
    Rvals = np.array([r_of_y(y) for y in YY]).reshape(Y1.shape)

    fig1, ax1 = plt.subplots(figsize=(6, 6))
    levels = [0.5, 1.0, 2.0, 4.0, 8.0]
    cs1 = ax1.contour(Y1, Y2, Rvals, levels=levels)
    ax1.set_aspect('equal', 'box')
    ax1.set_title("r(y) = 1/2 * y^T L y  (y = Q^T(x - A^{-1}b))")
    ax1.set_xlabel('y₁ (rotated coord)')
    ax1.set_ylabel('y₂ (rotated coord)')
    ax1.clabel(cs1, inline=True, fontsize=8)
    ax1.grid(True)

    # --- Plot q(x): rotated ellipses centered at x_star ---
    x1 = np.linspace(x_star[0] - 3, x_star[0] + 3, 400)
    x2 = np.linspace(x_star[1] - 3, x_star[1] + 3, 400)
    X1, X2 = np.meshgrid(x1, x2)
    XX = np.stack([X1.ravel(), X2.ravel()], axis=1)

    # Compute y = Q^T(x - x_star)
    Y = (Q.T @ (XX.T - x_star.reshape(2, 1))).T
    Qvals = np.array([0.5 * y.T @ L @ y for y in Y]).reshape(X1.shape)

    fig2, ax2 = plt.subplots(figsize=(6, 6))
    cs2 = ax2.contour(X1, X2, Qvals, levels=[0.5, 1.0, 2.0, 4.0, 8.0])
    ax2.set_aspect('equal', 'box')
    ax2.set_title("q(x) = 1/2 x^T A x - b^T x")
    ax2.set_xlabel('x₁')
    ax2.set_ylabel('x₂')
    ax2.plot(x_star[0], x_star[1], 'rx', label='x* = A⁻¹b')
    ax2.legend()
    ax2.clabel(cs2, inline=True, fontsize=8)
    ax2.grid(True)

    plt.show()

    # --- Print key results ---
    print("A =\n", A)
    print("Q =\n", Q)
    print("L =\n", L)
    print("b =", b)
    print("x* = A^{-1}b =", x_star)



def q3(x_range, y_range):
    '''
    make contour plot of Rosenbrock function
    '''

    # delta = 0.025
    # x = np.arange(-x_range, x_range, delta)
    # y = np.arange(-y_range, y_range, delta)
    # X, Y = np.meshgrid(x, y)

    # # Rosenbrock function 
    # Z = (1 - X)**2 + 100*(Y-X**2)**2

    # fig, ax = plt.subplots()
    # CS = ax.contour(X, Y, Z)
    # ax.clabel(CS, fontsize=10)
    # ax.set_title('Rosenbrock Function Contour')

    # plt.show()

    b = np.array([1, 1])
    theta = np.pi/2.0
    Q = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    Lambda = np.array([[1, 0],[0, 16]])


    delta = 0.025
    x = np.arange(-x_range, x_range, delta)
    y = np.arange(-y_range, y_range, delta)
    X, Y = np.meshgrid(x, y)

    # Rosenbrock function 
    Z = 0.5 * np.array([X, Y]).T@Q@Lambda@Q.T@np.array([X, Y]) - b.T@np.array([X, Y])

    fig, ax = plt.subplots()
    CS = ax.contour(X, Y, Z)
    ax.clabel(CS, fontsize=10)
    ax.set_title('Rosenbrock Function Contour')

    plt.show()

def rosenbrock(x, y) -> float:
    '''
    Calculate the rosenbrock function
    x: float
    y: float
    '''
    return (1.0 - x)**2.0 + 100.0*(y - x**2.0)**2.0

def rosenbrock_grad(X: np.array) -> np.array:
    '''
    Calculate the gradient of the rosenbrock function
    x: float
    y: float
    '''
    x = X[0]
    y = X[1]
    return np.array([-2.0*(1-x) - 400*x*(y-x**2), 200*(y-x**2)])

def inverse_Hessian(X: np.array) -> np.array:
    '''
    We can calculate the inverse of the Hessian matrix:

    (D^2f)^{-1} = 
    det(D^2f)^{-1} * 
    [[200, 400x], 
    [400x ,1200x^2 - 400y + 2]]

    note that we found that the det(D^2f) = 0 if on the equation y = x^2 + 1/200
    '''

    x = X[0]
    y = X[1]

    det = 1/(400 * (1 + 200*x**2 - 200*y)) # calculated from part 5

    a = 200
    b_c = 400*x
    d = 1200*x**2 - 400*y + 2

    inv_hessian = np.array([[a, b_c], [b_c, d]])

    return det * inv_hessian

def q6(iters: int):

    def newton_method(iters: int, x0: np.array):
        '''
        Newtons method for Rosenbrock function
        
        From class, for the Newton's Method we have:

        x_{n+1} = x_n - f''(x,y)^{-1}grad(f(x_n))
        '''

        x_cur = np.copy(x0)
        x_new = None

        res = [None for _ in range(iters + 1)]
        res[0] = x_cur

        for i in range(1, iters+1):
            x_new = x_cur - inverse_Hessian(x_cur)@rosenbrock_grad(x_cur)
            res[i] = np.copy(x_new)
            x_cur = x_new

        return res
    
    def plot(iters, Xs, plot):
        if plot:
            iters = [num for num in range(iters + 1)]
            x_coords = [v[0] for v in Xs]
            y_coords = [v[1] for v in Xs]

            # Create the plot
            plt.figure(figsize=(10, 6))
            plt.plot(iters, x_coords, label='x-coord', linewidth=2, color='blue')
            plt.plot(iters, y_coords, label='y-coord', linewidth=2, color='red')
            
            plt.xlabel('Iteration Number', fontsize=12)
            plt.ylabel('Newton Method Result', fontsize=12)
            plt.title(f'Iteration Number vs Newton Method, X0 = [{Xs[0][0]}, {Xs[0][1]}]', fontsize=14)
            plt.legend(fontsize=11)
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            
            plt.show()
    
    # test different starting vectors
    starts = [
        [2, 2],
        [1, 1],
        [3, -3],
        [1, 2],
        [0.5, 0.5],
        [1, -16],
        [100000, -100000]
    ]
    
    # loop to test different starting values
    for V in starts:
        res = newton_method(iters, np.array(V))
        # print("V, res: ", V, res, "\n")

        plot(iters, res, V == [3, -3])

if __name__ == "__main__":
    q2_1()
    q2_4()
    q3(5.0, 5.0)
    q6(10)
