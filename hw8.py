import numpy as np
import matplotlib.pyplot as plt

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

def main():

    # user specified fault tolerance
    Tau_g = 10**(-6) # some default for now
    # user specified step size tolerance
    Tau_f = 10**(-10) # some default for now

    # implement backtracking line search with parameters
    def backtracking_line_search(f, x0, p_k, alpha, rho, c=1e-4 , N_b = 128):
        '''
        Backtracking line search algorithm to find step size alpha_k
        Args:
            f: function to be minimized
            x0: current point
            p_k: search direction
            alpha: initial step size
            rho: reduction factor for step size
            c: constant for sufficient decrease condition, default is 10**(-4)
            N_b: maximum number of backtracking iterations, default is 128
        Returns:
            x_{k+1} : next point after taking step in direction p_k
            boolean flag, true if x_{k+1} which satisfies the Armijo condition is found, false otherwise
        '''

        x0 = np.asarray(x0)
        p_k = np.asarray(p_k)

        # epsilon = 1e-8
        epsilon = np.sqrt(np.finfo(float).eps)  # ~1.49e-8

        grad_f = np.zeros_like(x0)
        for i in range(len(x0)):
            ei = np.zeros_like(x0)
            ei[i] = epsilon
            grad_f[i] = (f(x0 + ei) - f(x0)) / epsilon

        directional_deriv = np.dot(grad_f, p_k)
        x_next = x0

        for i in range(N_b):
            x_next = x0 + alpha*p_k
            if (f(x_next) <= f(x0) + c*alpha*directional_deriv):
                return x_next, True
            else:
                alpha*=rho

        return x_next, False
    
    # different implementation, with gradient provided
    def bt_line_search(f, grad_f, x0, p_k, alpha, rho, c=1e-4 , N_b = 128):
        '''
        Backtracking line search algorithm to find step size alpha_k
        Args:
            f: function to be minimized
            grad_f: gradient of f
            x0: current point
            p_k: search direction
            alpha: initial step size
            rho: reduction factor for step size
            c: constant for sufficient decrease condition, default is 10**(-4)
            N_b: maximum number of backtracking iterations, default is 128
        Returns:
            x_{k+1} : next point after taking step in direction p_k
            boolean flag, true if x_{k+1} which satisfies the Armijo condition is found, false otherwise
        '''

        x0 = np.asarray(x0)
        p_k = np.asarray(p_k)
        grad = grad_f(x0)
        f_x0 = f(x0)
        directional_deriv = np.dot(grad, p_k)

        if directional_deriv >= 0:
            print("not in descent direction")
            return x0, False

        x_next = x0

        for i in range(N_b):
            x_next = x0 + alpha*p_k

            if (f(x_next) <= f_x0 + c*alpha*directional_deriv):
                return x_next, True
            
            alpha*=rho

        return x_next, False
    
    def plot(xs, grads, x_min, func: str):  
        '''
        generate plots for question 1 given vector sequence, gradient sequence, and the minimizing vector
        '''     
        ks = [i for i in range(len(xs))] # iterations 

        # Create the plot
        plt.figure(figsize=(10, 6))
        plt.plot(ks, np.log10(grads), label='log_10(grad(x_k))', linewidth=2)
        plt.plot(ks, np.log10([np.linalg.norm(x, ord=np.inf) for x in [xs[i] - x_min for i in ks]]), label='log_10(||x_k - A^{-1}F||_{inf})', linewidth=2)
        
        plt.xlabel('k values', fontsize=12)
        plt.title(f'{func}: k values versus log_10 of gradient and $||x_k - A^(-1)F||_inf', fontsize=14)
        plt.legend(fontsize=11)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
    
    def steepest_descent_backtracking(f, grad_f, x0, Tau_g, Tau_f, alpha, rho, N_sd, c=10**(-4), N_b=128):
        '''
        I chose to implement this with f and grad_f as separate functions
        
        Args:
            f: function to be minimized
            grad_f: gradient of the function to be minimized
            x0: initial point
            Tau_g: gradient tolerance for stopping criterion
            Tau_f: function value tolerance for stopping criterion
            alpha: initial step size for backtracking line search
            rho: reduction factor for step size in backtracking line search
            c: constant for sufficient decrease condition, default is 10**(-4)
            N_sd: maximum number of steepest descent iterations
            N_b: maximum number of backtracking iterations
        Returns:
            xs: np.array of points visited during optimization
            fs: np.array of function values at points visited during optimization
            grads: np.array of norm of gradients at points visited during optimization
            int flag, indicates stopping criterion met)
        '''
        x_cur = np.copy(x0)
        
        xs = [x_cur]
        fs = []
        grads = []
        
        for i in range(N_sd):
            f_cur = f(x_cur)
            grad_cur = grad_f(x_cur)
            grad_norm_cur = np.linalg.norm(grad_cur, ord=np.inf)
            
            fs.append(f_cur)
            grads.append(grad_norm_cur)
            
            # Check gradient tolerance
            if grad_norm_cur <= Tau_g:
                print("Gradient norm tolerance met.")
                return np.array(xs), np.array(fs), np.array(grads), 1
            
            # Search direction
            p_k = -grad_cur
            
            # Line search
            x_next, found = bt_line_search(f, grad_f, x_cur, p_k, alpha, rho, c, N_b)
            
            # Check if Armijo condition was satisfied
            if not found:
                xs.append(x_next)
                f_next = f(x_next)
                grad_next = grad_f(x_next)
                fs.append(f_next)
                grads.append(np.linalg.norm(grad_next, ord=np.inf))
                print("Armijo condition not fulfilled")
                return np.array(xs), np.array(fs), np.array(grads), 0
            
            # Compute next function value
            f_next = f(x_next)
            
            # Check function decrease tolerance
            if f_cur - f_next <= Tau_f:
                xs.append(x_next)
                grad_next = grad_f(x_next)
                fs.append(f_next)
                grads.append(np.linalg.norm(grad_next, ord=np.inf))
                print("Objective function value tolerance met.")
                return np.array(xs), np.array(fs), np.array(grads), 2
            
            # Update for next iteration
            xs.append(x_next)
            x_cur = x_next
        
        # Add final point's function and gradient
        f_cur = f(x_cur)
        grad_cur = grad_f(x_cur)
        fs.append(f_cur)
        grads.append(np.linalg.norm(grad_cur, ord=np.inf))
        
        print("Maximum number of iterations reached.")
        return np.array(xs), np.array(fs), np.array(grads), 3
    
    def q1_1():
        '''
        Steepest descent with mass spring system
        '''
        # dim
        N = 64
        
        A = np.array(spring_matrix(1, N + 1))

        identity_matrix = np.eye(N)

        e_31 = identity_matrix[30]
        e_32 = identity_matrix[31]

        F = e_31 - e_32

        x_min = np.linalg.solve(A, F)
        

        # quadratic function
        q = lambda x : 0.5 * x.T@A@x - F.T@x

        def grad_q(x):
            '''
            Calculate the gradient of q

            From the previous hw, we know that the gradient of q = 1/2 x.TAx - F.Tx should be Ax - F
            '''
            return A@x - F


        for alpha, rho in [
                           [1.0, 0.5], 
                           [2, 0.75], 
                           [0.5, 0.5], 
                           [1.0, 0.25]
                           ]:

            xs, fs, grads, flag = steepest_descent_backtracking(q, grad_q, np.zeros(64), 10**(-6), Tau_f, alpha, rho, 100000)

            print("Mass spring xs: ", xs)
            print("Mass spring fs: ", fs)
            print("Mass spring grads: ", grads)
            print("Mass spring flag: ", flag, "\n")

            plot(xs, grads, x_min, "q")


    def q1_2():
        '''
        Steepest descent with Rosenbrock function 
        '''
        def rosenbrock(x):
            """Rosenbrock function: f(x,y) = (1-x)^2 + 100(y-x^2)^2"""
            return (1 - x[0])**2 + 100*(x[1] - x[0]**2)**2
    
        def rosenbrock_grad(x):
            """Gradient of Rosenbrock function"""
            return np.array([
                -2*(1 - x[0]) - 400*x[0]*(x[1] - x[0]**2),
                200*(x[1] - x[0]**2)
            ])

        for x0 in [
            [1, 1], 
            [0, 0], 
            [-1, -1], 
            [0, -2], 
            [10, 10]
            ]:
            for alpha, rho in [
                [1, 0.5], 
                [2, 0.75], 
                [0.5, 0.5], 
                [1, 0.25]
                ]:

                x0 = np.array(x0) # chose some random starting vector 

                xs, fs, grads, flag = steepest_descent_backtracking(rosenbrock, rosenbrock_grad, x0, 10**(-6), Tau_f, alpha, rho, 100000)

                print("Rosenbrock xs: ", xs)
                print("Rosenbrock fs: ", fs)
                print("Rosenbrock grads: ", grads)
                print("Rosenbrock flag: ", flag, "\n")

                plot(xs, grads, np.array([1,1]), "Rosenbrock") # we know that [1,1] is the minimum from the previous hw


    #########################

    # start question 2

    # In your tests, take De = 1, re = 1/N, and a = 1, where N is the number of springs.
    def morse_potential_and_derivative(r, D_e, r_e, a):
        """
        Calculate the Morse potential energy and its derivative. (Draft
        version was generated by Claude! I’ll have you work with some
        more of Claude’s code next week.)
        Parameters:
        -----------

        r : float or array-like
            Distance between atoms
        D_e : float
            Dissociation energy (well depth)
        r_e : float
            Equilibrium bond distance
        a : float
            Parameter controlling the width of the potential well
            a = sqrt(k / (2*D_e)) where k is the force constant

        Returns:
        --------
        phi : float or array
            Morse potential energy
            phi(r) = D_e * (1 - exp(-a*(r - r_e)))^2
        
        dphi_dr : float or array
            Derivative of Morse potential with respect to r
            dphi/dr = 2 * D_e * a * (1 - exp(-a*(r - r_e))) * exp(-a*(r - r_e))
        """
        exp_term = np.exp(-a * (r - r_e))
        one_minus_exp = 1 - exp_term
        phi = D_e * one_minus_exp**2
        dphi_dr = 2 * D_e * a * one_minus_exp * exp_term
        
        return phi, dphi_dr

    def q2():
        '''
        Question 2
        '''

        '''
        We want the gradient of the term phi_ij at spot i or j

        Suppose we want to find spot i, can use the chain rule, and get dr_phi_ij * d(r_ij)/dx_i

        d(r_ij)/dx_i = d(np.abs(x[j] - x[i]))/dx_i 

        We know that the derivative of an absolute value function is thw sign of the function,

        so np.sign(x[j] - x[i]), then chain rule and multiply by d(-x[i])/dx[i] = -1, or d(x[j])/dx[j] = 1

        So dr_phi_ij * np.sign(x[j] - x[i]) * -1
        Or dr_phi_ij * np.sign(x[j] - x[i]) * 1
        '''

        def V(x, phi):
            N = len(x) - 1
            V_val = 0.0
            V_grad = np.zeros((N + 1,))
            for i in range(N):
                for j in range(i+1,N+1):
                    r_ij = np.abs(x[j] - x[i])
                    phi_ij, dr_phi_ij = phi(r_ij, 1, 1/N, 1) # match morse potential and derivative function
                    V_val += phi_ij
                    V_grad[j] += dr_phi_ij * np.sign(x[j] - x[i])
                    V_grad[i] += -dr_phi_ij * np.sign(x[j] - x[i])
                
            V_grad[0] = 0.
            V_grad[N] = 0.
            
            return V_val, V_grad
        
        # function from hw
        def random_test_vectors(N_masses, N_test_points):
            rng = np.random.default_rng(seed = 666)
            unit_vectors = rng.normal(size = (N_test_points, N_masses+1))
            unit_vectors[:, 0] = 0
            unit_vectors[:, N_masses] = 0
            lengths = np.linalg.norm(unit_vectors, axis = 1).reshape((N_test_points, 1))
            unit_vectors = unit_vectors/lengths
            x_evenly_spaced = np.linspace(0.,1.,N_masses+1)
            h = 1./(N_masses-1)
            xs = (rng.uniform(size = (N_test_points,N_masses+1))-0.5)*(h/16.)
            xs = xs + x_evenly_spaced
            xs[:,0] = 0.
            xs[:,N_masses] = 1.
            return unit_vectors, xs
        
        h = 1e-7 # tested to find this value
        unit_vecs, xs = random_test_vectors(65, 16) # 65 masses, get 16 vectors

        directional_dirs = np.zeros((16,))
        fin_approx_dirs = np.zeros((16,))
        
        for i in range(16):
            x = xs[i]
            unit_vec = unit_vecs[i]
            v_val, v_grad = V(x, morse_potential_and_derivative)

            directional_dirs[i] = np.dot(unit_vec, v_grad)
            fin_approx_dirs[i] = (V(x + h*unit_vec, morse_potential_and_derivative)[0] - v_val)/h # use forward difference
        
        print("finite approximation derivatives, forward difference: ", fin_approx_dirs)
        print("dir dirs: ", directional_dirs)

        print(np.max(np.abs(np.array([directional_dirs[i] - fin_approx_dirs[i] for i in range(16)]))))

    #########################

    # run code for hw
    q1_1()
    q1_2()
    q2()

if __name__ == "__main__":
    main()
