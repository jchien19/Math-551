import numpy as np
import matplotlib.pyplot as plt

def q1_and_3_main():

    # problem defined values
    k = 1
    m = 1

    v_0 = 0
    q_0 = 1
    
    def euler_method_iteration(q, v, h):
        '''
        Compute one iteration of the euler method

        :param q: position at the start of euler's method iteration
        :param v: velocity at the start of euler's method iteration
        :param h: step size
        '''

        # gradient of q is v
        q_next = q + v*h
        v_next = v - (h* k * q / m)

        return q_next, v_next
    
    def euler_method(q0, v0):

        q_n = q0
        v_n = v0

        qs = [q0]
        vs = [v0]

        for _ in range(round(30/0.02)):
            q_n, v_n = euler_method_iteration(q_n, v_n, 0.02)
            # print(q_n)
            # print(v_n)
            qs.append(q_n)
            vs.append(v_n)

        return qs, vs
    
    def Stromer_Verlet_method_iteration(x, v, v_h, h):
        '''
        compute one iteration of the stromer verlet method

        Formulas
        v_n+1/2 = v_n - 1/2 delta_t M^{-1} gradV(x_n)
        x_n+1 = x_n + delta_t v_n+1/2
        v_n+1 = v_n+1/2 - 1/2 delta_t M^{-1}gradV(x_n+1)
        '''
        
        # returns gradient of V
        grad_V = lambda x: k*x

        v_h = v - 0.5 * h * 1/M * grad_V(x)
        x_next = x + h * v_h
        v_next = v_h - 0.5 * h * 1/M * grad_V(x_next)
        
        return v_next, x_next, v_h
    
    def Stromer_Verlet_method(x0, v0, h):

        v_n = v0
        x_n = x0
        v_h = None

        xs = [x0]
        vs = [v0]

        for _ in range(round(30/0.02)):
            v_n, x_n, v_h = Stromer_Verlet_method_iteration(x_n, v_n, v_h, h)
            xs.append(x_n)
            vs.append(v_n)

        return xs, vs
    
    def plot(xs, vs, method):
        plt.figure(figsize=(10, 6))
        plt.plot(xs, vs, linewidth=2, label="position")
        
        plt.xlabel('x values', fontsize=12)
        plt.ylabel('v values', fontsize=12)
        plt.title(f'{method} Method', fontsize=14)
        plt.legend(fontsize=11)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
    
    def plot_hamiltionian(xs, vs, method):
        # Hamiltionian function for our specific case
        hamiltonian = lambda x, v: 0.5*x**2 + 0.5*v**2

        num_of_iters = len(xs)
        time = [_ for _ in range(num_of_iters)]
        hamiltionian_vals = [hamiltonian(xs[i], vs[i]) for i in range(num_of_iters)]

        # Create the plot
        plt.figure(figsize=(10, 6))
        plt.plot(time, hamiltionian_vals, linewidth=2, label="Hamiltonian")
        
        plt.xlabel('time values', fontsize=12)
        plt.ylabel('ham values', fontsize=12)
        plt.title(f'Hamiltonian vs Time for {method} Method', fontsize=14)
        plt.legend(fontsize=11)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    M = 1 # generally, this should be a matrix, but for only one mass we can let it be scalar

    # Euler
    xs, vs = euler_method(q_0, v_0)

    # Create the plot
    plot(xs, vs, "Euler")
    plot_hamiltionian(xs, vs, "Euler")

    # Stromer Verlet
    xs, vs = Stromer_Verlet_method(1, 0, 0.02)
    # Create the plot
    plot(xs, vs, "Stromer Verlet")
    plot_hamiltionian(xs, vs, "Stromer Verlet")

# ============================================================================
# PART 1: MODEL AND RESIDUAL FUNCTIONS
# ============================================================================

def morse_potential(r, params):
    """
    Compute Morse potential: V(r) = D_e * (1 - exp(-a*(r - r_e)))^2
    
    Parameters:
    -----------
    r : array, shape (n,)
        Distances
    params : array, shape (3,)
        Parameters [D_e, r_e, a]
    
    Returns:
    --------
    V : array, shape (n,)
        Potential energies
    """
    D_e, r_e, a = params
    exp_term = np.exp(-a * (r - r_e))
    return D_e * (1 - exp_term)**2


def residual(params, r_data, V_data):
    """
    Compute residual vector: r_i = V_model(r_i; params) - V_data_i
    
    Parameters:
    -----------
    params : array, shape (3,)
        Parameters [D_e, r_e, a]
    r_data : array, shape (n,)
        Distance measurements
    V_data : array, shape (n,)
        Energy measurements
    
    Returns:
    --------
    res : array, shape (n,)
        Residuals
    """
    V_model = morse_potential(r_data, params)
    return V_model - V_data


def objective_function(params, r_data, V_data):
    """
    Compute objective function: f(θ) = (1/2) * ||r(θ)||^2
    
    This is the sum of squared residuals (SSR).
    """
    res = residual(params, r_data, V_data)
    return 0.5 * np.sum(res**2)


# ============================================================================
# PART 2: JACOBIAN COMPUTATION 
# ============================================================================

def jacobian(params, r_data):
    """
    Compute the Jacobian matrix of the residual vector.
    
    J[i, j] = ∂r_i / ∂θ_j = ∂V_model(r_i; θ) / ∂θ_j
    
    DERIVATION:
    -----------
    Starting with V(r) = D_e * (1 - exp(-a*(r - r_e)))^2
    
    Let u = 1 - exp(-a*(r - r_e))
    Then V = D_e * u^2
    
    1. ∂V/∂D_e = u^2 = (1 - exp(-a*(r - r_e)))^2
    
    2. ∂V/∂r_e:
       ∂u/∂r_e = -a * exp(-a*(r - r_e))
       ∂V/∂r_e = 2*D_e*u*(∂u/∂r_e) = 2*D_e*a*(1 - exp(-a*(r - r_e)))*exp(-a*(r - r_e))
    
    3. ∂V/∂a:
       ∂u/∂a = (r - r_e) * exp(-a*(r - r_e))
       ∂V/∂a = 2*D_e*u*(∂u/∂a) = 2*D_e*(r - r_e)*(1 - exp(-a*(r - r_e)))*exp(-a*(r - r_e))
    
    Parameters:
    -----------
    params : array, shape (3,)
        Parameters [D_e, r_e, a]
    r_data : array, shape (n,)
        Distance measurements
    
    Returns:
    --------
    J : array, shape (n, 3)
        Jacobian matrix
    """
    D_e, r_e, a = params
    n = len(r_data)
    J = np.zeros((n, 3))
    
    # Compute common terms (efficient implementation)
    exp_term = np.exp(-a * (r_data - r_e))
    one_minus_exp = 1 - exp_term
    
    # Partial derivatives
    J[:, 0] = one_minus_exp**2  # ∂V/∂D_e
    J[:, 1] = -2 * D_e * a * one_minus_exp * exp_term  # ∂V/∂r_e # missing negative
    J[:, 2] = 2 * D_e * (r_data - r_e) * one_minus_exp * exp_term  # ∂V/∂a
    
    return J


# ============================================================================
# PART 3: BACKTRACKING LINE SEARCH 
# ============================================================================

def backtracking_line_search(params, direction, r_data, V_data, 
                             alpha=1.0, rho=0.5, c=1e-4, max_iter=50):
    """
    Backtracking line search with Armijo condition.
    
    ALGORITHM:
    ----------
    The Armijo condition ensures sufficient decrease:
        f(θ + α*d) ≤ f(θ) + c*α*∇f(θ)^T*d
    
    where:
    - f(θ) is the objective function
    - d is the search direction
    - α is the step size
    - c is a constant (typically 1e-4)
    - ∇f(θ)^T*d is the directional derivative
    
    For nonlinear least squares:
        f(θ) = (1/2)||r(θ)||^2
        ∇f(θ) = J^T*r
        ∇f(θ)^T*d = r^T*J*d
    
    Parameters:
    -----------
    params : array, shape (3,)
        Current parameters
    direction : array, shape (3,)
        Search direction
    r_data, V_data : arrays
        Data
    alpha : float
        Initial step size (usually 1.0)
    rho : float
        Step size reduction factor (0 < rho < 1, typically 0.5)
    c : float
        Armijo constant (typically 1e-4)
    max_iter : int
        Maximum backtracking iterations
    
    Returns:
    --------
    alpha : float
        Selected step size
    """
    f_current = objective_function(params, r_data, V_data)
    
    # Compute directional derivative: ∇f^T * d = r^T * J * d
    res = residual(params, r_data, V_data)
    J = jacobian(params, r_data)
    grad_dot_dir = np.dot(res, np.dot(J, direction))
    
    # Backtracking loop
    for i in range(max_iter):
        # Compute candidate parameters
        params_new = params + alpha * direction
        
        # Evaluate objective at new point
        f_new = objective_function(params_new, r_data, V_data)
        
        # Check Armijo condition for sufficient decrease
        if f_new <= f_current + c * alpha * grad_dot_dir:
            return alpha  # Accept step size
        
        # Reduce step size and try again
        alpha *= rho
    
    # If we reach here, return the smallest step size tried
    return alpha


# ============================================================================
# Steepest Descent Method
# ============================================================================

def steepest_descent(r_data, V_data, params_init, tol=1e-6,
                     max_iter=100, verbose=True):
    """
    Similar to Gauss Newton
    """
    params = params_init.copy()

    history = {
        "params": [params.copy()],
        "objective": [objective_function(params, r_data, V_data)],
        "gradient_norm": [],
        "step_size": [],
        "condition_number": []
    }

    if verbose:
        print(f"{'Iter':<6} {'Objective':<15} {'||Gradient||':<15} "
              f"{'Step Size':<12} {'Cond(JTJ)':<12}")
        print("-" * 85)

    for k in range(max_iter):
        # Residuals and Jacobian
        res = residual(params, r_data, V_data)
        J = jacobian(params, r_data)

        # Gradient = J^T@r
        gradient = np.dot(J.T, res)
        grad_norm = np.linalg.norm(gradient)

        history["gradient_norm"].append(grad_norm)

        # Check convergence
        if grad_norm < tol:
            if verbose:
                print(f"\n✓ Converged! Gradient norm: {grad_norm:.2e} < {tol:.2e}")
            break

        # Steepest descent direction (negative gradient)
        direction = -gradient

        # maybe don't need to check condition number
        try:
            cond_num = np.linalg.cond(J.T @ J)
        except np.linalg.LinAlgError:
            cond_num = np.inf
        history["condition_number"].append(cond_num)

        alpha = backtracking_line_search(params, direction, r_data, V_data)
        history["step_size"].append(alpha)

        params = params + alpha * direction

        history["params"].append(params.copy())
        f_val = objective_function(params, r_data, V_data)
        history["objective"].append(f_val)

        if verbose:
            print(f"{k:<6} {f_val:<15.6e} {grad_norm:<15.6e} "
                  f"{alpha:<12.6f} {cond_num:<12.2e}")
    else:
        if verbose:
            print(f"\nReached max iterations ({max_iter}). Gradient norm: {grad_norm:.2e}")

    return params, history

# ============================================================================
# PART 4: GAUSS-NEWTON METHOD
# ============================================================================

def gauss_newton(r_data, V_data, params_init, tol=1e-6, max_iter=100, verbose=True):
    """
    Gauss-Newton method with backtracking line search.
    
    ALGORITHM:
    ----------
    The Gauss-Newton method solves nonlinear least squares problems:
        min f(θ) = (1/2) * ||r(θ)||^2
    
    where r(θ) is the residual vector.
    
    At each iteration k:
    1. Compute residual r_k and Jacobian J_k
    2. Solve normal equations for search direction:
       (J_k^T J_k) d_k = -J_k^T r_k
    3. Use line search to find step size α_k
    4. Update: θ_{k+1} = θ_k + α_k * d_k
    
    KEY INSIGHTS:
    - Normal equations come from linearizing r(θ) around current point
    - J^T J approximates the Hessian of f (Gauss-Newton approximation)
    - Works well when residuals are small at the solution
    - Cheaper than full Newton's method (no second derivatives needed)
    
    Parameters:
    -----------
    r_data, V_data : arrays
        Measurement data
    params_init : array, shape (3,)
        Initial parameter guess
    tol : float
        Convergence tolerance (on gradient norm)
    max_iter : int
        Maximum iterations
    verbose : bool
        Print iteration info
    
    Returns:
    --------
    params : array, shape (3,)
        Optimized parameters
    history : dict
        Optimization history (objective, gradient norm, step sizes, etc.)
    """
    params = params_init.copy()
    
    # Initialize history tracking
    history = {
        'params': [params.copy()],
        'objective': [objective_function(params, r_data, V_data)],
        'gradient_norm': [],
        'step_size': [],
        'condition_number': []  # Track conditioning of J^T J
    }
    
    if verbose:
        print(f"{'Iter':<6} {'Objective':<15} {'||Gradient||':<15} {'Step Size':<12} "
              f"{'Cond(JTJ)':<12}")
        print("-" * 85)
    
    for k in range(max_iter):
        # Step 1: Compute residual and Jacobian at current parameters
        res = residual(params, r_data, V_data)
        J = jacobian(params, r_data)
        
        # Step 2: Compute gradient: ∇f = J^T * r
        gradient = np.dot(J.T, res)
        grad_norm = np.linalg.norm(gradient)
        
        # Check convergence based on gradient norm
        if grad_norm < tol:
            if verbose:
                print(f"\n✓ Converged! Gradient norm: {grad_norm:.2e} < {tol:.2e}")
            break
        
        # Step 3: Solve normal equations for search direction
        # (J^T J) d = -J^T r
        JTJ = np.dot(J.T, J)
        JTr = np.dot(J.T, res)
        
        # Check condition number (for analysis)
        try:
            cond_num = np.linalg.cond(JTJ)
        except:
            cond_num = np.inf
        
        # Solve with small Tikhonov regularization for numerical stability
        # This helps when J^T J is nearly singular
        regularization = 1e-8 * np.eye(len(params))
        direction = np.linalg.solve(JTJ + regularization, -JTr)
        
        # Step 4: Backtracking line search to find step size
        alpha = backtracking_line_search(params, direction, r_data, V_data)
        
        # Step 5: Update parameters
        params = params + alpha * direction
        
        # Store history
        history['params'].append(params.copy())
        history['objective'].append(objective_function(params, r_data, V_data))
        history['gradient_norm'].append(grad_norm)
        history['step_size'].append(alpha)
        history['condition_number'].append(cond_num)
        
        if verbose:
            print(f"{k:<6} {history['objective'][-1]:<15.6e} {grad_norm:<15.6e} "
                  f"{alpha:<12.6f} {cond_num:<12.2e}")
    else:
        if verbose:
            print(f"\n Maximum iterations ({max_iter}) reached. Gradient norm: {grad_norm:.2e}")
    
    return params, history


# ============================================================================
# PART 5: DATA GENERATION AND ANALYSIS TOOLS
# ============================================================================

def generate_data(params_true, r_range=(0.5, 3.0), n_points=30, noise_level=0.05):
    """Generate synthetic 'experimental' data with Gaussian noise."""
    np.random.seed(42)  # For reproducibility
    r_data = np.linspace(r_range[0], r_range[1], n_points)
    V_true = morse_potential(r_data, params_true)
    
    # Add Gaussian noise
    noise = noise_level * np.random.randn(n_points)
    V_data = V_true + noise
    
    return r_data, V_data


def compute_fit_statistics(r_data, V_data, params_fit):
    """Compute goodness-of-fit statistics."""
    res = residual(params_fit, r_data, V_data)
    n = len(r_data)
    p = 3  # number of parameters
    
    # Sum of squared residuals
    SSR = np.sum(res**2)
    
    # Mean squared error
    MSE = SSR / (n - p)
    
    # Root mean squared error
    RMSE = np.sqrt(MSE)
    
    # R-squared
    V_mean = np.mean(V_data)
    SST = np.sum((V_data - V_mean)**2)
    R_squared = 1 - SSR / SST
    
    # Adjusted R-squared
    R_squared_adj = 1 - (1 - R_squared) * (n - 1) / (n - p - 1)
    
    return {
        'SSR': SSR,
        'MSE': MSE,
        'RMSE': RMSE,
        'R_squared': R_squared,
        'R_squared_adj': R_squared_adj
    }


def plot_results(r_data, V_data, params_true, params_fit, history, plt_title):
    """Plot comprehensive fitting results and convergence history."""
    fig = plt.figure(figsize=(16, 10))
    fig.suptitle(plt_title) # add title to differentiate between Gauss Newton and Steepest Descent
    
    # Plot 1: Data and fitted curve
    ax1 = plt.subplot(2, 4, 1)
    r_fine = np.linspace(r_data.min(), r_data.max(), 200)
    V_true = morse_potential(r_fine, params_true)
    V_fit = morse_potential(r_fine, params_fit)
    V_data_fit = morse_potential(r_data, params_fit)
    
    ax1.plot(r_data, V_data, 'ko', label='Data', markersize=6, alpha=0.7)
    ax1.plot(r_fine, V_true, 'b--', label='True potential', linewidth=2)
    ax1.plot(r_fine, V_fit, 'r-', label='Fitted potential', linewidth=2)
    ax1.set_xlabel('Distance r (Å)', fontsize=11)
    ax1.set_ylabel('Potential Energy (eV)', fontsize=11)
    ax1.set_title('Morse Potential Fit', fontsize=12, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Residuals
    ax2 = plt.subplot(2, 4, 2)
    residuals = V_data - V_data_fit
    ax2.scatter(r_data, residuals, c='red', alpha=0.6, s=50)
    ax2.axhline(y=0, color='k', linestyle='--', linewidth=1)
    ax2.set_xlabel('Distance r (Å)', fontsize=11)
    ax2.set_ylabel('Residual (eV)', fontsize=11)
    ax2.set_title('Residual Plot', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Objective function
    ax3 = plt.subplot(2, 4, 3)
    ax3.semilogy(history['objective'], 'b-o', linewidth=2, markersize=4)
    ax3.set_xlabel('Iteration', fontsize=11)
    ax3.set_ylabel('Objective Function', fontsize=11)
    ax3.set_title('Convergence: Objective', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Gradient norm
    ax4 = plt.subplot(2, 4, 4)
    ax4.semilogy(history['gradient_norm'], 'r-o', linewidth=2, markersize=4)
    ax4.set_xlabel('Iteration', fontsize=11)
    ax4.set_ylabel('||Gradient||', fontsize=11)
    ax4.set_title('Convergence: Gradient Norm', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: Step sizes
    ax5 = plt.subplot(2, 4, 5)
    ax5.plot(history['step_size'], 'g-o', linewidth=2, markersize=4)
    ax5.set_xlabel('Iteration', fontsize=11)
    ax5.set_ylabel('Step Size α', fontsize=11)
    ax5.set_title('Line Search Step Sizes', fontsize=12, fontweight='bold')
    ax5.grid(True, alpha=0.3)
    
    # Plots 6-8: Parameter evolution
    param_names = ['D_e (eV)', 'r_e (Å)', 'a (Å⁻¹)']
    params_history = np.array(history['params'])
    
    for i in range(3):
        ax = plt.subplot(2, 4, 6 + i)
        ax.plot(params_history[:, i], 'purple', linewidth=2, marker='o', 
                markersize=4, label='Estimate')
        ax.axhline(y=params_true[i], color='b', linestyle='--', 
                   linewidth=2, label='True value')
        ax.axhline(y=params_fit[i], color='r', linestyle=':', 
                   linewidth=2, alpha=0.7, label='Final')
        ax.set_xlabel('Iteration', fontsize=11)
        ax.set_ylabel(param_names[i], fontsize=11)
        ax.set_title(f'Parameter: {param_names[i]}', fontsize=12, fontweight='bold')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('gauss_newton_results.png', dpi=300, bbox_inches='tight')
    plt.show()


def analyze_convergence(history):
    """Analyze convergence rate."""
    print("\n" + "="*80)
    print("CONVERGENCE ANALYSIS")
    print("="*80)
    
    obj = np.array(history['objective'])
    grad_norm = np.array(history['gradient_norm'])
    
    # Estimate convergence rate for objective
    if len(obj) > 3:
        errors = obj[1:] - obj[-1]  # Distance from final value
        errors = np.maximum(errors, 1e-16)  # Avoid log(0)
        
        # Linear convergence: e_{k+1} ≈ C * e_k → log(e_{k+1}) ≈ log(C) + log(e_k)
        # Fit log(e_{k+1}) vs log(e_k)
        if len(errors) > 2 and np.all(errors[:-1] > 0):
            log_errors = np.log10(errors)
            valid_idx = np.isfinite(log_errors[:-1]) & np.isfinite(log_errors[1:])
            if np.sum(valid_idx) > 2:
                coeffs = np.polyfit(log_errors[:-1][valid_idx], 
                                   log_errors[1:][valid_idx], 1)
                rate = coeffs[0]
                
                print(f"\nObjective convergence rate estimate: {rate:.3f}")
                if rate < 0.5:
                    print("  → Superlinear convergence")
                elif rate < 1.0:
                    print("  → Between linear and superlinear")
                else:
                    print("  → Approximately linear convergence")
    
    # Analyze step sizes
    step_sizes = np.array(history['step_size'])
    full_steps = np.sum(step_sizes == 1.0)
    print(f"\nStep size statistics:")
    print(f"  Full steps (α=1): {full_steps}/{len(step_sizes)} "
          f"({100*full_steps/len(step_sizes):.1f}%)")
    print(f"  Mean step size: {np.mean(step_sizes):.4f}")
    print(f"  Min step size: {np.min(step_sizes):.4f}")


# ============================================================================
# VERIFICATION FUNCTIONS
# ============================================================================

def verify_jacobian(params, r_data, eps=1e-7):
    """Verify Jacobian computation using finite differences."""
    print("\n" + "="*80)
    print("JACOBIAN VERIFICATION")
    print("="*80)
    
    J_analytical = jacobian(params, r_data)
    J_numerical = np.zeros_like(J_analytical)
    
    # Compute numerical Jacobian using central differences
    for j in range(len(params)):
        params_plus = params.copy()
        params_plus[j] += eps
        params_minus = params.copy()
        params_minus[j] -= eps
        
        # Finite difference approximation
        res_plus = residual(params_plus, r_data, np.zeros_like(r_data))
        res_minus = residual(params_minus, r_data, np.zeros_like(r_data))
        J_numerical[:, j] = (res_plus - res_minus) / (2 * eps)
    
    # Compare
    abs_error = np.linalg.norm(J_analytical - J_numerical)
    rel_error = abs_error / np.linalg.norm(J_numerical)
    
    print(f"\nAnalytical Jacobian vs Numerical (finite differences):")
    print(f"  Absolute error: {abs_error:.2e}")
    print(f"  Relative error: {rel_error:.2e}")
    
    if rel_error < 1e-5:
        print("  ✓ Jacobian implementation is CORRECT!")
    else:
        print("  ✗ Jacobian implementation may have ERRORS.")
    
    # Show sample comparison
    print(f"\nSample comparison (first 3 rows):")
    print(f"{'':20} {'∂V/∂D_e':>15} {'∂V/∂r_e':>15} {'∂V/∂a':>15}")
    print("-" * 70)
    for i in range(min(3, len(r_data))):
        print(f"Analytical (row {i}): {J_analytical[i,0]:15.6e} "
              f"{J_analytical[i,1]:15.6e} {J_analytical[i,2]:15.6e}")
        print(f"Numerical  (row {i}): {J_numerical[i,0]:15.6e} "
              f"{J_numerical[i,1]:15.6e} {J_numerical[i,2]:15.6e}")
        print()
    
    return rel_error

# brute force check if the problem has a global minimum within some subset of the parameter space
def brute_force_search(r_data, V_data,
                       D_bounds=(0.0, 10.0),
                       r_bounds=(0.0, 5.0),
                       a_bounds=(0.0, 5.0),
                       n_grid=80):
    grid_D = np.linspace(*D_bounds, n_grid)
    grid_r = np.linspace(*r_bounds, n_grid)
    grid_a = np.linspace(*a_bounds, n_grid)

    best_val = np.inf
    best_param = None
    for D in grid_D:
        for r_e in grid_r:
            for a in grid_a:
                params = np.array([D, r_e, a])
                val = objective_function(params, r_data, V_data)
                if val < best_val:
                    best_val = val
                    best_param = params.copy()

    # Identify points near the minimum
    near = []
    eps = 1e-4  # tolerance
    for D in grid_D:
        for r_e in grid_r:
            for a in grid_a:
                params = np.array([D, r_e, a])
                val = objective_function(params, r_data, V_data)
                if abs(val - best_val) < eps:
                    near.append(params)
    return best_param, best_val, near

def is_local_minima(r_data, V_data, params, step,
                    tol_obj=1e-8):
    """
    Check whether a parameter vector is a discrete local minima w.r.t
    26 neighboring points in a cubic space

    Params
    ----------
    params : array_like, shape (3,); [D_e, r_e, a]
    step : float or array_like shape (3,)
    tol_obj : float, allowable tolerance for minima

    Returns
    -------
    bool
        True if no 26-neighbor point has objective value < center_val - tol_obj.
    """
    params = np.asarray(params, dtype=float)
    step = np.asarray(step, dtype=float)
    if step.shape == ():
        step = np.full_like(params, step)

    neighbor_offsets = []
    for s0 in (-step[0], 0., step[0]):
        for s1 in (-step[1], 0., step[1]):
            for s2 in (-step[2], 0., step[2]):
                if s0 == 0.0 and s1 == 0.0 and s2 == 0.0:
                    # skip middle point
                    continue
                neighbor_offsets.append(np.array([s0, s1, s2]))

    center_val = objective_function(params, r_data, V_data)

    for delta in neighbor_offsets:
        neighbor = params + delta

        try:
            neighbor_val = objective_function(neighbor, r_data, V_data)
            if neighbor_val < center_val - tol_obj:
                return False
        except Exception as e:
            # just in case the neighboring values breaks the objective function somehow
            print('Error calculating objective_function: ', e)
            continue

    return True

# ============================================================================
# MAIN EXECUTION AND EXPERIMENTS
# ============================================================================

def main():
    """Main execution with comprehensive analysis."""
    print("="*80)
    print("COMPLETE SOLUTION: GAUSS-NEWTON METHOD FOR MORSE POTENTIAL FITTING")
    print("="*80)
    
    # True parameters (H2-like molecule)
    params_true = np.array([4.75, 0.74, 1.44])  # [D_e, r_e, a]
    print(f"\nTrue parameters: D_e={params_true[0]:.3f} eV, "
          f"r_e={params_true[1]:.3f} Å, a={params_true[2]:.3f} Å⁻¹")
    
    # Generate synthetic data
    r_data, V_data = generate_data(params_true, n_points=30, noise_level=0.05)
    print(f"Generated {len(r_data)} data points with 5% noise\n")
    
    # Verify Jacobian implementation
    params_test = np.array([4.0, 0.8, 1.2])
    verify_jacobian(params_test, r_data[:10])
    
    # Initial guess
    params_init = np.array([3.0, 1.0, 1.0])
    print(f"\n{'='*80}")
    print("OPTIMIZATION")
    print("="*80)
    print(f"\nInitial guess: D_e={params_init[0]:.3f}, "
          f"r_e={params_init[1]:.3f}, a={params_init[2]:.3f}\n")
    
    # Run Gauss-Newton optimization
    params_fit, history = gauss_newton(r_data, V_data, params_init, 
                                       tol=1e-6, max_iter=50, verbose=True)
    
    # Print results
    print("\n" + "="*80)
    print("RESULTS")
    print("="*80)
    
    print(f"\n{'Parameter':<15} {'True':<12} {'Initial':<12} {'Fitted':<12} "
          f"{'Abs Error':<12} {'Rel Error (%)':<15}")
    print("-"*80)
    param_names = ['D_e', 'r_e', 'a']
    for i, name in enumerate(param_names):
        abs_error = abs(params_fit[i] - params_true[i])
        rel_error = 100 * abs_error / abs(params_true[i])
        print(f"{name:<15} {params_true[i]:<12.6f} {params_init[i]:<12.6f} "
              f"{params_fit[i]:<12.6f} {abs_error:<12.6f} {rel_error:<15.4f}")
    
    # Fit statistics
    stats = compute_fit_statistics(r_data, V_data, params_fit)
    print(f"\nGoodness of fit:")
    print(f"  RMSE: {stats['RMSE']:.6f} eV")
    print(f"  R²: {stats['R_squared']:.6f}")
    print(f"  Adjusted R²: {stats['R_squared_adj']:.6f}")
    print(f"  Final objective: {history['objective'][-1]:.6e}")
    print(f"  Iterations: {len(history['objective']) - 1}")
    
    # Convergence analysis
    analyze_convergence(history)
    
    # Plot results
    plot_results(r_data, V_data, params_true, params_fit, history, "Gauss Newton")

    # == INSERTED CODE ==

    print("\n================================")

    print("\nparams_true is a local minima: ", is_local_minima(r_data, V_data, params_true, [1e-2, 1e-2, 1e-2], 1e-2))

    print("\n================================\n")

    print("\n================================\nSteepest Descent\n================================\n")

    # run steepest_descent method for comparison
    params_fit_sd, history_sd = steepest_descent(r_data, V_data, params_init, 
                                       tol=1e-6, max_iter=50, verbose=True)
    
    plot_results(r_data, V_data, params_true, params_fit_sd, history_sd, "Steepest Descent")

    print("\n================================\nHessian\n================================\n")

    J = jacobian(params_true, r_data)
    Hessian = J.T@J
    Hessian_inv = np.linalg.inv(Hessian)

    print("\nCondition number of Hessian Matrix : ", np.linalg.norm(Hessian, ord=np.inf) * np.linalg.norm(Hessian_inv, ord=np.inf))
    
    # Check if the problem has a global minimum
    best_param, best_val, near = brute_force_search(r_data, V_data)
    print("\nBest param: ", best_param)
    print("Best val: ", best_val)
    print("near points: ", near)

    # ===================
    
    return params_fit, history


if __name__ == "__main__":
    q1_and_3_main()
    params_fit, history = main()
