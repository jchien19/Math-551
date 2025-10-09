import numpy as np

def q2_2():
    '''
    Derived expression calculation for root
    '''
    r = 6/(9.0**12 + np.sqrt((9.0**12)**2) + 12)
    print("r value, method 1 (alternative): ", r)
    return r

def q2_3():
    '''
    Std quadratic formula calculation for root
    '''
    r = (-1*9.0**12 + np.sqrt((9.0**12)**2 +12 ))/2
    print("r value, method 2 (std form): ", r)
    return r

# 2.4, backwards error is |f(x_a)|
def q2_4(x):
    '''
    Calculate backwards error for question 2 part 4
    '''
    def p(x):
        return x**2 + (9.0**12)*x -3
    
    res = p(x)
    print("Backwards error (q2): ", np.abs(res))

    return res

def p_exp(x):
    '''
    Expanded Wilkinson Polynomial
    '''
    x=float(x)
    return (x**20 - 210*x**19 + 20615*x**18 - 1256850*x**17 + 
        53327946*x**16 - 1672280820*x**15 + 40171771630*x**14 -
        756111184500*x**13 + 11310276995381*x**12 -
        135585182899530*x**11 + 1307535010540395*x**10 -
        10142299865511450*x**9 + 63030812099294896*x**8 -
        311333643161390640*x**7 + 1206647803780373360*x**6 -
        3599979517947607200*x**5 + 8037811822645051776*x**4 -
        12870931245150988800*x**3 + 13803759753640704000*x**2 -
        8752948036761600000*x + 2432902008176640000)

def bisect_once(f,a,b):
    """
    Perform one step of bisection method.
    Inputs
    ------
    * f: continuous function
    * a: real number (float) so that f(a)<0
    * b: real number (float) so that f(b)>0
    Outputs
    -------
    * c, d: floats between a and b so that f(c)<0, f(d)>0, and |c-d| <=␣
    ↪|a-b|/2
    * root: boolean, True if an exact root has been found, False otherwise
    """
    if f(a) >= 0:
        raise Exception("f(a) must be less than zero!")
    elif f(b) <= 0:
        raise Exception("f(b) must be greater than zero!")
    if f((a+b)/2) < 0:
        c = (a+b)/2
        d = b
        root = False
    elif f((a+b)/2) >0:
        c = a
        d = (a+b)/2
        root = False
    else:
        c = (a+b)/2
        d = (a+b)/2
        root = True

    return c,d,root

def wilkinson_bisect(p_func):
    '''
    Bisection algorithm for the Wilkinson algorithm
    a = 15.91
    b = 16.1
    '''
    a = 15.91
    b = 16.1
    root = None
    i = 0

    a, b, root = bisect_once(p_func, a, b)
    while not root and i < 51:
        a, b, root = bisect_once(p_func, a, b)
        # print("iteration: ", i)
        i+=1

    if root:
        print("\nroot found")
    else:
        print("\n50 iterations hit")

    print("Bisection method root: ", a)
    return b

# Wilkinson polynomial 
def wilk_poly(x):
    '''
    For loop implementation of the wilkinson polynomial
    '''
    a = 1
    for i in range(1, 21):
        a *= (x - i)
        # print(i)

    # print("p: ", a)
    return a

if __name__ == "__main__":
    # alternative form of root
    alt_form = q2_2()

    # std quad form calculation
    std_form = q2_3()

    # backwards error for q2
    q2_4(std_form)

    # true forward error for q2
    print("True forward error (q2): ", np.abs(alt_form - std_form))

    # expanded wilkinson polynomial evaluated at 16
    expd_form_16  = p_exp(16)
    print("Expanded wilk poly at 16: ", expd_form_16)
    # print("^ expanded form")

    # expanded wilkinson polynomial root using bisection method
    expd_form_root = wilkinson_bisect(p_exp)

    # prints backwards error using the expanded formula and the bisection algorithm to find a root
    print("For loop implementation, p_exp implementation: ", wilk_poly(expd_form_root), p_exp(expd_form_root))

    # calculate wilkinson polynomial root using root function
    loop_form_root = wilkinson_bisect(wilk_poly)
