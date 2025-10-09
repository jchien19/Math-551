import numpy as np
from matplotlib import pyplot as plt

# Question 1

# q1 part 1
def S(x: int, n: int) -> np.array:
    ans = np.zeros(n+1, float)
    factorial = float(1)
    x = float(x)
    ans[0] = 1.0

    for k in range(1, n + 1):
        factorial *= k
        ans[k] = ans[k - 1] + (x**k / factorial)
    return ans

# q1 part 2 => compute S(-20, m)
def q1p2():
    res: np.array = None
    res = S(-20, 200)
    for i in range(200):
        res[i] = np.log10(abs(res[i] - np.exp(-20))/np.exp(-20))
    n = np.arange(0, 201)

    print("calculated e^-20: ", res[200])

    # plot
    plt.scatter(n, res, s=5)
    plt.xlabel("n")
    plt.ylabel('$\log_{10}(|S(-20, n)- exp(-20)|/exp(-20))$')
    plt.title("Partial Sum Approximation Relative Error")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.show()

def q1p5():
    # calculate e^-20 by calculating e^20 and taking reciprocal
    res: np.array = None
    res = S(20, 200)
    for i in range(200):
        res[i] = np.log10(abs(1/res[i] - np.exp(-20))/np.exp(-20))

    print("calculated e^-20: ", 1/res[200])

# ===============================

# Question 2

def d_c(f, x, h):
    return (f(x + h) - f(x - h)) / (2*h)

def g(x):
    return 1/(1 + x**2)

def g_deriv(x):
    return (-2 * x) / ((1 + x**2)**2)

def q2p3():
    hs = [2.0**-h for h in range(1, 16)]
    err = [abs(d_c(g, 2.7, h) - g_deriv(2.7)) for h in hs]

    log_err = np.log10(err)
    log_hs = np.log10(hs)

    print("log error: ", log_err)
    print("log hs: ", log_hs)
    
    plt.scatter(log_hs, log_err, s=5)
    plt.xlabel("$\log_{10}(h)$")
    plt.ylabel('$\log_{10}(error)$')
    plt.title("Error of Centered Difference and Derivative at x=2.7")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.show()

def q2p4():
    h1s = [2.0**(-h-1) for h in range(0, 15)]
    h2s = [2.0**-h for h in range(0, 15)]
    err1 = [abs(d_c(g, 2.7, h) - g_deriv(2.7)) for h in h1s]
    err2 = [abs(d_c(g, 2.7, h) - g_deriv(2.7)) for h in h2s]

    slopes = (np.log10(err2) - np.log10(err1)) / (np.log10(h2s) - np.log10(h1s))
    plt.scatter(range(0, 15), slopes, s=5)
    plt.xlabel("k")
    plt.ylabel('slope')
    plt.title("Slopes of log Error versus log H at x=2.7")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.show()

def q2p6():
    hs = [2.0**-h for h in range(1, 16)]
    log_hs = np.log10(hs)
    g_deriv_res = g_deriv(1)

    d_cs_1 = [d_c(g, 1, h) for h in hs]
    err_1 = [abs(d_cs_1[i] - g_deriv_res) for i in range(len(hs))]
    
    # plot log_10(err) vs log_10(h) centered at x=1
    log_err_1 = np.log10(err_1)
    plt.scatter(log_hs, log_err_1, s=5)
    plt.xlabel("$\log_{10}(h)$")
    plt.ylabel('$\log_{10}(error)$')
    plt.title("Error of Centered Difference and Derivative at x=1")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.show()

    # # plot slopes
    # h1s = [2.0**(-h-1) for h in range(0, 15)]
    # h2s = [2.0**-h for h in range(0, 15)]
    # err1 = [abs(d_c(g, 1, h) - g_deriv(1)) for h in h1s]
    # err2 = [abs(d_c(g, 1, h) - g_deriv(1)) for h in h2s]

    # slopes = (np.log10(err2) - np.log10(err1)) / (np.log10(h2s) - np.log10(h1s))
    # plt.scatter(range(0, 15), slopes, s=5)
    # plt.xlabel("k")
    # plt.ylabel('slope')
    # plt.title("Slopes for Error of Centered Difference and Derivative at x=1")
    # plt.grid(True, linestyle="--", alpha=0.6)
    # plt.show()


if __name__ == "__main__":
    pass
