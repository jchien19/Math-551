import numpy as np
from matplotlib import pyplot as plt

'''
Question 1

1. I do think that using floating point numbers to compute without x being a floating point number. Although there are alternating signs and 
the terms get closer to each other, consecutive terms are not both floats, and thus condition 1 for cancellation error is not met.
There could be room for error due to lack of significant digits, but if the polynomial is chosen such that exact arithmetic is accurate, then the result should be accurate.

2. I Since x is now computed in floating point arithmetic, I think now there is a chance for cancellation error to occur.
Due to the alternating signs, consecutive terms in the series will have different signs, and at large values of k, the terms will be close in distance to each other.
Since x is now a floating point number, the conditions for potential cancellation error are met.

3. ** Alternating series bound?

'''

# plot taylor polynomial approximation of log(1+x)

def taylor_approx(x, N):
    res = np.zeros(N)
    for k in range(1, N+1):
        if k - 2 >= 0:
            res[k-1] += res[k - 2] + (-1)**(k+1) * (1+x) **k / k
        else:
            res[0] = (-1)**(k+1) * (1+x)**k / k
    
    return res

if __name__ == "__main__":
    N = 200
    x = 0

    res = taylor_approx(x, N)

    Ns = [_ for _ in range(1, N+1)]


    plt.plot(Ns, res)
    plt.xlabel("Ns")
    plt.ylabel('res')
    plt.title("Taylor Approx")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.show()

