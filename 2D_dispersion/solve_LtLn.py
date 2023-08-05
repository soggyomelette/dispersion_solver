import numpy as np
import matplotlib.pyplot as plt
import mpmath
import random

from dispexac import dispersion

# increase m if dispersion relation doesn't converge
# decrease m if root finding takes too much time


def omegasolve(omega, omegastarT, omegastarN):
    return dispersion(omega, omegastarT, omegastarN, 1, 1, 40)


n = 2

LtLnlist = mpmath.linspace(0.01, 0.8, 80)
LbLtlist = mpmath.linspace(5, 11, n**2)
roots = dict()
for LbLt in LbLtlist:
    rootslist = list()
    for LtLn in LtLnlist:
        root = 0
        while root == 0:
            try:
                if rootslist:
                    x0 = rootslist[-1]
                else:
                    # use a good initial guess or
                    # a random value in the complex plane
                    x0 = random.random() * 10j + random.uniform(-1, 1) * 20
                omegafind = lambda omega: omegasolve(omega, LbLt, LtLn * LbLt)
                root = mpmath.findroot(omegafind, x0, tol=1e-10)
                if rootslist:
                    if (
                        mpmath.fabs(mpmath.fabs(root) - mpmath.fabs(rootslist[-1]))
                        / mpmath.fabs(rootslist[-1])
                        > 0.6
                    ):
                        root = 0
            except Exception as e:
                print(e)
                root = 0
        print(complex(root))
        rootslist.append(complex(root))
    roots[LbLt] = rootslist


fig, axes = plt.subplots(n, n, squeeze=0)
axes = axes.flatten()
for i in range(n**2):
    ax = axes[i]
    ax.plot(
        LtLnlist,
        np.imag(roots[LbLtlist[i]]) / 2,
        label="LbLt" + str(float(LbLtlist[i])),
    )
    ax.legend()
fig.suptitle("kperp = 0.01, tau = 1")
plt.show()
