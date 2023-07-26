import numpy as np
import matplotlib.pyplot as plt
import mpmath
import random

from dispexac import dispersion


def omegasolve(omega, omegastarT, omegastarN):
    return dispersion(omega, omegastarT, omegastarN, 4, 1, 100)


templist = mpmath.linspace(100, 1000, 100)
LtLnlist = mpmath.linspace(0.2, 0.76, 9)
roots = dict()
for LtLn in LtLnlist:
    rootslist = list()
    for temp in templist:
        root = 0
        while root == 0:
            try:
                if rootslist:
                    x0 = rootslist[-1]
                else:
                    x0 = random.random() * 10j + random.uniform(-1, 1) * 20
                omegafind = lambda omega: omegasolve(omega, temp, LtLn * temp)
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
    roots[LtLn] = rootslist


fig, axes = plt.subplots(3, 3, squeeze=1)
axes = axes.flatten()
for i in range(9):
    ax = axes[i]
    ax.plot(
        templist,
        np.imag(roots[LtLnlist[i]]) / 2,
        label="LtLn" + str(float(LtLnlist[i])),
    )
    ax.legend()
fig.suptitle("kperp = 4")
plt.show()
