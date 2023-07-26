import numpy as np
import matplotlib.pyplot as plt
import mpmath
import random

from dispexac import dispersion


def omegasolve(omega, omegastarT, omegastarN, kperp):
    return dispersion(omega, omegastarT, omegastarN, kperp, 1, 150)


# def omegamin(omega_array, kperp):
#     return mpmath.fabs(dispersion(omega_array[0]+omega_array[1]*1j, 100, 20, kperp, 1, 6))

# def omegader(omega, kperp):
#     return dispersionder(omega, 100, 30, kperp, 1, 6)

kperplist = mpmath.linspace(1, 6, 150)
LbLtlist = mpmath.linspace(100, 500, 9)
LtLn = 0.5
roots = dict()
for LbLt in LbLtlist:
    rootslist = list()
    for kperp in kperplist:
        root = 0
        while root == 0:
            try:
                if rootslist:
                    x0 = rootslist[-1]
                else:
                    x0 = random.random() * 10j + random.uniform(-1, 1) * 20
                omegafind = lambda omega: omegasolve(omega, LbLt, LtLn * LbLt, kperp)
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
        rootslist.append((complex(root)))
    roots[LbLt] = rootslist


fig, axes = plt.subplots(3, 3, squeeze=1)
axes = axes.flatten()
for i in range(9):
    ax = axes[i]
    ax.plot(
        kperplist,
        np.imag(roots[LbLtlist[i]]) / 2,
        label="LbLt" + str(float(LbLtlist[i])),
    )
    ax.legend()
fig.suptitle("LtLn = 0.5")
plt.show()
