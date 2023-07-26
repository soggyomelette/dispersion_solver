import mpmath


mpmath.mp.dps = 120


def Z(x):
    return 1j * mpmath.sqrt(mpmath.mp.pi) * mpmath.exp(-x**2)*mpmath.erfc(-1j*x)

calculatederb = {}
def derbIab(zeta, m):
    if m == 1:
        return 0.5 * (
            (1 + zeta) * (Z(mpmath.sqrt(zeta / 2))) ** 2
            + 4 * mpmath.sqrt(zeta / 2) * Z(mpmath.sqrt(zeta / 2))
            + 2
        )
    elif m == 0:
        return -0.5 * (Z(mpmath.sqrt(zeta / 2))) ** 2
    elif m in calculatederb:
        return calculatederb[m]
    else:
        calculatederb[m] = -(
            (2 * m - 1 + zeta) * derbIab(zeta, m - 1)
            + 2 * (m - 1) * zeta * derbIab(zeta, m - 2)
            + ((-1) ** m) * mpmath.fac(m - 1)
        )
        return calculatederb[m]

calculatedera = {}
def deraIaa(zeta, m):
    if m == 0:
        return 1 / 4 * (Z(mpmath.sqrt(zeta / 2))) ** 2 + mpmath.sqrt(zeta / 2) * Z(mpmath.sqrt(zeta / 2)) * (
            1 + mpmath.sqrt(zeta / 2) * Z(mpmath.sqrt(zeta / 2))
        )
    elif m == 1:
        return (
            -(1 + zeta) * deraIaa(zeta, 0)
            - 0.5 * (Z(mpmath.sqrt(zeta / 2))) ** 2
            - 2 * mpmath.sqrt(zeta / 2) * Z(mpmath.sqrt(zeta / 2))
            - zeta * (1 + mpmath.sqrt(zeta / 2) * Z(mpmath.sqrt(zeta / 2)))
            - 1.5
        )
    elif m in calculatedera:
        return calculatedera[m]
    else:
        calculatedera[m] = -(
            derbIab(zeta, m)
            + zeta * derbIab(zeta, m - 1)
            + (2 * m - 1 + zeta) * deraIaa(zeta, m - 1)
            + 2 * (m - 1) * zeta * deraIaa(zeta, m - 2)
            + ((-1)) ** (m + 1) * mpmath.fac(m - 1) * (m - 0.5)
        )
        return calculatedera[m]


'''def derwderbIab(zeta, m):
    if m == 0:
        return Z(mpmath.sqrt(zeta / 2)) * (1 + mpmath.sqrt(zeta) * Z(mpmath.sqrt(zeta / 2))) / mpmath.sqrt(2 * zeta)
    elif m == 1:
        return 0.5 * (Z(mpmath.sqrt(zeta / 2))) ** 2 - 2 * (1 + zeta) / (mpmath.sqrt(2 * zeta)) * (
            1 + mpmath.sqrt(zeta / 2) * Z(mpmath.sqrt(zeta / 2))
        ) + mpmath.sqrt(2 / zeta) * Z(mpmath.sqrt(zeta / 2)) - 2 * (
            1 + mpmath.sqrt(zeta / 2) * Z(mpmath.sqrt(zeta / 2))
        )
    else:
        return -(
            (2 * m - 1 + zeta) * derwderbIab(zeta, m - 1)
            + derbIab(zeta, m - 1)
            + 2 * (m - 1) * derbIab(zeta, m - 2)
            + 2 * (m - 1) * zeta * derwderbIab(zeta, m - 2)
        )


def derwderaIaa(zeta, m):
    if m == 0:
        return -Z(mpmath.sqrt(zeta / 2)) * (1 + mpmath.sqrt(zeta / 2) * Z(mpmath.sqrt(zeta / 2))) / 2 / mpmath.sqrt(
            2 * zeta
        ) + Z(mpmath.sqrt(zeta / 2)) / 2 / mpmath.sqrt(2 * zeta) - (
            1 + mpmath.sqrt(zeta / 2) * Z(mpmath.sqrt(zeta / 2))
        ) / 2 + (
            Z(mpmath.sqrt(zeta / 2))
        ) ** 2 / 2 - mpmath.sqrt(
            zeta / 2
        ) * Z(
            mpmath.sqrt(zeta / 2)
        ) * (
            1 + mpmath.sqrt(zeta / 2) * Z(mpmath.sqrt(zeta / 2))
        )
    elif m == 1:
        return (
            -deraIaa(zeta, 1)
            - (1 + zeta) * derwderaIaa(zeta, 0)
            + Z(mpmath.sqrt(zeta / 2))
            * (1 + mpmath.sqrt(zeta / 2) * Z(mpmath.sqrt(zeta / 2)))
            / mpmath.sqrt(2 * zeta)
            - Z(mpmath.sqrt(zeta / 2)) / mpmath.sqrt(2 * zeta)
            - zeta
            * (
                Z(mpmath.sqrt(zeta / 2)) / 2 / mpmath.sqrt(2 * zeta)
                - (1 + mpmath.sqrt(zeta / 2) * Z(mpmath.sqrt(zeta / 2)) / 2)
            )
        )
    else:
        return -(
            derwderbIab(zeta, m)
            + derbIab(zeta, m - 1)
            + zeta * derwderbIab(zeta, m - 1)
            + deraIaa(zeta, m - 1)
            + (2 * m - 1 + zeta) * derwderaIaa(zeta, m - 1)
            + 2 * (m - 1) * deraIaa(zeta, m - 2)
            + 2 * (m - 1) * zeta * derwderaIaa(zeta, m - 2)
        )'''


def dispersion(omega, omegastarT, omegastarN, kperp, tau, m):
    disp = -1 - tau
    for i in range(0, m + 1):
        disp = disp + (
            (-omega + omegastarN - 1.5 * omegastarT) * derbIab(omega, i)
            - omegastarT * deraIaa(omega, i)
        ) * mpmath.fac(2 * i) / ((mpmath.fac(i)) ** 4) * (0.5 * kperp) ** (
            2 * i
        )
    calculatedera.clear()
    calculatederb.clear()
    return disp


'''def dispersionder(omega, omegastarT, omegastarN, kperp, tau, m):
    dispder = 0 
    for i in range(0, m + 1):
        dispder = dispder + (
            -derbIab(omega, i)
            + (-omega + omegastarN - 1.5 * omegastarT) * derwderbIab(omega, i)
            - omegastarT * derwderaIaa(omega, i)
        ) * mpmath.fac(2 * i) / ((mpmath.fac(i)) ** 4) * (0.5 * kperp) ** (
            2 * i
        )
    return dispder'''


'''reom = mpmath.linspace(0, 10, 300)
imom = mpmath.linspace(0, 10, 300)
X, Y = np.meshgrid(reom, imom)
displot = np.absolute(dispersion(X + Y * 1j, 300, 150, 0.8, 1, 50))
fig, ax = plt.subplots()
con =  ax.contour(X, Y, displot)
ax.clabel(con)
plt.show()'''


