# %%

import numpy as np
import matplotlib.pyplot as plt


def euler_flux(U, n: np.array, gamma=1.4):
    assert (
        U.shape[0] == 5
    ), f"Input P must have shape (5) or (5, N) but has shape: {U.shape}"
    rho, u, v, w, p, E = conserved_to_primitive(U, gamma)
    q = u * n[0] + v * n[1] + w * n[2]

    return np.vstack(
        (
            rho * q,
            rho * u * q + p * n[0],
            rho * v * q + p * n[1],
            rho * w * q + p * n[2],
            (E + p) * q,
        )
    )


def conserved_to_primitive(U, gamma=1.4):
    assert (
        U.shape[0] == 5
    ), f"Input P must have shape (5) or (5, N) but has shape: {U.shape}"
    rho, _u, _v, _w, E = U
    u, v, w = _u / rho, _v / rho, _w / rho
    p = (gamma - 1.0) * (E - 0.5 * (u**2 + v**2 + w**2) * rho)

    return np.vstack((rho, u, v, w, p, E))


def roe_averaged_speeds(Ul, Ur, nx, gamma=1.4):
    rhoL, uL, vL, wL, pL, EL = conserved_to_primitive(Ul, gamma)
    rhoR, uR, vR, wR, pR, ER = conserved_to_primitive(Ur, gamma)

    # enthalpies
    HL = (EL + pL) / rhoL
    HR = (ER + pR) / rhoR

    sqrtRhoL = np.sqrt(rhoL)
    sqrtRhoR = np.sqrt(rhoR)
    denom = sqrtRhoL + sqrtRhoR

    uT = (sqrtRhoL * uL + sqrtRhoR * uR) / denom
    vT = (sqrtRhoL * vL + sqrtRhoR * vR) / denom
    wT = (sqrtRhoL * wL + sqrtRhoR * wR) / denom
    HT = (sqrtRhoL * HL + sqrtRhoR * HR) / denom

    unL = uL * nx[0] + vL * nx[1] + wL * nx[2]
    unR = uR * nx[0] + vR * nx[1] + wR * nx[2]
    unT = uT * nx[0] + vT * nx[1] + wT * nx[2]

    vel2T = uT * uT + vT * vT + wT * wT
    cT = np.sqrt((gamma - 1.0) * (HT - 0.5 * vel2T))

    SL = np.minimum(
        unL - np.sqrt((gamma - 1.0) * (HL - 0.5 * (uL * uL + vL * vL + wL * wL))),
        unT - cT,
    )
    SR = np.maximum(
        unR + np.sqrt((gamma - 1.0) * (HR - 0.5 * (uR * uR + vR * vR + wR * wR))),
        unT + cT,
    )
    return SL, SR


def HLL(Ul, Ur, n, gamma=1.4):

    SL, SR = roe_averaged_speeds(Ul, Ur, n, gamma)

    Fl = euler_flux(Ul, n, gamma)
    Fr = euler_flux(Ur, n, gamma)

    t1 = (np.minimum(SR, 0) - np.minimum(SL, 0)) / (SR - SL)
    t2 = 1 - t1
    t3 = (SR * np.abs(SL) - SL * np.abs(SR)) / (2 * (SR - SL))

    return t1 * Fr + t2 * Fl - t3 * (Ur - Ul), SL, SR


def evolve1d(dx, U0, T, gamma=1.4, CFL=0.5):

    U = U0.copy()
    t = 0
    while t < T:

        F, SL, SR = HLL(U[:, :-1], U[:, 1:], np.array([1, 0, 0]), gamma)

        max_speed = np.max(np.abs(np.hstack((SL, SR))))
        dtCFL = CFL * dx / max_speed if max_speed > 1e-14 else 1e-6
        dt = min(T - t, dtCFL)

        U[:, 1:-1] -= (dt / dx) * (F[:, 1:] - F[:, :-1])
        t += dt
    return U


CFL = 0.5
T = 0.4
nx = 1000
x = np.linspace(0, 1, nx)
xh = 0.5 * (x[1:] + x[:-1])
dx = xh[1] - xh[0]

n = np.array([1, 0, 0])

# U0 pour Sod's shock tube
U = np.zeros((5, nx - 1))
U[0] = np.where(xh < 0.5, 1, 0.125)
U[4] = np.where(xh < 0.5, 1.0, 0.1)


U = evolve1d(dx, U, T, gamma=1.4, CFL=CFL)
plt.plot(xh, conserved_to_primitive(U)[0])
plt.title("Euler")
