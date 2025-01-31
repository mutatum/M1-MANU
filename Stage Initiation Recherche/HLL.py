# %%

import numpy as np
import matplotlib.pyplot as plt


def euler_flux(U, n: np.array, gamma=1.4):
    assert (
        U.shape[0] == 5
    ), f"Input P must have shape (5) or (5, N) but has shape: {U.shape}"
    rho, u, v, w, p, E, e = conserved_to_primitive(U, gamma)
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


def primitive_to_conserved(U, gamma=1.4):
    rho, u, v, w, p = U
    e = p / (gamma - 1)
    E = e + .5 * rho * (u**2+v**2+w**2)
    return np.vstack((rho, rho*u, rho*v, rho*w, E))

def conserved_to_primitive(U, gamma=1.4):
    """
        Takes U = ( density rho,
                    rho * velocity u, 
                    rho * vel v, 
                    rho *vel w, 
                    Energy)
        and gives U_cons= ( density rho,
                            velocity u,
                            velocity v,
                            velocity w,
                            pressure p,
                            total energy per unit volume E,
                            internal energy e,
        )
    """
    assert (
        U.shape[0] == 5
    ), f"Input P must have shape (5) or (5, N) but has shape: {U.shape}"
    rho, _u, _v, _w, E = U
    u, v, w = _u / rho, _v / rho, _w / rho
    e = E - .5 * (u**2 + v**2 + w**2) * rho
    p = (gamma - 1.0) * e

    return np.vstack((rho, u, v, w, p, E, e))


def roe_averaged_speeds(Ul, Ur, nx, gamma=1.4):
    rhoL, uL, vL, wL, pL, EL, eL = conserved_to_primitive(Ul, gamma)
    rhoR, uR, vR, wR, pR, ER, eR = conserved_to_primitive(Ur, gamma)

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

def HLLC(Ul, Ur, n, gamma=1.4):

    SL, SR = roe_averaged_speeds(Ul, Ur, n, gamma)
    rhol, ul, vl, wl, pl, El, el = conserved_to_primitive(Ul)
    rhor, ur, vr, wr, pr, Er, er = conserved_to_primitive(Ur)
    ql = ul * n[0] + vl * n[1] + wl * n[2]
    qr = ur * n[0] + vr * n[1] + wr * n[2]
    SM = (rhor * qr (SR - qr) - rhol*ql*(SL-ql) + pl - pr) / (rhor*(SR-qr) - rhol*(SL-ql))

    Fl = euler_flux(Ul, n, gamma)
    Fr = euler_flux(Ur, n, gamma)

    F = np.where(SR<0, Fr, np.where(SL>0, Fl, (SR * Fl - SL*Fr + SL*SR*(Ur-Ul))/(SR - SL) ))
    return F, SL, SR

def HLL(Ul, Ur, n, gamma=1.4):

    SL, SR = roe_averaged_speeds(Ul, Ur, n, gamma)

    Fl = euler_flux(Ul, n, gamma)
    Fr = euler_flux(Ur, n, gamma)

    F = np.where(SR<0, Fr, np.where(SL>0, Fl, (SR * Fl - SL*Fr + SL*SR*(Ur-Ul))/(SR - SL) ))
    return F, SL, SR


def evolve1d(dx, U0, T, gamma=1.4, CFL=0.5):

    U = U0.copy()
    t = 0
    while t < T:

        F, SL, SR = HLL(U[:, :-1], U[:, 1:], np.array([1, 0, 0]), gamma)

        max_speed = np.max(np.abs(np.hstack((SL, SR))))
        dtCFL = CFL * dx / (2*max_speed) if max_speed > 1e-14 else 1e-6
        dt = min(T - t, dtCFL)

        U[:, 1:-1] -= (dt / dx) * (F[:, 1:] - F[:, :-1])
        U[:, 0] = U[:,1]
        U[:, -1] = U[:,-2]
        t += dt

    return U


CFL = 0.5
T = 0.2
nx = 1000
x = np.linspace(-.2, .5, nx)
xh = 0.5 * (x[1:] + x[:-1])
dx = xh[1] - xh[0]

n = np.array([1, 0, 0])

# U0 pour Toro's shock tube problem
U = np.zeros((5, nx - 1))
U[0] = np.where(xh < 0, 1, 0.125) # density rho
U[1] = np.where(xh < 0, .75, 0.) # velocity u
U[4] = np.where(xh < 0, 1, 0.1) # pressure p
U = primitive_to_conserved(U)
print(U.shape)


U = evolve1d(dx, U, T, gamma=1.4, CFL=CFL)

plt.figure(figsize=(10,10))

ax = plt.subplot(2,2,1)
ax.plot(xh, U[0], "-", label="HLL")
ax.set_xlabel("space")
ax.set_ylabel("density (rho)")
ax.legend()

ax = plt.subplot(2,2,2)
ax.plot(xh, conserved_to_primitive(U)[1], "-", label="HLL")
ax.set_xlabel("space")
ax.set_ylabel("velocity (u)")
ax.legend()

ax = plt.subplot(2,2,3)
ax.plot(xh, conserved_to_primitive(U)[4], "-", label="HLL")
ax.set_xlabel("space")
ax.set_ylabel("Pressure (p)")
ax.legend()

ax = plt.subplot(2,2,4)
ax.plot(xh, conserved_to_primitive(U)[6]/U[0], "-", label="HLL")
ax.set_xlabel("space")
ax.set_ylabel("Energy (E)")
ax.legend()
