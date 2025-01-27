# %%

import numpy as np
import matplotlib.pyplot as plt


def burgers(u):
    return 0.5 * u * u


def dburgers(u):
    return u


def buckley(u):
    return 4 * u**2 / (4 * u**2 + (1 - u) ** 2)


def dbuckley(u):
    return 8 * u * (1 - u) / (5 * u**2 - 2 * u + 1) ** 2


def u0_cours(x):
    return np.where(x < 0.3, 0.0, np.where(x < 0.7, -1.0, 0.5))


def solution_exacte_cours(t, x):
    if t > 1:
        raise ValueError("t must be less than 1")
    return np.where(x < t, 1, np.where(x <= 1, (1 - x) / (1 - t), 0))


def solution_td(t, x):
    t_etoile = 0.8
    if t == 0:
        return u0_cours(x)
    if t < t_etoile:
        return np.where(
            x < 0.3 - t / 2,
            0,
            np.where(
                x < 0.7 - t,
                -1,
                np.where(
                    (-1 <= (x - 0.7) / t) & ((x - 0.7) / t <= 0.5), (x - 0.7) / t, 0.5
                ),
            ),
        )
    return np.where(
        x < -np.sqrt(0.8 * t) + 0.7, 0, np.where(x <= 0.7 + t / 2, (x - 0.7) / t, 0.5)
    )


def u0_buckley(x):
    x = np.array(x)
    return np.where((-0.5 < x) & (x < 0), 1.0, 0.0)


u0_tp  = lambda x: np.sin(np.pi * 2 * x)
du0_tp = lambda x: 2 * np.pi * np.cos(2 * np.pi * x)


def forward_euler_CIR(x, u0, df, T, CFL):
    dx = x[1]-x[0]
    u = u0.copy()
    t = 0
    while t < T:
        if np.max(np.abs(df(u))) < 1e-8:
            dt = CFL * dx / 2.34
        else:
            dt = CFL * dx / np.max(np.abs(df(u)))
        if T - t < dt:
            dt = T - t
        
        uB = np.pad(u, 1, mode='edge')

        uL = uB[:-1]
        uC = (uB[1:] + uB[:-1])/2
        uR = uB[1:]
        
        w = df(uC)
        h = uL * np.maximum(0, w) + uR * np.minimum(0, w)

        u -= dt * (h[1:] - h[:-1])/dx

        t += dt
    return u

x = np.linspace(0,1,num=100)
# u0 = np.where( (x>.3) & (x<.7), 1.0, 0.0)
U = forward_euler_CIR(x, u0_tp(x), dburgers, .09, .4)
plt.plot(x, u0_tp(x))
plt.plot(x,U)

# %%
# def forwar

def forward_euler_Harten_Hyman(x, u0, f, df, T, CFL):
    dx = x[1]-x[0]
    u = u0.copy()
    t = 0
    while t < T:
        if np.max(np.abs(df(u))) < 1e-8:
            dt = CFL * dx / 2.34
        else:
            dt = CFL * dx / np.max(np.abs(df(u)))
        if T - t < dt:
            dt = T - t
        
        uB = np.pad(u, 1, mode='edge')

        uL = uB[:-1]
        uC = (uB[1:] + uB[:-1])/2
        uR = uB[1:]
        
        s = np.maximum(0, np.maximum(df(uR)-df(uC), df(uC)-df(uL)))
        w = df(uC)
        h = uL * np.maximum(0, np.maximum(s, w)) + uR * np.minimum(0, np.minimum(-s, w))

        u -= dt * (h[1:] - h[:-1])/dx

        t += dt
    return u

x = np.linspace(-1,1,num=20)
t=1

u0 = lambda x: np.where(x<0, 1.0, 0.)
flux = burgers
dflux = dburgers

plt.plot(x, u0(x), label="u0")
U = forward_euler_Harten_Hyman(x, u0(x), flux, dflux, t, .4)
plt.plot(x,U, label="Harten&Hyman")
U = forward_euler_CIR(x, u0(x), dflux, t, .4)
plt.plot(x,U, label="CIR")

sol = lambda x,t: np.where(x/t > 1, 0., 1.)
plt.plot(x, sol(x,t), '--', label="solution")
plt.legend()

# %%

x = np.linspace(-1,2,num=100)
t = .2
# u0 = np.where( (x>.3) & (x<.7), 1.0, 0.0)
U = forward_euler_CIR(x, u0_cours(x), dburgers, t, .4)
# plt.plot(x, u0_cours(x))
plt.plot(x,U, label='CIR')
U = forward_euler_Harten_Hyman(x, u0_cours(x), burgers, dburgers, t, .4)
plt.plot(x, U, label='HH')
plt.plot(x,solution_td(t, x),label='solution')
plt.legend()
# %%