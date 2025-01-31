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


def forward_euler_CIR(x, u0, flux_derivative, T, CFL):

    dx = x[1]-x[0]
    u = u0.copy()

    t = 0
    while t < T:
        if np.max(np.abs(flux_derivative(u))) < 1e-8:
            dt = CFL * dx / 2.34
        else:
            dt = CFL * dx / np.max(np.abs(flux_derivative(u)))
        if T - t < dt:
            dt = T - t

        u_padded = np.pad(u, 1, mode="wrap")

        wave_speeds = flux_derivative((u_padded[1:] + u_padded[:-1]) / 2)
        numerical_flux = lambda wave_speed, u_left, u_right: (
            np.maximum(0, wave_speed) * u_left + np.minimum(0, wave_speed) * u_right
        )

        flux_left = numerical_flux(wave_speeds[:-1], u_padded[:-2], u_padded[1:-1])
        flux_right = numerical_flux(wave_speeds[1:], u_padded[1:-1], u_padded[2:])

        u -= dt * (flux_right - flux_left) / dx

        t += dt
    return u

x = np.linspace(0,1,num=100)
# u0 = np.where( (x>.3) & (x<.7), 1.0, 0.0)
U = forward_euler_CIR(x, u0_tp(x), dburgers, .02, .2)
plt.plot(x, u0_tp(x))
plt.plot(x,U)

# %%

def forward_euler_Harten_Hyman(x, u0, flux_derivative, T, CFL):

    dx = x[1]-x[0]
    u = u0.copy()

    t = 0
    while t < T:
        if np.max(np.abs(flux_derivative(u))) < 1e-8:
            dt = CFL * dx / 2.34
        else:
            dt = CFL * dx / np.max(np.abs(flux_derivative(u)))
        if T - t < dt:
            dt = T - t

        u_padded = np.pad(u, 1, mode="wrap")

        u_half = (u_padded[1:] + u_padded[:-1]) /2
        wave_speeds = flux_derivative(u_half)
        sigmas = np.maximum(0,np.maximum(flux_derivative(u_padded[1:])-flux_derivative(u_half),flux_derivative(u_half)-flux_derivative(u_padded[:-1])))
        numerical_flux = lambda wave_speed, sigma, u_left, u_right: (
            np.maximum(0, np.maximum(sigma, wave_speed)) * u_left + np.minimum(0, np.minimum(sigma,wave_speed)) * u_right
        )

        flux_left = numerical_flux(sigmas[:-1], wave_speeds[:-1], u_padded[:-2], u_padded[1:-1])
        flux_right = numerical_flux(sigmas[1:], wave_speeds[1:], u_padded[1:-1], u_padded[2:])

        u -= dt * (flux_right - flux_left) / dx

        t += dt
    return u

x = np.linspace(-1,1,num=100)
t=20

u0 = lambda x: np.where(x<0, 1.0, 0.)
dflux = dburgers

plt.plot(x, u0(x), label="u0")
U = forward_euler_Harten_Hyman(x, u0(x), dflux, t, .2)
plt.plot(x,U, label="Harten&Hyman")
U = forward_euler_CIR(x, u0(x), dflux, t, .2)
plt.plot(x,U, label="CIR")

# plt.plot(x,solution_td(t, x), '--', label="solution")
plt.legend()

# %%


def HLL(x, u0, T, flux, dflux, CFL):

    dx = x[1]-x[0] # assuming even space
    


