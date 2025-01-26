# %%

import numpy as np
import matplotlib.pyplot as plt


def forward_euler(x, u0, T, CFL):

    dx = x[1]-x[0]
    u = u0.copy()



    flux_derivative = lambda x: np.ones_like(x)
    t = 0
    while t < T:
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
u0 = np.where( (x>.3) & (x<.7), 1.0, 0.0)
U = forward_euler(x, u0, 70, .5)
plt.plot(U)
