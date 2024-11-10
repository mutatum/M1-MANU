# %%
import numpy as np
import matplotlib.pyplot as plt
import math

u = lambda x,y: np.cos(x**2 + y**2)

delta_u = lambda x,y: -(2*np.sin(x**2 + y**2) + 4*x**2 * np.cos(x**2 + y**2))-(2*np.sin(x**2 + y**2) + 4*y**2 * np.cos(x**2 + y**2))

def Lh(u,x,y,h):
    return (u(x+h,y) + u(x, y+h) - 4 * u(x,y) + u(x-h,y) + u(x,y-h))/h**2

def tau(L, Lh, u, x, y, h):
    return np.abs(L(x,y) - Lh(u,x,y,h))

def CoefDF(k, xbar, x):
    x = np.array(x)
    n = len(x)
    A = np.zeros((n, n))
    B = np.zeros((n, 1))
    h = min(x[1:n] - x[0:n-1])
    h2 = min(abs(x - xbar))
    if h2 > 0:
        h = min(h, h2)
    p = n - k
    for i in range(n):
        for j in range(n):
            A[i, j] = (x[j] - xbar) ** i / math.factorial(i)
    B[k] = 1
    coef = np.linalg.solve(A, B)
    coef = coef*h**k
    
    return coef

x0, y0 = (1,1)
H = np.array([5e-1,1e-1, 5e-2, 1e-2, 5e-3])
h = H[0]

plt.loglog(H, [tau(delta_u, Lh, u,x0, y0, h) for h in H], label=r"$\tau_h$")
plt.loglog(H, H**2, label="h**2")


# points = np.array([2, 1, 0, -1, -2])
# coefs = CoefDF(2, x0, points)
# print("Coefs",coefs)
def Lh_order3(u, x, y, h):
    global x0,y0
    discr = np.array([2*h, h, 0, -h, -2*h])
    points = np.ones(5) * x + discr
    # print("\n\npoints: ",points, x, discr,end='\n\n')
    coefs = CoefDF(2, x0, points)
    print("Coefs in function",np.array([u(x,y0) for x in points]), coefs)
    print("Dot: ", np.dot(np.array([u(x0,x) for x in points]), coefs))
    # print('here',np.dot(np.array([u(point,point) for point in points]), coefs))
    return (np.dot(np.array([u(x0,x) for x in points]), coefs) + np.dot(np.array([u(x,y0) for x in points]), coefs))/h**2
print(tau(delta_u, Lh_order3, u, x0, y0, h))
print(tau(delta_u, Lh, u, x0, y0, h))
plt.loglog(H, [tau(delta_u, Lh_order3, u,x0, y0, h) for h in H], label=r"$\tau_h CoefDF$")
plt.loglog(H, H**4, label="h**2")


plt.legend()