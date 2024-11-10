# %%
import numpy as np
from scipy.linalg import hilbert
import matplotlib.pyplot as plt

# Define the function (e.g., quadratic)

def scalar_product(p, q, NX= 100):
    S = np.linspace(-1,1, NX)
    return .5 * np.sum([p(t)*q(t)*2/NX for t in S])

def proj(Q):
    QBase= [lambda x: 1, lambda x: np.sqrt(3) * x, lambda x: 3*np.sqrt(5)/2 * x**2 + np.sqrt(5)/2]
    return lambda x: np.sum([scalar_product(Q, Qi) * Qi(x) for Qi in QBase])

def func(Q, b):
    return scalar_product(Q, proj(Q))

def fd(f, x, h):
    return (f(x+h)- f(x-h)) / (2*h)

# Define constants
rho = 0.1  # Initial step size (learning rate)
tol = 1e-6  # Tolerance for stopping criterion
max_iter = 500  # Maximum number of iterations

# Initialization
x = np.ones(ndim) * 3
Q = lambda x: x**3
print(x)

def gd(x, rho, tol , max_iter):
    hist_func = []  # Stores function values
    hist_grad_norm = []  # Stores norms of gradients
    x_history = []  # Store x values
    successive_grad_scalar = []
    iter_count = 0
    hist_func.append(func(x))
    hist_grad_norm.append(np.linalg.norm(funcp(x)))
    x_history.append(x.copy())
    grad0 = funcp(x)
    # Loop over descent iterations
    while np.linalg.norm(funcp(x))/np.linalg.norm(grad0) > tol and iter_count < max_iter:
        # Calculate new iterate xn+1 using gradient descent step


# Plot convergence results
plt.figure(figsize=(12, 6))

# Plot function value convergence
xfinal, x_history,hist_func, hist_grad_norm, S = gd(x, optimal_step, 1e-3, 5000)
plt.subplot(2, 2, 1)
# plt.plot(hist_func / hist_func[0], label="Functional")
# plt.title("Convergence of Functional")
# plt.xlabel("Iterations")
# plt.ylabel("J(x) / J(x0)")
# plt.legend()
plt.scatter(np.arange(1, 21), x_history[-1])
# print(x_history[-1])
plt.axhline(1, color="r")
plt.title("Convergence of Gradient Norm")
plt.xlabel("Iterations")
plt.ylabel("||∇J(x)|| / ||∇J(x0)||")
plt.legend()

ax = plt.subplot(2,2,2)
ax.set_yscale('log')
plt.plot(S)


# Plot gradient norm convergence
xfinal, x_history,hist_func, hist_grad_norm, S = gd(x, lambda x:0.1, 1e-3, 5000)
plt.subplot(2, 2, 3)
# plt.plot(hist_grad_norm / hist_grad_norm[0], label="Gradient Norm")
plt.scatter(np.arange(1, 21), x_history[-1])
# print(x_history[-1])
plt.axhline(1, color="r")
plt.title("Convergence of Gradient Norm")
plt.xlabel("Iterations")
plt.ylabel("||∇J(x)|| / ||∇J(x0)||")
plt.legend()

ax = plt.subplot(2, 2, 4)
# plt.plot(hist_grad_norm / hist_grad_norm[0], label="Gradient Norm")
ax.set_yscale('log')
plt.plot(S)
# print(x_history[-1])

plt.tight_layout()
plt.show()


# %%
