# %%
import numpy as np
import matplotlib.pyplot as plt


def fdb(f, x, h):
    return (f(x) - f(x - h)) / h


def fdf(f, x, h):
    return (f(x + h) - f(x)) / h


def fdc(f, x, h):
    return (f(x + h) - f(x - h)) / (2 * h)


def fd3(f, x, h):
    return (2 * f(x + h) + 3 * f(x) - 6 * f(x - h) + f(x - 2 * h)) / (6 * h)


def Lh(u, x, fd, h):
    return fd(u, x, h)


def Lu(x):
    return np.cos(x)


def tau(u, x, fd, h):
    return np.abs(Lu(x) - Lh(u, x, fd, h))


H = [1e-01, 5e-02, 1e-02, 5e-03, 1e-03]
methods = [fdf, fdb, fdc, fd3]
orders = [1, 1, 2, 3]
x = 1
print("h", *[m.__name__ for m in methods], sep="\t")
for h in H:
    print(np.round(-np.log10(h), 2), end="\t")
    for i, fd in enumerate(methods):
        print(np.round((tau(np.sin, x, fd, h)/h**orders[i]), 2), end="\t")
    print()

plt.figure(figsize=(10,6))
[plt.loglog(H, [h**o for h in H],label=str(o)) for o in set(orders)]
[plt.loglog(H, [tau(np.sin, x, fd, h) for h in H],label=fd.__name__) for fd in methods]
plt.axis()
plt.legend()

print()
print("h", *[m.__name__ for m in methods], sep="\t")
for i,h in enumerate(H[:-1]):
    print(np.round(-np.log10(h), 2), end="\t")
    for fd in methods:
        print(f"{np.log(tau(np.sin, x, fd, H[i])/tau(np.sin, x, fd, H[i+1]))/np.log(H[i]/H[i+1]):.10g}", end="\t")
    print()
print(np.round(-np.log10(H[-1]), 2))
# %%
