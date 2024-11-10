# %%
import numpy as np
import matplotlib.pyplot as plt

def poisson1D(f, cond1, cond2,nint=20):
    a,alpha = cond1
    b, beta = cond2
    h = (b-a)/nint
    x = np.linspace(a,b, nint+1)
    print(x.shape)
    print(x.size)
    Ah = (np.eye(nint+2) * -2 + np.eye(nint+2,k=-1) + np.eye(nint+2,k=1))
    Ah[0,0] = -3/2
    Ah[0,1] = 2
    Ah[0,2] = -1/2
    Ah[-1,-1] = 1
    Ah[-1,-2] = 0
    Ah *= 1/h**2
    F=np.zeros(nint+2)
    F = f(x)
    F[0] = 0
    F[-1] = beta
    Bc = np.zeros(nint+2)
    Bc[0] = alpha/h
    Bc[-1] = beta/h**2
    print(Ah,)
    
    Uh = np.linalg.solve(Ah, F + Bc)
    return x,Uh
    


x, U = poisson1D(np.exp, [0,-5], [3,3], nint=20)
print(x,U)

a,alpha = [0,-5]
b, beta = [3,3]
sol = lambda x: np.exp(x) + (x - b) * (alpha-np.exp(a)) + beta - np.exp(b)

plt.plot(x, U)
plt.plot(x, sol(x))
