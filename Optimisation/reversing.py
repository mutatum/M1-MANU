# %%
from scipy import optimize
import numpy as np
V=1000
S=10
#vector function : gradient of Lagrangian
def fun(x):
    global eps
    return [x[1] + x[2] + x[3]*x[2]*x[1] + x[4]*x[1],       #dL/dx1
            x[0] + x[2] + x[3]*x[0]*x[2] + x[4]*x[0],       #dL/dx2
            x[0] + x[1] + x[3]*x[0]*x[1] + 0,               #dL/dx3
            x[0]*x[1]*x[2] - V + eps,                             #dL/dp1 : contrainte C1=0
            x[0]*x[1] - S + eps                                  #dL/dp2 : contrainte C2=0
            ]

def functional(x):
    return x[0]*x[1] + x[1] * x[2] + x[0]*x[2]

# #gradient of vector function  (gradient of gradient of Lagrangian)
def jacobian(x):
    return np.array([[0 , 1+x[3]*x[2]+x[4], 1+x[3]*x[1], x[2]*x[1], x[1]],
                     [1+x[3]*x[2]+x[4], 0, 1+x[3]*x[0], x[2]*x[0], x[0]],
                     [1+x[3]*x[1], 1+x[3]*x[0], 0, x[0]*x[1], 0],
                     [x[1]*x[2] ,x[0]*x[2], x[0]*x[1] , 0, 0],
                     [ x[1] , x[0] ,  0 , 0, 0]                     
                    ])

xinit = np.ones(5)
eps = 0
sol = optimize.root(fun, xinit, jac=jacobian)
x = sol.x
eps = 1e-4
soleps = optimize.root(fun, xinit, jac=jacobian)
xeps = soleps.x

dJ = fun(x)[:3]
dC1 = jacobian(x)[3,:3]
dC2 = jacobian(x)[4,:3]
print(dJ)

dC = np.column_stack((dC1, dC2))
print(dC)
print(np.linalg.lstsq(dC, dJ, rcond=None))

# print(dJdCi(0, x))
print(x)
print(xeps)
print((functional(xeps) - functional(x))/eps)
print(x[0] * x[1] *x[2])
print(x[0] * x[1])