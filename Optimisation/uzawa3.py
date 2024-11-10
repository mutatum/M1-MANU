# %%
import numpy as np

# Définition de la fonctionnelle et de la contrainte
J = lambda x: x[0]**2 + 2 * x[1]**2 + 2 * x[0] + 3 * x[1]
C = lambda x: x[0] + x[1] - 1

def uzawa_algorithm(A, B, c, rho, max_iter=1000, tol=1e-6):
    # Initialisation des variables
    n = A.shape[0]   # Taille du vecteur x
    m = B.shape[0]   # Nombre de contraintes
    x = np.zeros(n)  # Initialisation de x
    p = np.array([0.1])  # Initialisation de p (les multiplicateurs de Lagrange)
    
    # Itérations de l'algorithme d'Uzawa
    for k in range(max_iter):
        grad = lambda x,p : np.array([2*x[0] + 2, 4*x[1] + 3]) + B.T@p
        # Mise à jour de x en résolvant A x = b - B^T p
        # print(x,rho*grad(x),"p=", p)
        x -= rho * grad(x,p)
        # print(x,rho)
       
        # Mise à jour de p
        # print(B@x-c, "here",rho, 'p=',p)
        p_next = p + rho * (B @ x - c)
        
        # Condition d'arrêt : norme de la différence entre p et p_next
        if np.linalg.norm(p_next - p) < tol:
            break
        
        # Mise à jour de p
        p = p_next
    
    return x, p

# Exemple d'utilisation
# Matrice A symétrique définie positive
A = np.array([[2, 0], [0, 4]])

# Matrice B et vecteurs b, c
B = np.array([[1, 1]])
c = np.array([1])

# Paramètre de l'algorithme
rho = 0.1


A = np.array([[2, 0], [0, 4]])  # Matrice dérivée des termes quadratiques de J(x)
B = np.array([[1, 1]])          # Matrice dérivée des contraintes
BBt = B @ B.T
eigvals_A = np.linalg.eigvalsh(A)  # Valeurs propres de A
eigvals_BBt = np.linalg.eigvalsh(BBt)  # Valeurs propres de BB^T

lambda_min_A = np.min(eigvals_A)
lambda_max_BBt = np.max(eigvals_BBt)

# Condition sur rho pour assurer la convergence de l'algorithme d'Uzawa
# rho = 2 * lambda_min_A / lambda_max_BBt -1
# Exécution de l'algorithme d'Uzawa
x_opt, p_opt = uzawa_algorithm(A, B, c, rho)
# Affichage des résultats
print("Solution optimale x:", x_opt)
print("Multiplicateurs de Lagrange p:", p_opt)
