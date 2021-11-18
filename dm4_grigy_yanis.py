"""td5_mwe_odeint.py"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#############################
# PARAMÈTRES DE LA RÉSOLUTION
#############################
t0 = 0                   # bornes de l'intervalle de résolution
tf = t0 + 0.03                   # en secondes
dt = 1e-5                # pas de temps en secondes
n  = int((tf-t0)/dt + 1) # nombre de points
t = np.linspace(t0,tf,n)
r = 6.
E = 12.
C = 10**(-6)
Q = 14
Q0 = C*E
I1 = E/r
mu = omega0/(2*Q)
omega = u*np.sqrt(4*(Q**2)-1)
X0 = [0,I1]               # conditions initiales : [q(0), dq(0)/dt]
omega0 = 1610   # pulsation propre en s^-1

###############################################
# FONCTION ASSOCIÉE À L'ÉQUATION DIFFÉRENTIELLE
###############################################
def F(V, t):
	x, y = V
	return [y, -y*(omega0/Q) + (omega0**2)*(Q0 - x)]
def G(t):
    return np.ext(-mu*t)*(((I1 - mu*Q0)/omega) * np.sin(omega*t) - Q0  * np.cos(omega*t) + Q0
########################################
# RÉSOLUTION ET REPRÉSENTATION GRAPHIQUE
########################################
X = odeint(F, X0, t) # résolution
q = X[:,0]           # récupération des données
plt.plot(t, q, label="numerique")
plt.xlabel("Temps (s)")
plt.ylabel("Chqrge (C)")
plt.plot(t, G(t), label="analytique"))
plt.legend()
plt.show()





