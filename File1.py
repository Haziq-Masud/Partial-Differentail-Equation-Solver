import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt

# Definning Variables and their values
P = 1.07e5  # Pa            Pressure
R = 8.3145  # m3 Pa/ mol K  Universal gas constant
T = 295.25  # K             Temperature
C = P / (R * T)
# epsilon = 0.4
# epsilon_p = 0.35
Dm = 1.6e-5
Rp = 0.0004
torque = 3
L = 0.064
ax_length = 0.0282
Dl = 1.2 * (10 ** -4)
y_CO2 = 1
# qi_star = 5.501
qi = 9.17
z_100 = 0.064
z_98_5 = (0.064 / 100) * 98.5
z_50 = 0.032
z_2_5 = (0.064 / 100) * 2.5
A = 0.0069191
# v_interstitial = 0.40
v_superficial = 0.16

m = GEKKO()  # Defining Model for the equation
m.time = np.linespace(0, 7200, 72000)
t = m.Param(value=m.time)
P = m.Var(value=1.07e5)
T = m.Var(value=295.15)
R = m.Const(value=8.3145)
v_interstitial = m.Const(value=0.40)
epsilon = m.Const(value=0.4)
epsilon_p = m.Const(value=0.35)
Dp = m.Const(value=(Dm / 3))
rp = m.Const(value=0.0004)
qi_star = m.Const(value=5.501)
Ci = m.Var(value=0)
# Defining Value of ki leaving the concentration part
constant_k1 = m.Const(value=(15 * epsilon_p * Dm) / (3 * (qi_star * rp ** 2)))
# qi_star_qi = [m.Var(value=qi_star-qi) for ]
m.Equation((1 / P * P.dt()) - (1 / T * T.dt()) == (-R * T / P) * ((1 - epsilon) / epsilon))
# 1/P x dP/dt = 1/P.dt()
# 1/T x dP/dt = 1/T.dt()
