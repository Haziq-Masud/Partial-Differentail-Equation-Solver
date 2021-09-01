from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt

m = GEKKO()
m.time = np.linspace(0,20,100)
k = 10
y = m.Var(value=5.0)
t = m.Param(value=m.time)
m.Equation(k*y.dt()==-t*y)
m.options.IMODE = 4
m.solve(disp=False)

plt.plot(m.time,y.value)
plt.xlabel('time')
plt.ylabel('y')
plt.show()