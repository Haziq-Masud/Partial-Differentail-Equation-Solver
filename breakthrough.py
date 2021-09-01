from numpy import arange, array, zeros, ones
from matplotlib import pyplot


# Breakthrough curve for CO2 
START_CONCENTRATION = 0.1
STEPSIZE_TIME = 0.001           # s
END_TIME = 0.2                  # s             Total time of experiment
STEPSIZE_Z = 0.0001             # m

#Parameters with values
L = 0.064                       # m             Total length of the column
c = 43.60446577                 # mol/m3        Total gas phase concentration
c_0 = 0.0                       # mol/m3        Inlet concentration of CO2.
D_L = 0.0020112                 # m2/s          Axial dispersion coefficient
ν = 0.40                        # m/s           Interstitial velocity 
ε = 0.4                         # - #           bed voidage
ε_p = 0.35                      # -             Particle void fraction
q_star = 10.83697               # mol/m3        Equilibirum adsoprtion concentration 
q0 = 18.0616                    # mol/m3        Inlet adsorbed concentration
D_m = 1.60e-5                   # m/s^2         Molecular diffusion
D_p = 0.001                     # m             diameter particle
r_p = 0.0005                    # m             radius particle
τ = 3                           # -             turtoisity
μ = 1.4779e-5                   # Pa s          Viscosity
v_s = 0.16                      # m/s           Supervisial velocity

# Gas law: T = c * R / P
P = 1.07e5                      # Pa            Pressure 
R = 8.3145                      # m3 Pa/ mol K  Universal gas constant
T = 295.25                      # K             Temperature 

# -∂P/∂z --> Pressure drop in the column. 
pressure_drop = 150/4 * 1/r_p**2 * ((1 - ε)/ε)**2 * μ * ν   # Pa / m

# ∂c_i/∂t = ∂/∂z (c * D_L * ∂y_i/∂z + c_i * ν) - (1 - ε)/ε * ∂q_i/∂t
# ∂q_i/∂t = k * (q_star - q)
# k = mass transfer coefficient

# k-i / c_i --> k_i is the mass transfer coefficient 
mass_transfer_coefficient = 15 * ε_p * D_m / (r_p ** 2 * τ * q_star)

def gradient(x):
    grad = zeros(len(x))
    grad[:-1] = x[1:] - x[:-1]
    return grad

pressure = arange(P, P - L * pressure_drop, -pressure_drop * STEPSIZE_Z)
temperature = c * R / pressure

z_length = int(L / STEPSIZE_Z)
c_i = ones(z_length) * START_CONCENTRATION #
y_i = zeros(z_length)           # mol fraction
q_i = ones(z_length) * q0
c_gradient = zeros(z_length)

dy = zeros(z_length)

# Plot concentration at different lengths over time in seconds. 
t = 0  # seconds

PLOT_HEIGHT_100  = z_length - 1  # at 100% height
PLOT_HEIGHT_97_5 = int(0.975 * z_length)  # at 97.5% height
PLOT_HEIGHT_50   = int(0.5 * z_length)  # at 50% height
PLOT_HEIGHT_2_5  = int(0.025 * z_length)  # at 2.5% height

y_i100 = []
y_i97_5 = []
y_i50 = []
y_i2_5 = []

while t < END_TIME:
    t += STEPSIZE_TIME
    c_gradient = gradient(c_i) / STEPSIZE_Z                                   # ∂c_i/∂z with boundary condition ∂c_i/∂z(z=L) = 0 
    c_gradient[0] = - v_s / ε * (c_0 - c) / D_L
    y_gradient = gradient(y_i) / STEPSIZE_Z                                   # ∂y_i/∂z 
    y_double_gradient = gradient(y_gradient) / STEPSIZE_Z   # ∂²y_i/∂z²
    
    dq = mass_transfer_coefficient * c_i * (q_star - q_i)                     # ∂q_i/∂t = k_i * (q_star - q_i)
    c_i += (c * D_L * y_double_gradient + ν * c_gradient - (1 - ε)/ε * dq) * STEPSIZE_TIME
    dy = D_L * gradient(dy) - (ν * y_gradient + R * T / P * (1 - ε)/ε * dq * (1 - y_i))
    y_i += dy * STEPSIZE_TIME
    q_i += dq * STEPSIZE_TIME
    plot_variable = c_i
    y_i100.append(plot_variable[PLOT_HEIGHT_100])
    y_i97_5.append(plot_variable[PLOT_HEIGHT_97_5])
    y_i50.append(plot_variable[PLOT_HEIGHT_50])
    y_i2_5.append(plot_variable[PLOT_HEIGHT_2_5])

n = len(y_i100)
time = arange(0, n) * STEPSIZE_TIME * ν / L
time = arange(0, n) * STEPSIZE_TIME

pyplot.figure()
pyplot.plot(time, y_i100, label='100%', color = 'blue', linewidth = 1)
pyplot.plot(time, y_i97_5, label='97.5%', color = 'green', linewidth = 1)
pyplot.plot(time, y_i50, label='50%', color = 'pink', linewidth = 1)
pyplot.plot(time, y_i2_5, label='2.5%', color = 'red', linewidth = 1)
pyplot.ylabel('CO2 Concentration [mol/$m^3$]', fontsize = 10)
pyplot.xlabel('Dimensionless time [-]', fontsize = 10)
pyplot.title('CO2 concentration at different lengths in the column', fontsize = 10)
pyplot.legend(loc=2)
pyplot.show()
