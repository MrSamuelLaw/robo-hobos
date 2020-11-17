get_ipython().run_line_magic("reset", " -f")


get_ipython().run_line_magic("matplotlib", " inline")


import sympy as sm
import sympy.physics.mechanics as me
import matplotlib.pyplot as plt
from IPython.display import Latex, display
from mathtools.root_locus import *
from mathtools.symbolic import sym2poly
from numpy import linspace, array, pi, sqrt, polyfit, vstack, sign
from scipy.integrate import RK45


def show(eq_lhs: str, eq_rhs: str):
    eq = f"${eq_lhs} = {eq_rhs}$"
    return display(Latex(eq))


# create the inertial frame
N = me.ReferenceFrame("N")
pN = me.Point("pN")
pN.set_vel(N, 0)


# define the constant g
g, gear_ratio = sm.symbols("g, gear_ratio")


# declare independent variables
q1, u1 = me.dynamicsymbols("q1, u1")
q1d, u1d = me.dynamicsymbols("q1, u1", 1)
q1dd = me.dynamicsymbols("q1", 2)


# declare dependent variables
q0 = -q1*gear_ratio
q0d = -q1d*gear_ratio
q0dd = -q1dd*gear_ratio
u0 = -u1*gear_ratio


# create the bodies
bodyA = me.Body(
    name = "bodyA",
    frame = N.orientnew("A", "axis", (q0,N.z)),
)
bodyB = me.Body(
    name = "bodyB",
    frame = N.orientnew("B", "axis", (q1, N.z))
)


# define the motor shaft diameter
bodyA.D = sm.symbols("D_A")
# define the link length to bodyB COM
bodyB.L = sm.symbols("L_B")


# define positions of the COMs
bodyA.masscenter = pN.locatenew("pA", 0)
bodyB.masscenter = pN.locatenew("pB", bodyB.L*bodyB.frame.x)


# define the velocities
bodyA.masscenter.set_vel(N, bodyA.masscenter.pos_from(pN).dt(N))
bodyB.masscenter.set_vel(N, bodyB.masscenter.pos_from(pN).dt(N))


# define accelerations
bodyA.masscenter.set_acc(N, bodyA.masscenter.vel(N).dt(N))
bodyB.masscenter.set_acc(N, bodyB.masscenter.vel(N).dt(N))


# apply kanes method
KM = me.KanesMethod(
    frame = N,
    q_ind = [q1],
    u_ind = [u1],
    kd_eqs = [q1d - u1]
)


# create a list of the bodies
bodies = [bodyA, bodyB]


# define the applied loads
tau1 = sm.symbols("tau_1")
loads = [
    # point/frame, force.moment
    (bodyA.masscenter, -g*bodyA.mass*N.z),
    (bodyB.masscenter, -g*bodyB.mass*N.z),
    (bodyB.frame, tau1*bodyA.frame.z)
]


# apply the loads and bodies to the KM object
KM.kanes_equations(bodies, loads);


# extract the mass and force matrix
MM = KM.mass_matrix
FM = KM.forcing


# seperate the force out into driving and non-driving
A, b = sm.linear_eq_to_matrix(FM, tau1)
TM = A.inv()*(MM*sm.Matrix([q1dd]) + b)
show("TM", me.vlatex(TM))


# create a dict to hold numeric values
vals = {
    # lengths
    bodyB.L: 0.02162775,         # m
    bodyB.mass:   0.00728203,    # kg 
    "bodyB_izz":  6.39e-06,      # kg*m^2
    "bodyA_izz":  3.3e-07,       # kg*m^2
    gear_ratio: 144,             # a.u
}


# create variables for desired output
q1_desired, q1d_desired = sm.symbols("q1_desired, q1d_desired")

# creatue variables for current output
q1_current, q1d_current = sm.symbols("q1_current, q1d_current")


# create the controller signal
Kp, Kv = sm.symbols("Kp, Kv")
u = (Kv*(q1d_desired - q1d_current)) + (Kp*(q1_desired - q1_current))
u


# create the constant C1
C1, *_ = TM.subs(vals)/q1dd  # the term that describes the 
C1


# create the constant C2
R = 7.5        # ohms
kt = 0.891514  # Nm/amp, from manufacture's page
C2 = R/kt
C2


# calculate max Kp value
V_min = 3                   # voltage required for motor to actually move based on experimentation
CPR = 1000                  # counts per revolution
rad_min = 2*pi/1000         # smallest detectable error 
Kp_max = V_min/(rad_min*2)  # the 2 allows it to be +- one count from desired position
Kp_max


# create all the callables
simulation_parameters = {
    q1_desired: 3,
    q1d_desired: 0,
    Kp: int(Kp_max),
    Kv: sqrt(477*8),
}


# create controller callable
controller_signal = sm.lambdify([q1_current, q1d_current], u.subs(simulation_parameters), modules="numpy")
controller_signal(12, 0)


# create motor torque callable to emulate q1dd
rpms = (0, 40, 50)
taus = (0.980665, 0.4413, 0)
slope1, intercept1 = np.polyfit(rpms[0:2], taus[0:2], 1)
slope2, intercept2 = np.polyfit(rpms[1:], taus[1:], 1)
def tau_motor(x):
    if rpms[0] <= x <= rpms[1]:
        return (slope1*x) + intercept1
    elif rpms[1] <= x <= rpms[2]:
        return (slope2*x) + intercept2
    else:
        return 0


# create the RK45 simulation
def model(t, x):    
    # unpack the vector
    q1, q1d = x
    
    # calculate the controller signal
    u = controller_signal(q1, q1d)
    
    # get globals
    global C1
    
    # calculate tau
    tau_desired = abs(u*C1)    
    
    # calculate RPM
    RPM = abs(q1d*30/pi)
        
    # check if desired torque is in the opposite direction of RPM
    if sign(u) get_ipython().getoutput("= sign(q1d):")
        tau_max = tau_motor(0)
    else:
        tau_max = tau_motor(RPM)
    
    # determine value of tau
    tau = tau_desired if tau_desired < tau_max else tau_max
        
    # calculate the acceleration due to torque
    q1dd = tau*sign(u)/C1
    
    # return the xdots
    x1d = q1d
    x2d = q1dd
    
    return (x1d, x2d)

# create the initial conditions
simulator = RK45(model, 0, [0, 0], 3.5, max_step=0.05)
results = array([0, 0])
while True:
    try:
        simulator.step()
    except:
        break
    else:
        results = vstack([results, (simulator.t, simulator.y[0])])


plt.plot(results[:, 0], results[:, 1])
plt.grid()
plt.xlabel("time [s]")
plt.ylabel("q1 [rad]")
plt.title("Joint Space Controller");


print(results[-50:])


# create all the callables
simulation_parameters = {
    q1_desired: 12,
    q1d_desired: 0,
    Kp: int(Kp_max),
    Kv: 100
}


# create controller callable
controller_signal = sm.lambdify([q1_current, q1d_current], u.subs(simulation_parameters), modules="numpy")
controller_signal(12, 0)


# create the RK45 simulation
def model(t, x):    
    # unpack the vector
    q1, q1d = x
    
    # calculate the controller signal
    u = controller_signal(q1, q1d)
    
    # calculate tau
    tau_desired = abs(u)    
    
    # calculate RPM
    RPM = abs(q1d*30/pi)
        
    # check if desired torque is in the opposite direction of RPM
    if sign(u) get_ipython().getoutput("= sign(q1d):")
        tau_max = tau_motor(0)
    else:
        tau_max = tau_motor(RPM)
    
    # determine value of tau
    tau = tau_desired if tau_desired < tau_max else tau_max
    
    # get globals
    global C1
        
    # calculate the acceleration due to torque
    q1dd = tau*sign(u)/C1
    
    # return the xdots
    x1d = q1d
    x2d = q1dd
    
    return (x1d, x2d)

# create the initial conditions
simulator = RK45(model, 0, [0, 0], 3, max_step=0.05)
results = array([0, 0])
while True:
    try:
        simulator.step()
    except:
        break
    else:
        results = vstack([results, (simulator.t, simulator.y[0])])


plt.plot(results[:, 0], results[:, 1])
plt.grid()
plt.xlabel("time [s]")
plt.ylabel("q1 [rad]")
plt.title("PD Controller");
