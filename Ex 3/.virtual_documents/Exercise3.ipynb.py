get_ipython().run_line_magic("reset", " -f")


get_ipython().run_line_magic("matplotlib", " inline")


import sympy as sm
import sympy.physics.mechanics as me
import matplotlib.pyplot as plt
from IPython.display import Latex, display
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
    bodyB.L: 0.02162775,     # m
    bodyB.mass: 0.00728203,  # kg 
    "bodyB_izz": 6.39e-06,   # kg*m^2
    "bodyA_izz": 3.3e-07,    # kg*m^2
    gear_ratio: 144,         # a.u
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
C1, *_ = TM.subs(vals)/q1dd  # kg*(m^2)
C1


# create the constant C2
R = 7.5        # ohms
kt = 0.891514  # Nm/amp, from manufacture's page
C2 = R/kt
C2


# determine the minimum voltage to get mapped to < 0.5
def arduino_map(x, in_min, in_max, out_min, out_max):
    return ((x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min)

vmin = 0.0235
print(arduino_map(vmin, 0, 12, 0, 255))

CPR = 1000
res = 2*pi/CPR
print(res)
Kp_max = vmin/(res*C2*C1)
print(Kp_max)


# create all the callables
simulation_parameters = {
    q1_desired: 3,
    q1d_desired: 0,
    Kp: int(Kp_max),
    Kv: sqrt(int(Kp_max)*8),
}


# create controller callable
controller_signal = sm.lambdify([q1_current, q1d_current], u.subs(simulation_parameters), modules="numpy")
controller_signal(3, 0)


# create motor torque callable to drive q1dd
rpms = (0, 40)
taus = (0.980665, 0.4413)
slope, intercept = polyfit(rpms, taus, 1)

def tau_motor(x):
    if rpms[0] <= x <= rpms[1]:
        return (slope*x) + intercept
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
        
# save the results to a file
with open("joint_space_3rad.csv", "w") as f:
    for r in results:
        f.write(f"{r[0]},{r[1]}\n")


# plot the results from the 3 second simulation and experimentation
fig, ax = plt.subplots(figsize=(7, 5))
with open("joint_space_3rad.csv", "r") as f:
    res = [l.split(",") for l in f.readlines()]
    times = [float(v[0]) for v in res]
    angles = [float(v[1]) for v in res]
    ax.plot(times, angles, label="3 rad simulated")
    
with open("joint_space_3rad_actual.csv", "r") as f:
    res = [l.split(",") for l in f.readlines()]
    times = [float(v[0])/1000 for v in res]
    angles = [float(v[1]) for v in res]
    ax.plot(times, angles, label="3 rad actual")
    
with open("joint_space_6rad.csv", "r") as f:
    res = [l.split(",") for l in f.readlines()]
    times = [float(v[0]) for v in res]
    angles = [float(v[1]) for v in res]
    ax.plot(times, angles, label="6 rad simulated")
    
with open("joint_space_6rad_actual.csv", "r") as f:
    res = [l.split(",") for l in f.readlines()]
    times = [float(v[0])/1000 for v in res]
    angles = [float(v[1]) for v in res]
    ax.plot(times, angles, label="6 rad actual")

with open("joint_space_12rad.csv", "r") as f:
    res = [l.split(",") for l in f.readlines()]
    times = [float(v[0]) for v in res]
    angles = [float(v[1]) for v in res]
    ax.plot(times, angles, label="12 rad simulated")
    
with open("joint_space_12rad_actual.csv", "r") as f:
    res = [l.split(",") for l in f.readlines()]
    times = [float(v[0])/1000 for v in res]
    angles = [float(v[1]) for v in res]
    ax.plot(times, angles, label="12 rad actual")


ax.grid()
ax.set_title("Joint Space Results")
ax.set_xlabel("time [s]")
ax.set_ylabel("angle [rad]")
ax.legend();


# create the points to connect
bodyB.W = bodyB.L/3
p1 = pN.locatenew("p1", 0.5*bodyB.W*bodyB.frame.y)
p2 = p1.locatenew("p2", bodyB.L*bodyB.frame.x)
p3 = p2.locatenew("p3", bodyB.W*-bodyB.frame.y)
p4 = p3.locatenew("p4", bodyB.L*-bodyB.frame.x)

# turn the points into a callable
points = [p1, p2, p3, p4, p1]
points = [p.pos_from(pN).to_matrix(N)[0:2] for p in points]
points = [[p.subs(vals) for p in tup] for tup in points]
points = [[p.subs(vals) for p in tup] for tup in points]
points_func = sm.lambdify([q1], points, modules="numpy")


for ang in results[0::5, 1]:
    points = points_func(ang)
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    plt.plot(x, y)
plt.grid()
plt.title("Motion Capture Plot 3 rad")
plt.xlabel("N.x [m]");
plt.ylabel("N.y [m]");


CPR = 1000
res = 2*pi/CPR
print(res)
Kp_max = vmin/(res*C2)
print(Kp_max)


# create all the callables
simulation_parameters = {
    q1_desired: 6,
    q1d_desired: 0,
    Kp: int(Kp_max),
    Kv: 0.1
}


# create controller callable
controller_signal = sm.lambdify([q1_current, q1d_current], u.subs(simulation_parameters), modules="numpy")
controller_signal(6, 0)


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
simulator = RK45(model, 0, [0, 0], 3, max_step=0.1)
results = array([0, 0])
while True:
    try:
        simulator.step()
    except:
        break
    else:
        results = vstack([results, (simulator.t, simulator.y[0])])
        
# save the results to a file
with open("PD_10.csv", "w") as f:
    for r in results:
        f.write(f"{r[0]},{r[1]}\n")


# plot the results from the 3 second simulation and experimentation
fig, ax = plt.subplots(figsize=(7, 5))
with open("PD_10.csv", "r") as f:
    res = [l.split(",") for l in f.readlines()]
    times = [float(v[0]) for v in res]
    angles = [float(v[1]) for v in res]
    ax.plot(times, angles, label="Kv=10 simulated")
    
with open("PD_10_actual.csv", "r") as f:
    res = [l.split(",") for l in f.readlines()]
    times = [float(v[0])/1000 for v in res]
    angles = [float(v[1]) for v in res]
    ax.plot(times, angles, label="Kv=10 actual")
    
with open("PD_100.csv", "r") as f:
    res = [l.split(",") for l in f.readlines()]
    times = [float(v[0]) for v in res]
    angles = [float(v[1]) for v in res]
    ax.plot(times, angles, label="Kv=100 simulated")
    
with open("PD_100_actual.csv", "r") as f:
    res = [l.split(",") for l in f.readlines()]
    times = [float(v[0])/1000 for v in res]
    angles = [float(v[1]) for v in res]
    ax.plot(times, angles, label="Kv=100 actual")
    
with open("PD_1000.csv", "r") as f:
    res = [l.split(",") for l in f.readlines()]
    times = [float(v[0]) for v in res]
    angles = [float(v[1]) for v in res]
    ax.plot(times, angles, label="Kv=1000 simulated")
    
with open("PD_1000_actual.csv", "r") as f:
    res = [l.split(",") for l in f.readlines()]
    times = [float(v[0])/1000 for v in res]
    angles = [float(v[1]) for v in res]
    ax.plot(times, angles, label="Kv=1000 actual")

ax.grid()
ax.set_title("PD Results")
ax.set_xlabel("time [s]")
ax.set_ylabel("angle [rad]")
ax.legend();


from mathtools.symbolic import sym2poly
from mathtools.root_locus import *
from numpy import linspace, cos, angle


Kp, Kv, x1, s = sm.symbols("Kp, Kv, x1, s")
PD = Kp + (Kv*s)
CLTF = PD/((s**2) - ((x1)*s) + PD)
CLTF
num, den = sm.fraction(CLTF)
num, den
vals = {Kp: 1, Kv: 0.2, x1: 2.29}
num, *_ = sym2poly(num.subs(vals), var=s)
den, *_ = sym2poly(den.subs(vals), var=s)
num, den
gains = linspace(0,20,1000)
quick_plot(num, den, gains)


tf = transfer_function(num, den)
res = compute_roots(tf, gains)
[print(g, cos(angle(r[0]))) for g, r in zip(gains, res)];
