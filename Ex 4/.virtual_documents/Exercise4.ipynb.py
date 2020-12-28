get_ipython().run_line_magic("reset", " -f")


get_ipython().run_line_magic("matplotlib", " inline")


import sympy as sm
import sympy.physics.mechanics as me
from scipy.integrate import RK45
from scipy.optimize import fsolve
from numpy import (array, vstack, pi, sqrt, arccos, 
                   sin, cos, linalg, linspace, sign, polyfit)
from math import degrees
import matplotlib.pyplot as plt


# setup N frame
N = me.ReferenceFrame("N")
pN = me.Point("pN")
pN.set_vel(N, 0)


# set up constants
g, gr, L = sm.symbols("g, gr, L")


# set up independent variables
q1, q2, q4, q5, u1, u2, u4, u5 = me.dynamicsymbols("q1, q2, q4, q5, u1, u2, u4, u5")
q1d, q2d, q4d, q5d, u1d, u2d, u4d, u5d = me.dynamicsymbols("q1, q2, q4, q5, u1, u2, u4, u5", 1)


# set up dependent variables which are the gears
q0, q3 = -q1*gr, -q4*gr
u0, u3 = -u1*gr, -u4*gr
q0d, q3d = -q1d*gr, -q4d*gr
u0d, u3d = -u1d*gr, -u4d*gr


# create dict for use when generating ccode
var2sym = {
    # positions
    q1: sm.Symbol("q1"),
    q2: sm.Symbol("q2"),
    q4: sm.Symbol("q4"),
    q5: sm.Symbol("q5"),
    # position derivatives
    q1d: sm.Symbol("q1d"),
    q2d: sm.Symbol("q2d"),
    q4d: sm.Symbol("q4d"),
    q5d: sm.Symbol("q5d"),
    # velocities
    u1: sm.Symbol("u1"),
    u2: sm.Symbol("u2"),
    u4: sm.Symbol("u4"),
    u5: sm.Symbol("u5"),
    # velocity derivatives
    u1d: sm.Symbol("u1d"),
    u2d: sm.Symbol("u2d"),
    u4d: sm.Symbol("u4d"),
    u5d: sm.Symbol("u5d"),
}


# set up the motor bodies
right_motor = me.Body(
    name = "motorR",
    frame = N.orientnew("R", "axis", (q0, N.z)),
    masscenter = pN.locatenew("motorR_COM", 0)
)
left_motor = me.Body(
    name = "motorL",
    frame = N.orientnew("L", "axis", (q3, N.z)),
    masscenter = pN.locatenew("motorL_COM", -L*N.x)
)


# set up the link bodies
bodyA = me.Body(
    name = "bodyA",
    frame = N.orientnew("A", "axis", (q1, N.z))
)
bodyB = me.Body(
    name = "bodyB",
    frame = bodyA.frame.orientnew("B", "axis", (q2, bodyA.frame.z))
)
bodyC = me.Body(
    name = "bodyC",
    frame = N.orientnew("C", "axis", (q4, N.z))
)
bodyD = me.Body(
    name = "bodyD",
    frame = bodyC.frame.orientnew("D", "axis", (q5, bodyC.frame.z))
)


# define the link lengths
bodyA.L1, bodyA.L2 = sm.symbols("LA_1, LA_2")
bodyB.L1, bodyB.L2 = sm.symbols("LB_1, LB_2")
bodyC.L1, bodyC.L2 = sm.symbols("LC_1, LC_2")
bodyD.L1, bodyD.L2 = sm.symbols("LD_1, LD_2")

# define the COM lengths
bodyA.L = sm.symbols("LA_COM")
bodyB.L = sm.symbols("LB_COM")
bodyC.L = sm.symbols("LC_COM")
bodyD.L = sm.symbols("LD_COM")


# set all the positions
bodyA.masscenter = pN.locatenew("pA", bodyA.L*bodyA.frame.x)
bodyB.masscenter = pN.locatenew("pB", (bodyA.L1*bodyA.frame.x) + (bodyB.L*bodyB.frame.x))
bodyC.masscenter = left_motor.masscenter.locatenew("pC", bodyC.L*bodyC.frame.x)
bodyD.masscenter = left_motor.masscenter.locatenew("pD", (bodyC.L1*bodyC.frame.x) + (bodyD.L*bodyD.frame.x))


# set all the velocities
right_motor.masscenter.set_vel(N, 0)
left_motor.masscenter.set_vel(N, 0)
bodyA.masscenter.set_vel(N, bodyA.masscenter.pos_from(pN).dt(N))
bodyB.masscenter.set_vel(N, bodyB.masscenter.pos_from(pN).dt(N))
bodyC.masscenter.set_vel(N, bodyC.masscenter.pos_from(pN).dt(N))
bodyD.masscenter.set_vel(N, bodyD.masscenter.pos_from(pN).dt(N))


# set all the accelerations
bodyA.masscenter.set_acc(N, bodyA.masscenter.vel(N).dt(N))
bodyB.masscenter.set_acc(N, bodyB.masscenter.vel(N).dt(N))
bodyC.masscenter.set_acc(N, bodyC.masscenter.vel(N).dt(N))


# define the independent coordinates
q_ind = [q1, q4]
q_dep = [q2, q5]
u_ind = [u1, u4]
u_dep = [u2, u5]

# and the kinematic differential equations
kde = [
    [q1d - u1],
    [q2d - u2],
    [q4d - u4],
    [q5d - u5],
]

# and create the vector loop equations
zero = pN.locatenew(
    "zero", 
    bodyA.L1*bodyA.frame.x + bodyB.L1*bodyB.frame.x - \
    bodyC.L1*bodyC.frame.x - bodyD.L1*bodyD.frame.x - L*N.x
)

pos_constraints = [zero.pos_from(pN).dot(N.x), 
                   zero.pos_from(pN).dot(N.y)]
vel_constraints = [zero.pos_from(pN).dt(N).dot(N.x), 
                   zero.pos_from(pN).dt(N).dot(N.y)]


KM = me.KanesMethod(
    frame = N, 
    q_ind = q_ind,
    u_ind = u_ind,
    kd_eqs = kde,
    q_dependent = q_dep,
    configuration_constraints = pos_constraints,
    u_dependent = u_dep,
    velocity_constraints = vel_constraints
)


# define the loads
tau1, tau2 = sm.symbols("tau1, tau2")
loads = [
    # point/frame, force/moment
    (right_motor.masscenter, right_motor.mass*g*-N.z),
    (left_motor.masscenter, left_motor.mass*g*-N.z),
    (bodyA.masscenter, bodyA.mass*g*-N.z),
    (bodyB.masscenter, bodyB.mass*g*-N.z),
    (bodyC.masscenter, bodyC.mass*g*-N.z),
    (bodyD.masscenter, bodyD.mass*g*-N.z),
    (bodyA.frame, tau1*right_motor.frame.z),
    (bodyC.frame, tau2*left_motor.frame.z),
]
# define the bodies
bodies = [right_motor, bodyA, bodyB, left_motor, bodyC, bodyD]


Fr, Fr_star = KM.kanes_equations(bodies, loads)
Fr.shape, Fr_star.shape


in2m = lambda x: x*25.4/1000
vals = {
    # constants
    g: 9.81,               # N*m
    L: in2m(4),            # m
    gr: 144,               # a.u.
    # joint c-c
    bodyA.L1: in2m(2.4),   # m
    bodyC.L1: in2m(2.4),   # m
    bodyB.L1: in2m(3.2),   # m
    bodyD.L1: in2m(3.2),   # m
    # com properties
    bodyA.L: 0.027,        # m
    bodyC.L: 0.027,        # m
    bodyB.L: 0.04,         # m
    bodyD.L: 0.04,         # m
    # masses
    bodyA.mass: 0.007,     # kg
    bodyC.mass: 0.007,     # kg 
    bodyB.mass: 0.008,     # kg
    bodyD.mass: 0.008,     # kg
    # moments of inertia 
    "bodyA_izz": 6.39e-6,  # kg m^2
    "bodyC_izz": 6.39e-6,  # kg m^2
    "bodyB_izz": 1.01e-5,  # kg m^2
    "bodyD_izz": 1.01e-5,  # kg m^2
    "motorR_izz": 3.3e-7,  # kg m^2
    "motorL_izz": 3.3e-7,  # kg m^2
}


# create python callable for qdd vec based on forcing matrix
FM = KM.forcing.subs(vals)
MM = KM.mass_matrix.subs(vals)
FM = me.msubs(FM, {q1d: u1, q2d: u2, q4d: u4, q5d: u5})
FM_func = sm.lambdify([q1, q2, q4, q5, u1, u2, u4, u5, tau1, tau2], FM, modules="numpy")
MM_func = sm.lambdify([q1, q2, q4, q5], MM, modules="numpy")
FM_func(*[0.1]*10), MM_func(*[0.1]*4)


# create container to capture auto generated c code
c_code = {}
# save off expression for tau1 and tau2 in c code
c_code["tau1"] = sm.ccode(me.msubs(Fr_star[0], var2sym))
c_code["tau2"] = sm.ccode(me.msubs(Fr_star[1], var2sym))


# create python callable for tau1 and tau2
tau_func = sm.lambdify([q1, q2, q4, q5, u1, u2, u4, u5, u1d, u2d, u4d, u5d], Fr_star.subs(vals))
tau_func(*[0.1]*12)


# create the position matrix
pE = pN.locatenew("pE", (bodyA.L1*bodyA.frame.x) + (bodyB.L1*bodyB.frame.x))
PNE = (pE.pos_from(pN).to_matrix(N))
PNE = PNE.row_del(-1)
c_code["x_current"] = sm.ccode(me.msubs(PNE[0], var2sym))
c_code["y_current"] = sm.ccode(me.msubs(PNE[1], var2sym))


# create callable for the simulation
PNE = PNE.subs(vals)
PNE_func = sm.lambdify([q1, q2], PNE, modules="sympy")
PNE_func(0, 0)


# create velocity matrix
VNE = pE.pos_from(pN).dt(N).to_matrix(N)
VNE = VNE.row_del(-1)
# save off the expressions in c-code
c_code["xd_current"] = sm.ccode(me.msubs(VNE[0], var2sym))
c_code["yd_current"] = sm.ccode(me.msubs(VNE[1], var2sym))


# create the python callable
VNE = VNE.subs(vals)
VNE_func = sm.lambdify([q1, q2, q1d, q2d], VNE, modules="sympy")
VNE_func(0, 0, 1, 0)


# create the task space controller
Kp, Kv = sm.symbols("Kp, Kv")
X_current = sm.Matrix(sm.symbols("x_current, y_current"))
Xd_current = sm.Matrix(sm.symbols("xd_current, yd_current"))
X_desired = sm.Matrix(sm.symbols("x_desired, y_desired"))
Xd_desired = sm.Matrix(sm.symbols("xd_desired, yd_desired"))
u = (Kv*(Xd_desired - Xd_current)) + (Kp*(X_desired - X_current))
u


# save off the controller code in c
c_code["xdd_desired"] = sm.ccode(me.msubs(u[0], var2sym))
c_code["ydd_desired"] = sm.ccode(me.msubs(u[1], var2sym))


# get the positions
pE_rhs = right_motor.masscenter.locatenew("pE_rhs", (bodyA.L1*bodyA.frame.x) + (bodyB.L1*bodyB.frame.x))
pE_lhs = left_motor.masscenter.locatenew("pE_lhs", (bodyC.L1*bodyC.frame.x) + (bodyD.L1*bodyD.frame.x))

# create the 2x2 matrix
M_rhs = pE_rhs.pos_from(pN).dt(N).to_matrix(N)
M_lhs = pE_lhs.pos_from(pN).dt(N).to_matrix(N)
M_rhs = M_rhs.subs(vals)
M_lhs = M_lhs.subs(vals)
M_rhs = M_rhs.row_del(-1)
M_lhs = M_lhs.row_del(-1)
M = sm.Matrix.vstack(M_rhs, M_lhs)

# create the jacobains
X_rhs = sm.Matrix([q1d, q2d])
X_lhs = sm.Matrix([q4d, q5d])
X = sm.Matrix([q1d, q2d, q4d, q5d])
J_rhs = M_rhs.jacobian(X_rhs)
J_lhs = M_lhs.jacobian(X_lhs)
J = M.jacobian(X)
J_rhs.shape, J_lhs.shape, J.shape


# extract the Jacobina c++ code
rows, cols = J.shape
for r in range(rows):
    for c in range(cols):
        c_code[f"J({r},{c})"] = J[r, c]


# create jocobain python callable
J_func = sm.lambdify([q1, q2, q4, q5], J, modules="sympy")
J_func(1, 1, 1, 1)


# create the jacobian derivative
Jd = sm.diff(J, 't')

# collect the jacobian derivative in c code
rows, cols = Jd.shape
for r in range(rows):
    for c in range(cols):
        c_code[f"J({r},{c})"] = Jd[r, c]


# create the jacobian derivative python callable
Jd_func = sm.lambdify([q1, q2, q4, q5, q1d, q2d, q4d, q5d], Jd, modules="sympy")
Jd_func(1, 1, 1, 1, 1, 1, 1, 1)


# create length vectors
p1 = right_motor.masscenter.locatenew("p1", bodyA.L1*bodyA.frame.x)
p2 = left_motor.masscenter.locatenew("p2", bodyC.L1*bodyC.frame.x)

# calculate gamma using the law of cosines
W = p2.pos_from(p1)
gamma = sm.acos(W.magnitude()/(2*bodyB.L1))

# collect the exprssion for gamma in c code
c_code["gamma"] = sm.ccode(me.msubs(gamma, var2sym))
gamma


# substitute in numerical values
gamma = gamma.subs(vals)

# create gamma python callable
gamma_func = sm.lambdify([q1, q4], gamma, modules="numpy")
gamma_func(pi/2, pi/2)


# create expression for alpha1
alpha1 = sm.acos(W.dot(-p1.pos_from(right_motor.masscenter))/(W.magnitude()*p1.pos_from(right_motor.masscenter).magnitude()) )
# save expression for alpha1 in c code
c_code["alpha1"] = sm.ccode(me.msubs(alpha1, var2sym))
alpha1


# sub in numbers
alpha1 = alpha1.subs(vals)
alpha1


# create expression for alpha2
alpha2 = q1 - q4 - alpha1 + sm.pi
# save off expression in c code
c_code["alpha2"] = sm.ccode(me.msubs(alpha2, var2sym))
alpha2


# create python callables for alpha1 and alpha2
alpha_func = sm.lambdify([q1, q4], [alpha1, alpha2], modules="numpy")
alpha_func(pi/2, pi/2)


# solve for q2d and q5d using vector loop constraint
M = zero.pos_from(pN).dt(N).to_matrix(N)
M = M.subs(vals)
sol = sm.solve([M], [q2d, q5d])
sol


# save off the expressions in c code
eq1 = me.msubs(sol[q2d], var2sym)
eq2 = me.msubs(sol[q5d], var2sym)
c_code["q2d"] = eq1
c_code["q5d"] = eq2


# create python callable that solves for q2d, and q5d
qd_func = sm.lambdify([q1, q2, q4, q5, q1d, q4d], [sol[q2d], sol[q5d]], modules="numpy")
qd_func(pi/2, 0.68, pi/2, -0.68, 1, 1)


# create motor torque python callable to drive q1dd and q4dd
rpms = (0, 40)
taus = (0.980665, 0.4413)
slope, intercept = polyfit(rpms, taus, 1)

def tau_motor(x):
    if rpms[0] <= x <= rpms[1]:
        return (slope*x) + intercept
    else:
        return 0


# create the python callable for the target point
R, r, d = in2m(0.5), in2m(0.1), in2m(0.65)
x_offset, y_offset = -in2m(2), in2m(4)
xh = lambda theta: float(((R-r)*cos(theta)) + (d*cos((R-r)*(theta/r))) + x_offset)
yh = lambda theta: float(((R-r)*sin(theta)) - (d*sin((R-r)*(theta/r))) + y_offset)


# plot the target
m2in = lambda x: x*1000/25.4
thetas = linspace(0, 2*pi, 75)
x_vals = [xh(t) for t in thetas]
y_vals = [yh(t) for t in thetas]
plt.plot([m2in(v) for v in x_vals], [m2in(v) for v in y_vals], '.')
plt.grid()
plt.title("control space target")
plt.xlabel("N.x [m]")
plt.ylabel("N.y [m]");


# create the ODE45 to run the simulation
# x1 = q1,    x1d = q1d
# x2 = q4,    x2d = q4d
# x3 = q1d,   x3d = q1dd
# x4 = q4d,   x4d = q4dd

# declare controller variables
KP = 100
KV = sqrt(KP*8)

# declare globals for point tracking
X_CURRENT, Y_CURRENT = 0, 0
X_DESIRED, Y_DESIRED = 0, 0

def model(t, x):
    # unpack the vector
    q1, q4, q1d, q4d = x
    
    # get gamma, alpha1, and alpha2
    gamma = gamma_func(q1, q4)
    alpha1, alpha2 = alpha_func(q1, q4)
    
    # solve for q2, and q5
    q2 = pi - (alpha1 + gamma)
    q5 = pi + (alpha2 + gamma)
     
    # solve for the velocities
    q2d, q5d = qd_func(q1, q2, q4, q5, q1d, q4d)
    
    # calculate point E position
    x, y = PNE_func(q1, q2)
    
    # assign to global namespace
    global X_CURRENT, Y_CURRENT
    X_CURRENT, Y_CURRENT = x, y
    
    # calculate the speeds
    xd, yd = VNE_func(q1, q2, q1d, q2d)
    
    # calculate the controller output
    global X_DESIRED, Y_DESIRED, KP, KV
    u = (KV*(sm.Matrix([0 - xd, 0 - yd]))) + \
        (KP*(sm.Matrix([X_DESIRED - x, Y_DESIRED - y])))
    
    # create task space acceleration vector
    Xdd = sm.Matrix.vstack(u, u)
    
    # calculate the jacobian
    J = J_func(q1, q2, q4, q5)
    
    # calculate the angular velocities
    Xd = sm.Matrix([xd, yd, xd, yd])
    q1d, q2d, q4d, q5d = J.inv()*Xd
    
    # calculate the desired accelerations
    Jd = Jd_func(q1, q2, q4, q5, q1d, q2d, q4d, q5d)
    q1dd, q2dd, q4dd, q5dd = J.inv()*(Xdd - (Jd*sm.Matrix([q1d, q2d, q4d, q5d])))
    
    # calculate desired torques
    TM = tau_func(q1, q2, q4, q5, q1d, q2d, q4d, q5d, q1dd, q2dd, q4dd, q5dd)
    tau1_desired, tau2_desired = abs(TM[0][0]), abs(TM[1][0])
    
    # calculate max torques possible from motors
    if sign(q1dd) get_ipython().getoutput("= sign(q1d):")
        tau1_max = tau_motor(0)
    else:
        tau1_max = tau_motor(abs(q1d*30/pi))
        
    if sign(q4dd) get_ipython().getoutput("= sign(q4d):")
        tau2_max = tau_motor(0)
    else:
        tau2_max = tau_motor(abs(q4d*30/pi))
        
    # constrain acclerations based on torques
    if (tau1_desired > tau1_max) or (tau2_desired > tau2_max):
        MM = MM_func(q1, q2, q4, q5)
        FM = FM_func(q1, q2, q4, q5, q1d, q2d, q4d, q5d, tau1_max*sign(q1dd), tau2_max*sign(q4dd))
        Qdd = linalg.inv(MM)*FM
        q1dd, q2dd, q4dd, q5dd = Qdd[0][0], Qdd[1][0], Qdd[2][0], Qdd[3][0]
                    
    # assign the xdots
    x1d = q1d
    x2d = q4d
    x3d = q1dd
    x4d = q4dd
    
    # return the values
    return (x1d, x2d, x3d, x4d)


conditions = array([pi/2, pi/2, 0, 0])  # stores state between targets
times = array([0])                      # records the time state between targets
results = array([0, 0, *conditions])    # records the x, y, q1, q4, q1d, & q4d points for the entire simulation
eps = 0.002                             # error threashold for event based simulation conclusion


def run_interval(theta, times, conditions, results):
    # set up this leg of the simulation
    simulator = RK45(model, times[0], conditions, t_bound=0.1, max_step = 0.1)
    # set the desired points 
    global X_DESIRED, Y_DESIRED
    X_DESIRED = xh(theta)
    Y_DESIRED = yh(theta)
    # run this leg of the simulation
    while True:
        # attempt to step
        try:
            simulator.step()
        # if the solver is finished
        except Exception as e:
            break
        # if it's close enough to the target
        if linalg.norm([float(X_DESIRED - X_CURRENT), float(Y_DESIRED - Y_CURRENT)]) <= eps:
            break
        # store the current state
        conditions = simulator.y
        # append step results to simulation results
        times = vstack([times, simulator.t])
        results = vstack([results, (X_CURRENT, Y_CURRENT, *conditions)])
            
    # return this leg of the simulation
    return theta, times, conditions, results


# run all the intervals
theta = 0
while theta < 7:
    # put in the last conditions
    theta, times, conditions, results = run_interval(theta, times, conditions, results)
    theta += 0.1


# curate the points
Ex = [m2in(x) for x in results[1::, 0]]
Ey = [m2in(y) for y in results[1::, 1]]
Rx = [m2in(vals[bodyA.L1]*cos(q)) for q in results[1::, 2]]
Ry = [m2in(vals[bodyA.L1]*sin(q)) for q in results[1::, 2]]
Lx = [m2in(vals[bodyC.L1]*cos(q) - vals[L]) for q in results[1::, 3]]
Ly = [m2in(vals[bodyC.L1]*sin(q)) for q in results[1::, 3]]

# plot the results
for i, _ in enumerate(Ex):
    plt.plot([0, Rx[i], Ex[i]],[0, Ry[i], Ey[i]])                # right arm
    plt.plot([m2in(-vals[L]), Lx[i], Ex[i]], [0, Ly[i], Ey[i]])  # left arm

plt.grid()
plt.title("control space simulation")
plt.xlabel("N.x [in]")
plt.ylabel("N.y [in]");


# define paths to end effector
rhs = right_motor.masscenter.locatenew("rhs", bodyA.L1*bodyA.frame.x + bodyB.L1*bodyB.frame.x)
lhs = left_motor.masscenter.locatenew("lhs", bodyC.L1*bodyC.frame.x + bodyD.L1*bodyD.frame.x)
# get the matricies
rhs_M = rhs.pos_from(pN).to_matrix(N).row_del(-1)
lhs_M = lhs.pos_from(pN).to_matrix(N).row_del(-1)
# sub in numerical values for constants
rhs_M = rhs_M.subs(vals)
lhs_M = lhs_M.subs(vals)
# create the equations
x_desired, y_desired = sm.symbols("x_desired, y_desired")
eq1 = rhs_M[0] - x_desired
eq2 = rhs_M[1] - y_desired
eq3 = lhs_M[0] - x_desired
eq4 = lhs_M[1] - y_desired
# create the callable
kin_eq_func = sm.lambdify([q1, q2, q4, q5, x_desired, y_desired], [eq1, eq2, eq3, eq4], modules="scipy")
kin_eq_func(0, 0, 0, 0, 1, 4)


# put in values for x and y
def model(x, args):
    # unpack the q vec
    q1, q2, q4, q5 = x
    # get the desired points
    x_des, y_des = args
    # return the output for the given inputs
    return kin_eq_func(q1, q2, q4, q5, x_des, y_des)

# loop over each desires point
q_desireds = array([0, 0, 0, 0])
conditions = array([pi/2, 0.1, pi/2, -0.1])
for x, y in zip(x_vals, y_vals):
    sol = fsolve(model, x0 = conditions, args=[x, y])
    q_desireds = vstack([q_desireds, sol])
    conditions = sol
    
# drop the initializer values
q_desireds = q_desireds[1::]


# create the ODE45 to run the simulation
# x1 = q1,    x1d = q1d
# x2 = q4,    x2d = q4d
# x3 = q1d,   x3d = q1dd
# x4 = q4d,   x4d = q4dd

# declare controller variables
KP = 100
KV = sqrt(KP*8)

# declare globals for point tracking
Q1_CURRENT, Q2_CURRENT, Q4_CURRENT, Q5_CURRENT = 0, 0, 0, 0
Q1_DESIRED, Q2_DESIRED, Q4_DESIRED, Q5_DESIRED = 0, 0, 0, 0

def model(t, x):
    # unpack the vector
    q1, q4, q1d, q4d = x
    
    # get gamma, alpha1, and alpha2
    gamma = gamma_func(q1, q4)
    alpha1, alpha2 = alpha_func(q1, q4)
    
    # solve for q2, and q5
    q2 = pi - (alpha1 + gamma)
    q5 = pi + (alpha2 + gamma)
    
    # assign to to globals
    global Q1_CURRENT, Q2_CURRENT, Q4_CURRENT, Q5_CURRENT
    Q1_CURRENT, Q2_CURRENT, Q4_CURRENT, Q5_CURRENT = q1, q2, q4, q5   
     
    # solve for the velocities
    q2d, q5d = qd_func(q1, q2, q4, q5, q1d, q4d)
    
    # calculate the controller output
    global Q1_DESIRED, Q2_DESIRED, Q4_DESIRED, Q5_DESIRED, KP, KV
    u = (KV*(sm.Matrix([0 - q1d, 0 - q2d, 0 - q4d, 0 - q5d]))) + \
        (KP*(sm.Matrix([Q1_DESIRED - q1, Q2_DESIRED - q2, Q4_DESIRED - q4, Q5_DESIRED - q5])))
    
    # assign to vars
    q1dd, q4dd = u[0], u[2]
    
    # calculate max torques possible from motors
    if sign(q1dd) get_ipython().getoutput("= sign(q1d):")
        tau1_max = tau_motor(0)
    else:
        tau1_max = tau_motor(abs(q1d*30/pi))
        
    if sign(q4dd) get_ipython().getoutput("= sign(q4d):")
        tau2_max = tau_motor(0)
    else:
        tau2_max = tau_motor(abs(q4d*30/pi))
        
    # constrain acclerations based on possible torque
    MM = MM_func(q1, q2, q4, q5)
    FM = FM_func(q1, q2, q4, q5, q1d, q2d, q4d, q5d, tau1_max*sign(q1dd), tau2_max*sign(q4dd))
    Qdd = linalg.inv(MM)*FM
    q1dd_max, q4dd_max = Qdd[0][0],  Qdd[2][0], 
    if abs(q1dd) > abs(q1dd_max):
        q1dd = abs(q1dd_max)*sign(q1dd)
    if abs(q4dd) > abs(q4dd_max):
        q4dd = abs(q4dd_max)*sign(q4dd)
           
    # assign the xdots
    x1d = q1d
    x2d = q4d
    x3d = q1dd
    x4d = q4dd
    
    # return the values
    return (x1d, x2d, x3d, x4d)


conditions = array([pi/2, pi/2, 0, 0])  # stores state between targets
times = array([0])                      # records the time state between targets
results = array([*conditions])          # records the, q1, q2, q4, & q5 for the entire simulation
eps = 1e-5                              # error threashold for event based simulation conclusion


def run_interval(count, times, conditions, results):
    # set up this leg of the simulation
    simulator = RK45(model, times[0], conditions, t_bound=5, max_step = 0.01)
    # set the desired points 
    global Q1_DESIRED, Q2_DESIRED, Q4_DESIRED, Q5_DESIRED, q_desireds
    Q1_DESIRED = q_desireds[count, 0]
    Q2_DESIRED = q_desireds[count, 1]
    Q4_DESIRED = q_desireds[count, 2]
    Q5_DESIRED = q_desireds[count, 3]
    # run this leg of the simulation
    while True:
        # attempt to step
        try:
            simulator.step()
        # if the solver is finished
        except Exception as e:
            break
        # if it's close enough to the target
        if (abs(float(Q1_DESIRED - Q1_CURRENT)) <= eps)  and (abs(float(Q4_DESIRED - Q4_CURRENT)) <= eps):
            break
        # store the current state
        conditions = simulator.y
        # append step results to simulation results
        times = vstack([times, simulator.t])
        results = vstack([results, (Q1_CURRENT, Q2_CURRENT, Q4_CURRENT, Q5_CURRENT)])
            
    # return this leg of the simulation
    return times, conditions, results


# run all the intervals
for count, _ in enumerate(q_desireds[1:], start=1):
    # put in the last conditions
    times, conditions, results = run_interval(count, times, conditions, results)


# curate the points
step = 100
Ex = [m2in(PNE_func(q1, q2)[0]) for q1, q2, *_ in results[1::step]]
Ey = [m2in(PNE_func(q1, q2)[1]) for q1, q2, *_ in results[1::step]]
Rx = [m2in(vals[bodyA.L1]*cos(q)) for q in results[1::step, 0]]
Ry = [m2in(vals[bodyA.L1]*sin(q)) for q in results[1::step, 0]]
Lx = [m2in(vals[bodyC.L1]*cos(q) - vals[L]) for q in results[1::step, 2]]
Ly = [m2in(vals[bodyC.L1]*sin(q)) for q in results[1::step, 2]]

# plot the results
for i, _ in enumerate(Ex):
    plt.plot([0, Rx[i], Ex[i]],[0, Ry[i], Ey[i]])                # right arm
    plt.plot([m2in(-vals[L]), Lx[i], Ex[i]], [0, Ly[i], Ey[i]])  # left arm

plt.grid()
plt.title("joint space simulation")
plt.xlabel("N.x [in]")
plt.ylabel("N.y [in]");


# print out c_code expressions
for key, value in c_code.items():
    print(f"{key} = {value}\n")


per_line = 5

print("x_vals")
for i, _ in enumerate(x_vals[::per_line]):
    i *= per_line
    print(*x_vals[i:i+5], "", sep=",")
    
print("y_vals")
for i, _ in enumerate(y_vals[::per_line]):
    i *= per_line
    print(*y_vals[i:i+5], "", sep=",")


print("q1_vals")
q1_vals = [q[0] for q in q_desireds]
for i, _ in enumerate(q1_vals[::per_line]):
    i *= per_line
    print(*q1_vals[i:i+5], "", sep=",")

print("q2_vals")
q2_vals = [q[1] for q in q_desireds]
for i, _ in enumerate(q2_vals[::per_line]):
    i *= per_line
    print(*q2_vals[i:i+5], "", sep=",")

print("q4_vals")
q4_vals = [q[2] for q in q_desireds]
for i, _ in enumerate(q4_vals[::per_line]):
    i *= per_line
    print(*q4_vals[i:i+5], "", sep=",")

print("q5_vals")
q5_vals = [q[3] for q in q_desireds]
for i, _ in enumerate(q5_vals[::per_line]):
    i *= per_line
    print(*q5_vals[i:i+5], "", sep=",")
