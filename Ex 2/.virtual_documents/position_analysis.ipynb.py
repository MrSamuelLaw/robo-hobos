get_ipython().run_line_magic("reset", " -f")
import sympy as sm
import sympy.physics.mechanics as me
from IPython.display import display, Latex


def show(eq_lhs: str, eq_rhs: str):
    eq = f"${eq_lhs} = {eq_rhs}$"
    return diplay(Latex(eq))


# define variables
q1 = me.dynamicsymbols("q1")

# define constants
LA = sm.symbols("LA")


# define the frames
N = me.ReferenceFrame("N")
A = N.orientnew("A", "axis", (q1, N.z))


# define the points
pN = me.Point("pN")
pA = pN.locatenew("pA", LA*A.x)


# create the expresson
print(me.vlatex(N.dcm(A).subs({q1: "q1"})))
