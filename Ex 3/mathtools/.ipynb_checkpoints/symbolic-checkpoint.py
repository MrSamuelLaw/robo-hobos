#!/usr/bin/env python3


import sympy as sym
import numpy as np
from IPython.display import Markdown as md


class sym2poly():
    """Sym2poly class used to convert polynomial coefficients
    to vectors of coefficients. turned into a class to allow
    for easy extension as needed."""

    def __new__(self, *polynomials, var=None, floater=False):
        """takes a symbolic polynomial and returns
        a vector of the leading coefficients

        args:
            polynomials: Expression(s) for which you wish to get the
                        coefficients
            var: Polynomial variable such as s in s**2 + s**1 + s**0
                must be specified if polynomial has symbolic coefficients
        """

        results = []
        for poly in polynomials:
            p = sym.simplify(poly)
            # if the variable is defined
            if var:
                p = sym.Poly(p, var)
                c = p.all_coeffs()
                results.append(c)
            # if variable is not defined
            else:
                try:
                    p = sym.Poly(p)
                    c = p.all_coeffs()
                # handles s**0 poly's
                except sym.polys.GeneratorsNeeded:
                    c = [sym.Float(p)]
                finally:
                    results.append(c)

        # try to convert to floats if floater set to true
        if floater:
            results = [[float(x) for x in item] for item in results]

        return results


def show(left, right, *the_rest):
    """returns a markdown object print out
    equations, both symbolic and numeric in
    Jupyter notebook"""
    if the_rest:
        return md(f"${left} = {right}{the_rest}$")
    return md(f"${left} = {right}$")


def round_expr(expr, num_digits):
    try:
        val = expr.xreplace({n: round(n, num_digits) for n in expr.atoms(sym.Number)})
    except AttributeError:
        val = expr
    finally:
        return val


class RouthArray():
    """Symilar to the RHtable func, however
    it provides additional functionality."""

    def __init__(self, CE, var):
        """create routh array and assign to self"""
        # define printout rounding
        self.print_rounding = 3

        CE = sym.simplify(CE)
        CE = sym.expand(CE)
        CE = sym.simplify(CE)
        CE = sym.collect(CE, var)

        self.CE = CE

        rows = [[] for i in range(sym.degree(CE, var)+1)]  # initialize vector
        for count, n in enumerate(range(sym.degree(CE, var), -1, -1)):
            if not count % 2:  # it's an odd number append to first row
                rows[0].append(CE.coeff(var, n))
            else:
                rows[1].append(CE.coeff(var, n))

        # append additional zeros
        while len(rows[0]) != len(rows[1]):
            rows[1].append(0)
        rows[0].append(0)
        rows[1].append(0)

        # fill up the rest of the rows
        for row, r in enumerate(rows[2:], start=2):  # start at 2nd row
            x = []
            for col in range(len(rows[row-1]) - 1):
                try:
                    val = ((rows[row-1][0] * rows[row-2][col+1]) -
                    (rows[row-2][col] * rows[row-1][col+1]))/rows[row-1][0]
                except ZeroDivisionError:
                    val = sym.nan
                finally:
                    x.append(val)
            rows[row].extend(x)

        self.array = rows

    def __str__(self):
        """Returns  a sympy matrix to
        print so everything looks super pretty"""

        rows = [[str(round_expr(y, self.print_rounding)) for y in x] for x in self.array]
        rows = [' | '.join(x) for x in rows]
        rows = '\n'.join(rows)
        return rows

    def get_domain(self):
        """Using the first column of the array
        it determines the domain assuming entry is
        univariable and the desired value of the cell
        to to be greater than zero"""

        return [(sym.solve(row[0] > 0)) for row in self.array]


def RHtable(CE, var):
    """Generates a routh horowitz table
    that can be called as a functor"""

    CE = sym.simplify(CE)
    CE = sym.expand(CE)
    CE = sym.simplify(CE)
    CE = sym.collect(CE, var)

    rows = [[] for i in range(sym.degree(CE, var)+1)]  # initialize vector
    for count, n in enumerate(range(sym.degree(CE, var), -1, -1)):
        if not count % 2:  # it's an odd number append to first row
            rows[0].append(CE.coeff(var, n))
        else:
            rows[1].append(CE.coeff(var, n))

    # append additional zeros
    while len(rows[0]) != len(rows[1]):
        rows[1].append(0)
    rows[0].append(0)
    rows[1].append(0)

    # fill up the rest of the rows
    for row, r in enumerate(rows[2:], start=2):  # start at 2nd row
        x = []
        for col in range(len(rows[row-1]) - 1):
            try:
                val = ((rows[row-1][0] * rows[row-2][col+1]) -
                (rows[row-2][col] * rows[row-1][col+1]))/rows[row-1][0]
            except ZeroDivisionError:
                val = sym.nan
            finally:
                x.append(val)
        rows[row].extend(x)

    return rows


def gain_for_zeta(gain_var, zeta, ce, start, stop, num=None):
    """Takes a characteristic equation and solves for the gain
    that produces the desired value of zeta.

    args: gain_var, variable to sub in the gain for
        zeta: the desired damping ratio
        ce: symbolic characteristic equation with only the gain_var and
            what ever the polynomial char is
        start: starting gain value
        stop: ending gain value
        num: number of equally spaced gain values
    returns:
        angle delta: how many degrees the closest value is
                     from
    """
    import warnings

    if num is None:
        num = (stop - start) + 1

    target = 180 - np.rad2deg(np.arccos(zeta))
    gains = np.linspace(start, stop, num=num)
    results = []
    for gain in gains:
        poly, *_ = sym2poly(ce.subs(gain_var, gain))
        poly = np.array(poly, dtype=np.float64)
        roots = np.roots(poly)
        angles = np.angle(roots, deg=True)
        results.append((gain, angles))

    closest = [(g, min(abs(a - target))) for g, a in results]
    closest = min(closest,  key=lambda x: x[1])

    # if angle delta is greater than 1 degree
    if closest[1] > 1:
        warnings.warn("answer has significant error, targen may be asymptotic")

    return closest
