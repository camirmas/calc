from sympy import *
from sympy import symbols, init_printing, evalf, nsimplify


def tangent_plane(eq, point):
    subs = {x: point[0], y: point[1]}

    f_z = point[2]
    f_x = nsimplify(eq.diff(x).evalf(subs=subs))
    f_y = nsimplify(eq.diff(y).evalf(subs=subs))

    expr = nsimplify(f_x*(x-point[0])+f_y*(y-point[1])+f_z)

    return expr


def differential(eq, point, approx_point):
    subs = {x: point[0], y: point[1]}
    f_x = nsimplify(eq.diff(x).evalf(subs=subs))
    f_y = nsimplify(eq.diff(y).evalf(subs=subs))

    subs_approx = {x: approx_point[0], y: approx_point[1]}

    expr = (f_x*(x-point[0])+f_y*(y-point[1])).evalf(subs=subs_approx)

    return expr


if __name__ == "__main__":
    x, y, z = symbols('x y z')
    init_printing(use_unicode=True)
