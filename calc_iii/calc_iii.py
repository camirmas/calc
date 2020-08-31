from sympy import *
from sympy import (
    symbols, init_printing, nsimplify, solve, S, pprint
)


x, y, z = symbols('x y z', real=True)


def get_tangent_plane(eq, point):
    subs = {x: point[0], y: point[1]}

    f_z = point[2]
    f_x = nsimplify(eq.diff(x).evalf(subs=subs))
    f_y = nsimplify(eq.diff(y).evalf(subs=subs))

    expr = nsimplify(f_x*(x-point[0])+f_y*(y-point[1])+f_z)

    return expr


def get_differential(eq, point, approx_point):
    subs = {x: point[0], y: point[1]}
    f_x = nsimplify(eq.diff(x).evalf(subs=subs))
    f_y = nsimplify(eq.diff(y).evalf(subs=subs))

    subs_approx = {x: approx_point[0], y: approx_point[1]}

    expr = (f_x*(x-point[0])+f_y*(y-point[1])).evalf(subs=subs_approx)

    return expr


def find_critical_points(eq, verbose=True):
    if verbose:
        pprint(eq)
    f_x = nsimplify(eq.diff(x))
    f_y = nsimplify(eq.diff(y))

    result = solve([f_x, f_y], [x, y], domain=S.Reals)

    if isinstance(result, dict):
        result = [(result[x], result[y])]

    return result


def second_derivative_test(eq):
    pprint(eq)
    f_xx = nsimplify(eq.diff(x, x))
    f_yy = nsimplify(eq.diff(y, y))
    f_xy = nsimplify(eq.diff(x, y))

    critical_points = find_critical_points(eq, False)
    pprint(critical_points)

    for point in critical_points:
        subs = {x: point[0], y: point[1]}
        determinant = (f_xx*f_yy-f_xy**2).evalf(subs=subs)

        if determinant < 0:
            print("(%s, %s) is a saddle point" % point)
        elif determinant > 0:
            if f_xx.evalf(subs=subs) < 0:
                print("(%s, %s) is a local maximum" % point)
            else:
                print("(%s, %s) is a local minimum" % point)
        else:
            print("(%s, %s) is inconclusive" % point)


if __name__ == "__main__":
    init_printing(use_unicode=True)
