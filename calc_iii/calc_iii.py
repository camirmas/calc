from sympy import *
from sympy.vector import matrix_to_vector, CoordSys3D
from sympy import (
    symbols, init_printing, nsimplify, solve, S, pprint, pi,
    Matrix, integrate, sqrt, sympify
)


x, y, z, t = symbols('x y z t', real=True)
C = CoordSys3D('C')


def project(u, v):
    """Projects vector u onto vector v."""
    u = matrix_to_vector(Matrix(u), C)
    v = matrix_to_vector(Matrix(v), C)
    m_v = v.magnitude()

    scal = u.dot(v)/m_v
    pprint(scal)

    return scal*(v/m_v)


def dot(v1, v2):
    return Matrix(v1).dot(Matrix(v2))


def cross(u, v):
    u = matrix_to_vector(Matrix(u), C)
    v = matrix_to_vector(Matrix(v), C)

    return u.cross(v)


def magnitude(v):
    v = matrix_to_vector(Matrix(v), C)

    return v.magnitude()


def unit(v):
    v = matrix_to_vector(Matrix(v), C)

    return v.normalize()


def plane(v, point):
    (x_p, y_p, z_p) = point

    return nsimplify(v[0]*(x - x_p) + v[1]*(y - y_p) + v[2]*(z - z_p))


def arc_length(vec, start, end):
    pprint(vec)
    derivs = 0

    for i in vec:
        derivs += i.diff(t)**2

    integrand = sqrt(derivs)
    result = integrate(integrand, (t, start, end))

    return nsimplify(result, [pi])


def gradient(eq, point=None, verbose=True):
    f_x = eq.diff(x)
    f_y = eq.diff(y)
    f_z = eq.diff(z)

    if verbose:
        pprint(eq)
        pprint((f_x, f_y, f_z))

    if point:
        subs = {x: point[0], y: point[1], z: point[2]}
        f_x = f_x.subs(subs)
        f_y = f_y.subs(subs)
        f_z = f_z.subs(subs)

    return (nsimplify(f_x, [pi]), nsimplify(f_y, [pi]), nsimplify(f_z, [pi]))


def directional_derivative(eq, point, direction):
    grad = gradient(eq, point, False)
    m1 = Matrix(grad)
    m2 = Matrix(direction)

    return nsimplify(m1.dot(m2), [pi])


def tangent_plane(eq, point, implicit=False):
    pprint(eq)
    subs = {x: point[0], y: point[1], z: point[2]}

    if implicit:
        f_z = eq.diff(z).subs(subs)
    else:
        f_z = point[2]

    f_x = eq.diff(x).subs(subs)
    f_y = eq.diff(y).subs(subs)

    if implicit:
        expr = f_x*(x-point[0])+f_y*(y-point[1])+f_z*(z-point[2])
    else:
        expr = f_x*(x-point[0])+f_y*(y-point[1])+f_z

    return nsimplify(expr, [pi])


def differential(eq, point, approx_point):
    subs = {x: point[0], y: point[1]}
    f_x = eq.diff(x).subs(subs)
    f_y = eq.diff(y).subs(subs)

    subs_approx = {x: approx_point[0], y: approx_point[1]}

    expr = (f_x*(x-point[0])+f_y*(y-point[1])).evalf(subs=subs_approx)

    return expr


def linear_approx(eq, point, approx_point):
    subs = {x: point[0], y: point[1]}
    point = (point[0], point[1], nsimplify(eq.evalf(subs=subs), [pi]))
    plane = tangent_plane(eq, point)
    pprint(plane)

    approx_subs = {x: approx_point[0], y: approx_point[1]}

    return nsimplify(plane.evalf(subs=approx_subs), [pi])


def critical_points(eq, verbose=True):
    if verbose:
        pprint(eq)
    f_x = eq.diff(x)
    f_y = eq.diff(y)

    result = solve([f_x, f_y], [x, y], domain=S.Reals)

    if isinstance(result, dict):
        result = [(result[x], result[y])]

    return result


def second_derivative_test(eq):
    pprint(eq)
    f_xx = eq.diff(x, x)
    f_yy = eq.diff(y, y)
    f_xy = eq.diff(x, y)

    crit_points = critical_points(eq, False)
    pprint(crit_points)

    for point in crit_points:
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
