def s2(phi1, phi2, theta1, theta2):
    from math import acos, sin, cos

    return acos(sin(theta1) * sin(theta2) + cos(theta1) * cos(theta2) * cos(phi1 - phi2))


def r2(x1, x2, y1, y2):
    from math import sqrt
    return sqrt(x1 * x2 + y1 * y2)


def cross(phi1a, theta1a, phi1b, theta1b, phi2a, theta2a, phi2b, theta2b):
    from math import sin, cos, asin, atan, sqrt

    a1 = [cos(theta1a) * cos(phi1a), cos(theta1a) * sin(phi1a), sin(theta1a)]
    b1 = [cos(theta1b) * cos(phi1b), cos(theta1b) * sin(phi1b), sin(theta1b)]
    c1 = [a1[1] * b1[2] - a1[2] * b1[1], a1[2] * b1[0] - a1[0] * b1[2], a1[0] * b1[1] - a1[1] * b1[0]]

    a2 = [cos(theta2a) * cos(phi2a), cos(theta2a) * sin(phi2a), sin(theta2a)]
    b2 = [cos(theta2b) * cos(phi2b), cos(theta2b) * sin(phi2b), sin(theta2b)]
    c2 = [a2[1] * b2[2] - a2[2] * b2[1], a2[2] * b2[0] - a2[0] * b2[2], a2[0] * b2[1] - a2[1] * b2[0]]

    a = [c1[1] * c2[2] - c1[2] * c2[1], c1[2] * c2[0] - c1[0] * c2[2], c1[0] * c2[1] - c1[1] * c2[0]]
    mod = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])

    return atan(a[1] / a[0]), asin(a[2] / mod)


def restore_value_11(phi, theta, f, x, y, i, j):

    return f[i][j] + (f[i + 1][j] - f[i][j]) * s2(x[i][j], phi, theta, theta) / \
                        s2(x[i][j], x[i + 1][j], theta, theta) + (f[i][j + 1] - f[i][j]) * \
                        s2(phi, phi, y[i][j], theta) / s2(phi, phi, y[i][j], y[i][j + 1])


def restore_value_12(phi, theta, f, x, y, i, j):

    return f[i][j + 1] + (f[i + 1][j + 1] - f[i][j + 1]) * s2(x[i][j + 1], phi, theta, theta) / \
                        s2(x[i][j + 1], x[i + 1][j + 1], theta, theta) + (f[i][j] - f[i][j + 1]) * \
                        s2(phi, phi, y[i][j + 1], theta) / s2(phi, phi, y[i][j + 1], y[i][j])


def restore_value_21(phi, theta, f, x, y, i, j):

    return f[i + 1][j] + (f[i][j] - f[i + 1][j]) * s2(x[i + 1][j], phi, theta, theta) / \
                        s2(x[i + 1][j], x[i][j], theta, theta) + (f[i + 1][j + 1] - f[i + 1][j]) * \
                        s2(phi, phi, y[i + 1][j], theta) / s2(phi, phi, y[i + 1][j], y[i + 1][j + 1])


def restore_value_22(phi, theta, f, x, y, i, j):

    return f[i + 1][j + 1] + (f[i][j + 1] - f[i + 1][j + 1]) * s2(x[i + 1][j + 1], phi, theta, theta) / \
                        s2(x[i + 1][j + 1], x[i][j + 1], theta, theta) + (f[i + 1][j] - f[i + 1][j + 1]) * \
                        s2(phi, phi, y[i + 1][j + 1], theta) / s2(phi, phi, y[i + 1][j + 1], y[i + 1][j])


def restore_value_3(phi, theta, f, x, y, i, j):

    return f[i][j] + ((f[i + 1][j] - f[i][j]) + (f[i + 1][j + 1] - f[i][j + 1])) * s2(x[i][j], phi, theta, theta) / \
                        s2(x[i][j], x[i + 1][j], theta, theta) + ((f[i][j + 1] - f[i][j]) +
                        (f[i + 1][j + 1] - f[i + 1][j])) * s2(phi, phi, y[i][j], theta) / \
                        s2(phi, phi, y[i][j], y[i][j + 1])
