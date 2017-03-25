def t_solver(q1, q2, u1, u2):
    import numpy as np
    return np.roots([u2, u1 + 2 * q2, 2 * q1 - u2, -u1])


def mean_t_solver(q1, q2, u1, u2, x):
    return u2 * x**3 + (u1 + 2 * q2) * x**2 + (2 * q1 - u2) * x + (-u1)


def d_solver(q1, q2, u1, u2):
    import numpy as np
    return np.roots([3 * u2, 2 * (u1 + 2 * q2), (2 * q1 - u2)])


def det_d_solver(q1, q2, u1, u2):
    return (2 * (u1 + 2 * q2))*(2 * (u1 + 2 * q2)) - 4 * 3 * u2 * (2 * q1 - u2)


def det_func(q1, q2, u1, u2):
    return q1 * u2 - q2 * u1


def type_solver(qx, qy, ux, uy, param=False):

    if not param:
        import numpy as np

        root1_d, root2_d = d_solver(qx, qy, ux, uy)
        mean_root1_d = mean_t_solver(qx, qy, ux, uy, root1_d)
        mean_root2_d = mean_t_solver(qx, qy, ux, uy, root2_d)

        if det_func(qx, qy, ux, uy) < 0:
            return 1
        elif (mean_root1_d > 0 > mean_root2_d) or (mean_root1_d < 0 < mean_root2_d):
            return 2
        else:
            return 3

    else:
        import numpy as np
        parameter = 1e-5

        root1, root2, root3 = t_solver(qx, qy, ux, uy)

        if det_func(qx, qy, ux, uy) < 0:
            return 1
        elif (-parameter < np.imag(root1) < parameter and
                      -parameter < np.imag(root2) < parameter and
                      -parameter < np.imag(root3) < parameter):
            return 2
        else:
            return 3


def types_search(f, qx, qy, ux, uy, whitelist, up_bounds=False, down_bounds=False):

    n = qx.shape[0] - 1

    for i in xrange(0, n):
        for j in xrange(1, n / 2):

            if whitelist[i][j] != 0:

                if down_bounds == False and up_bounds == False:
                    return type_solver(qx[i][j], qy[i][j], ux[i][j], uy[i][j])
                elif down_bounds != False or up_bounds != False and \
                        ((down_bounds <= f[i][j] <= up_bounds) or
                             (down_bounds >= f[i][j] >= up_bounds)):
                    return type_solver(qx[i][j], qy[i][j], ux[i][j], uy[i][j])
