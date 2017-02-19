def area(y, field, level=0.0):
    from math import cos

    a = 0.0  # Area without normalization
    na = 0.0  # Normalization

    n = field.shape[0] - 1

    for i in xrange(0, n):
        for j in xrange(1, n / 2):

            mean = (field[i][j] + field[i + 1][j + 1] + field[i + 1][j] + field[i][j + 1]) / 4.0

            if mean >= level:
                a += cos(y[i][j])
            na += cos(y[i][j])

    return a / na


def length(x, y, field, level=0.0):
    from math import fabs, pi
    from distance_1 import s2

    n = field.shape[0] - 1

    l = 0.0

    f = field - level

    for i in xrange(0, n):
        for j in xrange(1, n / 2):

            h_theta = fabs(y[n / 2 + 1][n / 4])
            h_phi = fabs(x[n / 4][j] - x[n / 4 + 1][j])

            if f[i][j] == 0.0 and (f[i][j + 1] == 0.0 or f[i + 1][j] == 0.0):
                l += 0.0

            if field[i][j] * field[i][j + 1] < 0.0:

                if f[i][j] * f[i + 1][j] < 0.0:

                    phi1 = x[i][j]
                    theta1 = y[i][j] + h_theta * fabs(f[i][j]) / (fabs(f[i][j]) + fabs(f[i][j + 1]))

                    phi2 = x[i][j] + h_phi * fabs(f[i][j]) / (fabs(f[i][j]) + fabs(f[i + 1][j]))
                    theta2 = y[i][j]

                    l += s2(phi1, phi2, theta1, theta2)

                elif f[i + 1][j] * f[i + 1][j + 1] < 0.0:

                    phi1 = x[i][j]
                    theta1 = y[i][j] + h_theta * fabs(f[i][j]) / (fabs(f[i][j]) + fabs(f[i][j + 1]))

                    phi2 = x[i + 1][j]
                    theta2 = y[i + 1][j] + h_theta * fabs(f[i + 1][j]) / (fabs(f[i + 1][j]) + fabs(f[i + 1][j + 1]))

                    l += s2(phi1, phi2, theta1, theta2)

                elif f[i][j + 1] * f[i + 1][j + 1] < 0.0:

                    phi1 = x[i][j]
                    theta1 = y[i][j] + h_theta * fabs(f[i][j]) / (fabs(f[i][j]) + fabs(f[i][j + 1]))

                    phi2 = x[i][j + 1] + h_phi * fabs(f[i][j + 1]) / (fabs(f[i][j + 1]) + fabs(f[i + 1][j + 1]))
                    theta2 = y[i][j + 1]

                    l += s2(phi1, phi2, theta1, theta2)

            elif f[i][j] * f[i + 1][j] <= 0.0:

                if f[i + 1][j] * f[i + 1][j + 1] < 0.0:

                    phi1 = x[i][j] + h_phi * fabs(f[i][j]) / (fabs(f[i][j]) + fabs(f[i + 1][j]))
                    theta1 = y[i][j]

                    phi2 = x[i + 1][j]
                    theta2 = y[i + 1][j] + h_theta * fabs(f[i + 1][j]) / (fabs(f[i + 1][j]) + fabs(f[i + 1][j + 1]))

                    l += s2(phi1, phi2, theta1, theta2)

                elif f[i][j + 1] * f[i + 1][j + 1] < 0.0:

                    phi1 = x[i][j] + h_phi * fabs(f[i][j]) / (fabs(f[i][j]) + fabs(f[i + 1][j]))
                    theta1 = y[i][j]

                    phi2 = x[i][j + 1] + h_phi * fabs(f[i][j + 1]) / (fabs(f[i][j + 1]) + fabs(f[i + 1][j + 1]))
                    theta2 = y[i][j + 1]

                    l += s2(phi1, phi2, theta1, theta2)

            elif f[i + 1][j] * f[i + 1][j + 1] < 0.0:

                if f[i][j + 1] * f[i + 1][j + 1] < 0.0:
                    phi1 = x[i + 1][j]
                    theta1 = y[i + 1][j] + h_theta * fabs(f[i + 1][j]) / (fabs(f[i + 1][j]) + fabs(f[i + 1][j + 1]))

                    phi2 = x[i][j + 1] + h_phi * fabs(f[i][j + 1]) / (fabs(f[i][j + 1]) + fabs(f[i + 1][j + 1]))
                    theta2 = y[i][j + 1]

                    l += s2(phi1, phi2, theta1, theta2)

    return l / (4 * pi)


def condition_1(xx, yy, xy):
    if (xx * yy - xy * xy >= 0.0 and xx >= 0.0) or (xx * yy - xy * xy >= 0.0 and yy > 0.0):
        return 0
    else:
        if (xx * yy - xy * xy >= 0.0 >= xx) or (xx * yy - xy * xy >= 0.0 >= yy):
            return 2
        else:
            return 1


def type_points(x, y, f, fx, fy, fxx, fyy, fxy, sigma_0, sigma_1, sigma_2, my_file=False, my_map=False,
                up_bounds=False, down_bounds=False, whitelist=False, whitelist_flag=False):
    from numpy import zeros
    from math import fabs
    from lib.distance_1 import s2, cross

    n = f.shape[0] - 1

    g = 0.0
    n_min = 0.0
    n_max = 0.0
    n_sad = 0.0

    z_x = zeros((n + 1, n / 2 + 1))
    z_y = zeros((n + 1, n / 2 + 1))

    phi1a = 0.0
    phi1b = 0.0

    phi2a = 0.0
    phi2b = 0.0

    theta1a = 0.0
    theta1b = 0.0

    theta2a = 0.0
    theta2b = 0.0

    for i in xrange(0, n):
        for j in xrange(1, n / 2):

            if ((whitelist_flag != False) and (whitelist[i][j] == 0)) or (whitelist_flag == False):

                h_theta = fabs(y[n / 2 + 1][n / 4])
                h_phi = fabs(x[n / 4][j] - x[n / 4 + 1][j])

                if fx[i][j] * fx[i][j + 1] < 0.0:

                    if fx[i][j] * fx[i + 1][j] < 0.0:

                        phi1a = x[i][j]
                        theta1a = y[i][j] + h_theta * fabs(fx[i][j]) / (fabs(fx[i][j]) + fabs(fx[i][j + 1]))

                        phi1b = x[i][j] + h_phi * fabs(fx[i][j]) / (fabs(fx[i][j]) + fabs(fx[i + 1][j]))
                        theta1b = y[i][j]

                        z_x[i][j] = 1

                    elif fx[i + 1][j] * fx[i + 1][j + 1] < 0.0:

                        phi1a = x[i][j]
                        theta1a = y[i][j] + h_theta * fabs(fx[i][j]) / (fabs(fx[i][j]) + fabs(fx[i][j + 1]))

                        phi1b = x[i + 1][j]
                        theta1b = y[i + 1][j] + h_theta * fabs(fx[i + 1][j]) / (fabs(fx[i + 1][j])
                                                                                + fabs(fx[i + 1][j + 1]))

                        z_x[i][j] = 1

                    elif fx[i][j + 1] * fx[i + 1][j + 1] < 0.0:

                        phi1a = x[i][j]
                        theta1a = y[i][j] + h_theta * fabs(fx[i][j]) / (fabs(fx[i][j]) + fabs(fx[i][j + 1]))

                        phi1b = x[i][j + 1] + h_phi * fabs(fx[i][j + 1]) / (fabs(fx[i][j + 1]) + fabs(fx[i + 1][j + 1]))
                        theta1b = y[i][j + 1]

                        z_x[i][j] = 1

                elif fx[i][j] * fx[i + 1][j] < 0.0:

                    if fx[i + 1][j] * fx[i + 1][j + 1] < 0.0:

                        phi1a = x[i][j] + h_phi * fabs(fx[i][j]) / (fabs(fx[i][j]) + fabs(fx[i + 1][j]))
                        theta1a = y[i][j]

                        phi1b = x[i + 1][j]
                        theta1b = y[i + 1][j] + h_theta * fabs(fx[i + 1][j]) / (
                            fabs(fx[i + 1][j]) + fabs(fx[i + 1][j + 1]))

                        z_x[i][j] = 1

                    elif fx[i][j + 1] * fx[i + 1][j + 1] < 0.0:

                        phi1a = x[i][j] + h_phi * fabs(fx[i][j]) / (fabs(fx[i][j]) + fabs(fx[i + 1][j]))
                        theta1a = y[i][j]

                        phi1b = x[i][j + 1] + h_phi * fabs(fx[i][j + 1]) / (fabs(fx[i][j + 1]) + fabs(fx[i + 1][j + 1]))
                        theta1b = y[i][j + 1]

                        z_x[i][j] = 1

                elif fx[i + 1][j] * fx[i + 1][j + 1] < 0.0:

                    if fx[i][j + 1] * fx[i + 1][j + 1] < 0.0:

                        phi1a = x[i + 1][j]
                        theta1a = y[i + 1][j] + h_theta * fabs(fx[i + 1][j]) / (
                            fabs(fx[i + 1][j]) + fabs(fx[i + 1][j + 1]))

                        phi1b = x[i][j + 1] + h_phi * fabs(fx[i][j + 1]) / (fabs(fx[i][j + 1]) + fabs(fx[i + 1][j + 1]))
                        theta1b = y[i][j + 1]

                        z_x[i][j] = 1

                if fy[i][j] * fy[i][j + 1] < 0.0:

                    if fy[i][j] * fy[i + 1][j] < 0.0:

                        phi2a = x[i][j]
                        theta2a = y[i][j] + h_theta * fabs(fy[i][j]) / (fabs(fy[i][j]) + fabs(fy[i][j + 1]))

                        phi2b = x[i][j] + h_phi * fabs(fy[i][j]) / (fabs(fy[i][j]) + fabs(fy[i + 1][j]))
                        theta2b = y[i][j]

                        z_y[i][j] = 1

                    elif fy[i + 1][j] * fy[i + 1][j + 1] < 0.0:

                        phi2a = x[i][j]
                        theta2a = y[i][j] + h_theta * fabs(fy[i][j]) / (fabs(fy[i][j]) + fabs(fy[i][j + 1]))

                        phi2b = x[i + 1][j]
                        theta2b = y[i + 1][j] + h_theta * fabs(fy[i + 1][j]) / (
                            fabs(fy[i + 1][j]) + fabs(fy[i + 1][j + 1]))

                        z_y[i][j] = 1

                    elif fy[i][j + 1] * fy[i + 1][j + 1] < 0.0:

                        phi2a = x[i][j]
                        theta2a = y[i][j] + h_theta * fabs(fy[i][j]) / (fabs(fy[i][j]) + fabs(fy[i][j + 1]))

                        phi2b = x[i][j + 1] + h_phi * fabs(fy[i][j + 1]) / (fabs(fy[i][j + 1]) + fabs(fy[i + 1][j + 1]))
                        theta2b = y[i][j + 1]

                        z_y[i][j] = 1

                elif fy[i][j] * fy[i + 1][j] < 0.0:

                    if fy[i + 1][j] * fy[i + 1][j + 1] < 0.0:

                        phi2a = x[i][j] + h_phi * fabs(fy[i][j]) / (fabs(fy[i][j]) + fabs(fy[i + 1][j]))
                        theta2a = y[i][j]

                        phi2b = x[i + 1][j]
                        theta2b = y[i + 1][j] + h_theta * fabs(fy[i + 1][j]) / (
                            fabs(fy[i + 1][j]) + fabs(fy[i + 1][j + 1]))

                        z_y[i][j] = 1

                    elif fy[i][j + 1] * fy[i + 1][j + 1] < 0.0:

                        phi2a = x[i][j] + h_phi * fabs(fy[i][j]) / (fabs(fy[i][j]) + fabs(fy[i + 1][j]))
                        theta2a = y[i][j]

                        phi2b = x[i][j + 1] + h_phi * fabs(fy[i][j + 1]) / (fabs(fy[i][j + 1]) + fabs(fy[i + 1][j + 1]))
                        theta2b = y[i][j + 1]

                        z_y[i][j] = 1

                elif fy[i + 1][j] * fy[i + 1][j + 1] < 0.0:

                    if fy[i][j + 1] * fy[i + 1][j + 1] < 0.0:

                        phi2a = x[i + 1][j]
                        theta2a = y[i + 1][j] + h_theta * fabs(fy[i + 1][j]) / (
                            fabs(fy[i + 1][j]) + fabs(fy[i + 1][j + 1]))

                        phi2b = x[i][j + 1] + h_phi * fabs(fy[i][j + 1]) / (fabs(fy[i][j + 1]) + fabs(fy[i + 1][j + 1]))
                        theta2b = y[i][j + 1]

                        z_y[i][j] = 1

                if (z_y[i][j] != 0 and z_x[i][j] != 0) and (down_bounds == False and up_bounds == False):

                    flag = 0

                    phi_precision = 0.0
                    theta_precision = 0.0

                    phi_a, theta_a = cross(phi1a, theta1a, phi1b, theta1b, phi2a, theta2a, phi2b, theta2b)

                    if (x[i][j] <= phi_a <= x[i + 1][j]) and (y[i][j] <= theta_a <= y[i][j + 1]):

                        phi_precision = phi_a
                        theta_precision = theta_a
                        flag = 1

                    if (x[i][j] <= phi_a <= x[i + 1][j]) and (y[i][j] <= - theta_a <= y[i][j + 1]):

                        phi_precision = phi_a
                        theta_precision = - theta_a
                        flag = 1

                    if flag == 1:

                        fxx_precision = fxx[i][j] + (fxx[i + 1][j] - fxx[i][j]) * \
                                        s2(phi_precision, x[i][j], theta_precision, theta_precision) / \
                                        s2(x[i + 1][j], x[i][j], theta_precision, theta_precision) + \
                                        (fxx[i][j + 1] - fxx[i][j]) * \
                                        s2(phi_precision, phi_precision, theta_precision, y[i][j]) / \
                                        s2(phi_precision, phi_precision, y[i][j], y[i][j + 1])

                        fyy_precision = fyy[i][j] + (fyy[i + 1][j] - fyy[i][j]) * \
                                        s2(phi_precision, x[i][j], theta_precision, theta_precision) / \
                                        s2(x[i + 1][j], x[i][j], theta_precision, theta_precision) + \
                                        (fyy[i][j + 1] - fyy[i][j]) * \
                                        s2(phi_precision, phi_precision, theta_precision, y[i][j]) / \
                                        s2(phi_precision, phi_precision, y[i][j], y[i][j + 1])

                        fxy_precision = fxy[i][j] + (fxy[i + 1][j] - fxy[i][j]) * \
                                        s2(phi_precision, x[i][j], theta_precision, theta_precision) / \
                                        s2(x[i + 1][j], x[i][j], theta_precision, theta_precision) + \
                                        (fxy[i][j + 1] - fxy[i][j]) * \
                                        s2(phi_precision, phi_precision, theta_precision, y[i][j]) / \
                                        s2(phi_precision, phi_precision, y[i][j], y[i][j + 1])

                        cond_answ = condition_1(fxx_precision, fyy_precision, fxy_precision)

                        if cond_answ == 0:

                            my_type = 'o'
                            g += 1
                            n_min += 1
                            ms = 15

                        if cond_answ == 2:

                            my_type = '+'
                            g += 1
                            n_max += 1
                            ms = 100

                        if cond_answ == 1:

                            my_type = 'x'
                            g -= 1
                            n_sad -= 1
                            ms = 100

                        if my_map:

                            from lib.cmbplot import point

                            point(my_map, phi_precision, theta_precision, ms, my_type)

                        if my_file:

                            my_file.write(repr(f[i][j]) + '    ' + repr(phi_precision) + '    ' +
                            repr(theta_precision) + '   ' + repr(cond_answ) + '    ' + repr(fxx_precision) + '    ' +
                            repr(fyy_precision) + '    ' + repr(fxy_precision) + '    ' +
                            repr(fxx_precision * fyy_precision - fxy_precision * fxy_precision) + '\n')

                            my_file.write(repr(f[i][j]) + '    ' + repr(cond_answ) + '    ' + repr(sigma_0) + '    ' +
                            repr(sigma_1) + '    ' + repr(sigma_2) + '\n')

                if ((z_y[i][j] != 0 and z_x[i][j] != 0) and (down_bounds != False or up_bounds != False)) and (
                            (down_bounds <= f[i][j] <= up_bounds) or (down_bounds >= f[i][j] >= up_bounds)):

                    flag = 0

                    phi_precision = 0.0
                    theta_precision = 0.0

                    phi_a, theta_a = cross(phi1a, theta1a, phi1b, theta1b, phi2a, theta2a, phi2b, theta2b)

                    if (x[i][j] <= phi_a <= x[i + 1][j]) and (y[i][j] <= theta_a <= y[i][j + 1]):

                        phi_precision = phi_a
                        theta_precision = theta_a
                        flag = 1

                    if (x[i][j] <= phi_a <= x[i + 1][j]) and (y[i][j] <= - theta_a <= y[i][j + 1]):

                        phi_precision = phi_a
                        theta_precision = - theta_a
                        flag = 1

                    if flag != 0:

                        fxx_precision = fxx[i][j] + (fxx[i + 1][j] - fxx[i][j]) * \
                                        s2(phi_precision, x[i][j], theta_precision, theta_precision) / \
                                        s2(x[i + 1][j], x[i][j], theta_precision, theta_precision) + \
                                        (fxx[i][j + 1] - fxx[i][j]) * \
                                        s2(phi_precision, phi_precision, theta_precision, y[i][j]) / \
                                        s2(phi_precision, phi_precision, y[i][j], y[i][j + 1])

                        fyy_precision = fyy[i][j] + (fyy[i + 1][j] - fyy[i][j]) * \
                                        s2(phi_precision, x[i][j], theta_precision, theta_precision) / \
                                        s2(x[i + 1][j], x[i][j], theta_precision, theta_precision) + \
                                        (fyy[i][j + 1] - fyy[i][j]) * \
                                        s2(phi_precision, phi_precision, theta_precision, y[i][j]) / \
                                        s2(phi_precision, phi_precision, y[i][j], y[i][j + 1])

                        fxy_precision = fxy[i][j] + (fxy[i + 1][j] - fxy[i][j]) * \
                                        s2(phi_precision, x[i][j], theta_precision, theta_precision) / \
                                        s2(x[i + 1][j], x[i][j], theta_precision, theta_precision) + \
                                        (fxy[i][j + 1] - fxy[i][j]) * \
                                        s2(phi_precision, phi_precision, theta_precision, y[i][j]) / \
                                        s2(phi_precision, phi_precision, y[i][j], y[i][j + 1])

                        cond_answ = condition_1(fxx_precision, fyy_precision, fxy_precision)

                        if cond_answ == 0:

                            my_type = 'o'
                            g += 1
                            n_min += 1
                            ms = 15

                        if cond_answ == 2:

                            my_type = '+'
                            g += 1
                            n_max += 1
                            ms = 100

                        if cond_answ == 1:

                            my_type = 'x'
                            g -= 1
                            n_sad -= 1
                            ms = 100

                        if my_map:
                            from lib.cmbplot import point

                            point(my_map, phi_precision, theta_precision, ms, my_type)

                        if my_file:
                            my_file.write(repr(f[i][j]) + '    ' + repr(phi_precision) + '    ' +
                                          repr(theta_precision) + '   ' + repr(cond_answ) + '    ' + repr(
                                fxx_precision) + '    ' +
                                          repr(fyy_precision) + '    ' + repr(fxy_precision) + '    ' +
                                          repr(fxx_precision * fyy_precision - fxy_precision * fxy_precision) + '\n')

                            my_file.write(repr(f[i][j]) + '    ' + repr(cond_answ) + '    ' + repr(sigma_0) + '    ' +
                                          repr(sigma_1) + '    ' + repr(sigma_2) + '\n')

    if down_bounds != False or up_bounds != False:
        return g, n_max, n_min, n_sad
