def moll(x, y, field):
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    from math import pi

    rad = 180.0 / pi

    plt.figure(figsize=(8, 4))
    cmbmap = Basemap(projection='moll', lon_0=0, resolution='l')
    cmbmap.contourf(x * rad, y * rad, field, 512, cmap=plt.cm.jet, latlon=True)

    return cmbmap


def ortho(x, y, field):
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    from math import pi

    rad = 180.0 / pi

    plt.figure(figsize=(8, 4))
    cmbmap = Basemap(projection='ortho', lat_0=45, lon_0=0, resolution='l')
    cmbmap.contourf(x * rad, y * rad, field, 512, cmap=plt.cm.jet, latlon=True)

    return cmbmap


def flat(x, y, field):
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    from math import pi

    rad = 180.0 / pi

    plt.figure(figsize=(8, 4))
    cmbmap = Basemap(lon_0=0, resolution='l')
    cmbmap.contourf(x * rad, y * rad, field, 512, cmap=plt.cm.jet, latlon=True)

    return cmbmap


def point(in_cmbmap, phi, theta, size, marker):
    from math import pi

    rad = 180.0 / pi

    lons, lats = in_cmbmap(rad * phi, rad * theta)
    in_cmbmap.scatter(lons, lats, size, marker=marker)


def contour(in_cmbmap, field, x, y, number):
    from math import pi

    rad = 180.0 / pi

    in_cmbmap.contour(x * rad, y * rad, field, number, colors='k', latlon=True)


def level_plot(in_cmbmap, field, x, y, level):
    from math import pi

    rad = 180.0 / pi

    in_cmbmap.contour(x * rad, y * rad, field - level, 1, colors='k', latlon=True)


def polarization(size, in_cmbmap, q, u, x, y):
    from numpy import arctan, reciprocal
    from math import pi, fabs, cos, sin

    rad = 180.0 / pi

    phi_map = 0.5 * arctan(u * reciprocal(q))
    phi_normal_map = 0.5 * arctan(u * reciprocal(q))

    for j in xrange(1, size / 2):

        h_phi = fabs(x[size / 4][j] - x[size / 4 + 1][j])

        for i in xrange(1, size):

            if (q[i][j] >= 0) and (u[i][j] >= 0):
                phi_normal_map[i][j] = phi_map[i][j]
            elif (q[i][j] <= 0) and (u[i][j] >= 0):
                phi_normal_map[i][j] = phi_map[i][j] + pi / 2
            elif (q[i][j] >= 0) and (u[i][j] <= 0):
                phi_normal_map[i][j] = phi_map[i][j]
            elif (q[i][j] <= 0) and (u[i][j] <= 0):
                phi_normal_map[i][j] = phi_map[i][j] - pi / 2

            in_cmbmap.drawgreatcircle(rad*x[i][j], rad*y[i][j], rad*(x[i][j] + cos(phi_normal_map[i][j])*0.2*h_phi),
                                  rad*(y[i][j]+sin(phi_normal_map[i][j])*0.2*h_phi), linewidth=1, color='k')
            in_cmbmap.drawgreatcircle(rad * x[i][j], rad * y[i][j], rad*(x[i][j] - cos(phi_normal_map[i][j])*0.2*h_phi),
                                  rad*(y[i][j]-sin(phi_normal_map[i][j])*0.2*h_phi), linewidth=1, color='k')


def show():
    import matplotlib.pyplot as plt

    plt.show()
