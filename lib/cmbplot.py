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


def show():
    import matplotlib.pyplot as plt

    plt.show()
