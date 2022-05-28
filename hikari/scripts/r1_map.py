from hikari.scripts.angular_explorer import angular_property_explorer_factory


def r1_map(a, b, c, al, be, ga,
           space_group='P1',
           axis='',
           fix_scale=False,
           histogram=True,
           opening_angle=35,
           orientation=None,
           path='~/sortav.lst',
           output_quality=3,
           resolution=1.2,
           wavelength='MoKa'):
    """
    Calculate and draw a r1 map for a given crystal in diamond anvil cell
    (DAC) with a given opening angle, as a function of crystal orientation.

    The script accepts unit cell & space group information, runs SHELXL,
    and reads value of R1 after single refinement. Results are logged into text
    files and drawn with gnuplot or matplotlib, depending on settings.

    For further detail concerning r1_map, its basis and uses, refer to
    :py:func:`hikari.scripts.potency_map`, as well as selected terminology
    described in `this paper <https://doi.org/10.1107/S2052252521009532>`_.
    """
    kwargs = locals()
    ape = angular_property_explorer_factory.create(prop='r1')
    ape.set_up(**kwargs)
    ape.explore()
    ape.write_hist_file()
    ape.draw_matplotlib_map()
    ape.draw_gnuplot_map()


if __name__ == '__main__':
    r1_map(5.64109, 5.64109, 5.64109, 90, 90, 90, space_group='Fm-3m',
           path='~/_/NaCl/NaCl.hkl', output_quality=5, wavelength='MoKa')
