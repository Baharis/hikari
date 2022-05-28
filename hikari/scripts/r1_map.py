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
    kwargs = locals()
    ape = angular_property_explorer_factory.create(prop='r1')
    ape.set_up(**kwargs)
    ape.explore()
    ape.write_hist_file()
    ape.draw_matplotlib_map()
    ape.draw_gnuplot_map()


if __name__ == '__main__':
    # r1_map(5.64109, 5.64109, 5.64109, 90, 90, 90, space_group='Fm-3m',
    #       job_directory='~/_/NaCl', job_name='NaCl',
    #       output_quality=5, wavelength='MoKa')
    # r1_map(5.64109, 5.64109, 5.64109, 90, 90, 90, space_group='Fm-3m',
    #        job_directory='~/Documents/python_stubs/r1_map_tests2',
    #        job_name='NaCl', output_quality=2, wavelength='MoKa')
    r1_map(a=5.641087, b=5.641087, c=5.641087, al=90, be=90, ga=90,
           space_group='Fm-3m', fix_scale=False,
           path='~/Documents/python_stubs/r1_map_tests2/NaCl.hkl',
           output_quality=2)
    pass
