from kesshou.dataframes.hkl import HklFrame
from kesshou.utility import fibonacci_sphere
from kesshou.symmetry import PG


def dac_cplt_descriptive_statistics(a, b, c, al, be, ga,
                                    laue_group=PG['-1'],
                                    output_path='output.txt',
                                    opening_angle=35,
                                    precision=1000,
                                    random_seed=1337,
                                    resolution=None,
                                    wavelength='MoKa'):
    """Calculate max, min, avg completeness and its other statistics
    for a series of multiple random orientations of diamond anvil cell"""

    def make_reference_ball():
        hkl_frame = HklFrame()
        hkl_frame.crystal.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
        hkl_frame.edit_wavelength(wavelength)
        hkl_frame.make_ball(radius=hkl_frame.r_lim)
        hkl_frame.merge()
        hkl_frame.extinct('000')
        if not(resolution is None):
            hkl_frame.trim(resolution)
        return hkl_frame
    p = make_reference_ball()
    total_reflections = len(p)
    max_resolution = max(p.data['r'])
    out = open(output_path, 'w', buffering=1)
    out.write('total_reflections: ' + str(total_reflections) + '\n')
    out.write('maximum_r_in_reciprocal_coordinates: ' + str(max_resolution))
    vectors = fibonacci_sphere(samples=precision, seed=random_seed)
    reflections = list()
    for vector in vectors:
        q = p.duplicate()
        q.dac(opening_angle=opening_angle, vector=vector)
        q.resymmetrify(laue_group.chiral_operations, merge=True)
        reflections.append(len(q))
        out.write('\n' + str(vector) + ': ' + str(len(q)))
    out.write('\nmax_reflections: ' + str(max(reflections)))
    out.write('\nmin_reflections: ' + str(min(reflections)))
    out.write('\navg_reflections: ' + str(sum(reflections) / len(reflections)))
    out.close()


if __name__ == '__main__':
    dac_cplt_descriptive_statistics(a=10.0, b=10.0, c=10.0,
                                    al=90.0, be=90.0, ga=90.0,
                                    laue_group=PG['m-3'])
