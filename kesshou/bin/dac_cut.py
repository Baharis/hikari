from kesshou.dataframes.cif import CifFrame
from kesshou.dataframes.hkl import HklFrame
import numpy as np

p = HklFrame()
c = CifFrame('/home/dtchon/git/kesshou/test_data/exp_353.cif',
             file_data_block='exp_353')
p.read('/home/dtchon/git/kesshou/test_data/exp_353_neu.hkl', 2)
p.crystal.edit_cell(a=11.0, b=12.0, c=16.0)

p.crystal.orient_matrix = np.array(((0, -1, 0), (1, 0, 0), (0, 0, 1)))
p.crystal.import_from_frame(c)
p.drop_zero()
p.reduce()
p.place()
# p.calculate_uncertainty('I')
# p.dac(opening_angle=35)
# p.write('/home/dtchon/git/kesshou/test_data/output.hkl', 4)
# p.draw(projection=(0, 'k', 'l'), scale=1,
#        savepath='/home/dtchon/git/kesshou/test_data/output.png',
#        showfig=True, colored='r', alpha='u')
print(p.data['r'].max()/2)
