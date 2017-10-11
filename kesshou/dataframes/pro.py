import os
import pandas as pd
from collections import OrderedDict


class ProFrame:
    def __init__(self, file_path=None):
        self.data = pd.DataFrame()
        self.meta = dict()
        if file_path is not None:
            self.read(path=file_path)

    def read(self, path):
        """Read data from specified xd_pro.out file and return an DataFrame"""

        # SPECIFY SOME META
        self.meta['path'] = os.path.abspath(path)
        self.meta['name'] = str()
        self.meta['comment'] = str()

        # DECLARE SOUGHT .DATA FOR PANDAS DATAFRAME
        pro_keys_list = (
            # First line information :
            'ID', 'Label', 'iZ',
            # Atomic populations (Electrons) :
            'N', 'Q',
            # Atomic energies (Atomic Units) :
            'L', 'VNEO', 'VNET',
            # Atomic forces (Atomic Units) :
            'FAXA', 'FAYA', 'FAZA', 'FBXA', 'FBYA', 'FBZA',
            # Radial atomic expectation values (Atomic Units) :
            'R(-1)', 'R(+1)', 'R(+2)', 'R(+3)', 'R(+4)',
            'GR(-1)', 'GR(', 'GR(+1)', 'GR(+2)',
            # Atomic dipole moment vector (Debye) :
            'DX', 'DY', 'DZ', 'DM',
            # Atomic displacement vector (Angstroms) :
            'DCX', 'DCY', 'DCZ',
            # Coordinates of the centroid of negative charge (Angstroms) :
            'IUN', 'INO',
            # Atomic volumes and related properties (Angstroms^3) :
            'V001', 'N001', 'R001', 'V002', 'N002', 'R002', 'VTOT')
        pro_keys_dict, pro_data_dict = OrderedDict(), OrderedDict()
        for key in pro_keys_list:
            pro_keys_dict[key] = None
            pro_data_dict[key] = []
        pro_atom_dict = pro_keys_dict

        # READ THE FILE AND LOAD THE DATA
        file = open(path, "r", encoding='utf-8', errors='ignore')
        for line in file:
            # IGNORE EMPTY LINES
            if not line.strip():
                pass
            # ADD COMMENTS TO META
            elif line[0] == '!':
                self.meta['comment'] += line
            # FOR NEW ATOM CREATE NEW PRO_ATOM_DICT AND PUSH OLD IF EXISTS
            elif 'Atom ' in line.strip() and 'iZ =' in line.strip():
                if pro_atom_dict['ID'] is not None:
                    for k, v in pro_atom_dict.items():
                        if not v:
                            v = 0.0
                            pro_atom_dict[k] = v
                    for k, v in pro_atom_dict.items():
                        pro_data_dict[k].append(v)
                pro_atom_dict = pro_keys_dict
                pro_atom_dict['ID']     = int(line.strip().split()[1])
                pro_atom_dict['Label']  = line.strip().split()[3]
                pro_atom_dict['iZ']     = int(line.strip().split()[6])
            # ADD EXCEPTED VALUES TO THE ATOM DICTIONARY
            elif pro_atom_dict['ID'] and line.strip().split()[0] in pro_keys_list:
                if line.strip().split()[0] == 'GR(':
                    pro_atom_dict['GR('] = float(line.strip().split()[2])
                    continue
                try:
                    k, v = line.strip().split()[0:2]
                except IndexError:
                    pass
                else:
                    pro_atom_dict[k] = float(v)
        # ON THE END OF THE FILE PUSH LAST ATOM TO PRO_DATA_DICT
        for k, v in pro_atom_dict.items():
            if not v:
                v = 0.0
                pro_data_dict[k].append(v)
        for k, v in pro_atom_dict.items():
            pro_data_dict[k].append(v)
        file.close()

        # PREPARE LIST OF "@Other calculated information": '@N-Rsc' and '@Q-Rsc'
        pro_N_sum, pro_Z_sum = 0.0, 0.0
        for N, Z in zip(pro_data_dict['N'], pro_data_dict['iZ']):
            pro_N_sum += N
            pro_Z_sum += Z
        pro_data_dict['@N-Rsc'], pro_data_dict['@Q-Rsc'] = [], []
        for N, Z in zip(pro_data_dict['N'], pro_data_dict['iZ']):
            pro_data_dict['@N-Rsc'].append(N * pro_Z_sum / pro_N_sum)
            pro_data_dict['@Q-Rsc'].append(N * pro_Z_sum / pro_N_sum - Z)

        # SAVE .DATA FROM PRO_DATA_DICT
        self.data = pd.DataFrame.from_dict(pro_data_dict)

    def write(self, path):
        self.data.to_csv(path, decimal='.')


if __name__ == '__main__':
    pro = ProFrame(file_path='/home/dtchon/git/kesshou/test_data/xd_pro.out')
    pro.write('/home/dtchon/git/kesshou/test_data/xd_pro_digested.out')
