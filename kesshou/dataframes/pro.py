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
        pro_single_keys = (
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
            'GR(-1)', 'GR( 0)', 'GR(+1)', 'GR(+2)',
            # Atomic dipole moment vector (Debye) :
            'DX', 'DY', 'DZ', 'DM',
            # Atomic displacement vector (Angstroms) :
            'DCX', 'DCY', 'DCZ',
            # Coordinates of the centroid of negative charge (Angstroms) :
            'IUN', 'INO',
            # Atomic volumes and related properties (Angstroms^3) :
            'V001', 'N001', 'R001', 'V002', 'N002', 'R002', 'VTOT',
            # Not directly seeked, but necessary keys for nine SFs:
            '_SF1', '_SF2', '_SF3', '_SF4', '_SF5',
            '_SF6', '_SF7', '_SF8', '_SF9',
        )
        pro_keys_dict, pro_data_dict = OrderedDict(), OrderedDict()
        for key in pro_single_keys:
            pro_keys_dict[key] = None
            pro_data_dict[key] = []
        pro_atom_dict = pro_keys_dict

        # PREPARE A FUNCTION FOR PUSHING COMPLETE DICTIONARY
        def _push_dict():
            for k, v in pro_atom_dict.items():
                pro_data_dict[k].append(v)

        # READ THE FILE TO THE LIST OF STRINGS, DELETE EMPTY
        file = open(path, "r", encoding='utf-8', errors='ignore')
        lines = file.readlines()
        lines.reverse()

        while lines:
            line = lines.pop()
            # ADD COMMENTS TO META
            if line[0] == '!':
                self.meta['comment'] += line
            # IF ATOM SECTION IS FOUND
            elif 'Atom ' in line and 'iZ =' in line:
                # IF IT'S DUMMY ATOM, SKIP IT
                if line.split()[3][0] == 'X':
                    continue
                # LOAD ALL NECESSARY INFORMATION FROM THE FIRST LINE
                pro_atom_dict = pro_keys_dict
                pro_atom_dict['ID']     = int(line.strip().split()[1])
                pro_atom_dict['Label']  = line.strip().split()[3]
                pro_atom_dict['iZ']     = int(line.strip().split()[6])
                # WHILE IN THIS SECTION, LOAD REST OF THE DATA
                while '-'*20 not in line:
                    line = lines.pop()
                    # LOAD SINGLE KEYWORD VALUES
                    for key in pro_keys_dict.keys():
                        if ' '+key+' ' in line:
                            pro_atom_dict[key] = float(line.strip().split()[-1])
                    # LOAD SOURCE FUNCTION VALUES
                    if 'Source function :' in line:
                        line = lines.pop()
                        for index in range(1, 9, 1):
                            try:
                                _sf = float(lines[-1].strip().split()[-1])
                                pro_atom_dict['_SF'+str(index)] = _sf
                            except (ValueError, IndexError):
                                break
                            else:
                                line = lines.pop()
                # AFTER FINDING DASHED LINE, PUSH DATA
                _push_dict()
                print(pro_atom_dict)

        # AFTER THE FILE IS READ, CLOSE IT
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

        # PREPARE LIST OF "@Other calculated information": '@SF%'
        pro_SF_sum = OrderedDict()
        try:
            for index in range(1, 9, 1):
                key = '_SF' + str(index)
                try:
                    pro_data_dict[key][0]
                except IndexError:
                    break
                pro_SF_sum[key] = 0.0
                for sf in pro_data_dict[key]:
                    pro_SF_sum[key] += sf
                pro_data_dict['_%SF' + str(index)] = []
                for sf in pro_data_dict[key]:
                    pro_data_dict['_%SF' + str(index)].append(sf / pro_SF_sum[key])
        except TypeError:
            pass

        # DELETE EMPTY KEYWORDS
        keys_to_delete = []
        for key in pro_data_dict.keys():
            print(pro_data_dict[key])
            if not pro_data_dict[key][0] and not pro_data_dict[key][-1]:
                keys_to_delete.append(key)
        for key in keys_to_delete:
            del(pro_data_dict[key])

        # SAVE .DATA FROM PRO_DATA_DICT
        self.data = pd.DataFrame.from_dict(pro_data_dict)

    def write(self, path):
        self.data.to_csv(path, decimal='.')


if __name__ == '__main__':
    pro = ProFrame(file_path='/home/dtchon/git/kesshou/test_data/xd_pro.out')
    pro.write('/home/dtchon/git/kesshou/test_data/xd_pro_digested.out')
