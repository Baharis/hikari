import pandas as pd
import os


class CpFrame:
    def __init__(self, file_path=None):
        self.data = pd.DataFrame()
        self.meta = dict()
        if file_path is not None:
            self.read(path=file_path)

    def read(self, path):
        """Read Gcp_Vcp.dat file of MoPro as specified by path
        and update CPFrame's data and meta accordingly"""

        # SET DAT FILE READING FORMAT AND PARSING OBJECTS
        indices = (0, 2, 4, 5, 6, 7, 8, 9, 10, 13, 14, 18, 19)
        labels = ('atom1', 'atom2', 'cp', 'sym',
            'gcp_atomic', 'vcp_atomic', 'gcp', 'vcp',
            'len', 'den', 'lap', 'eli', 'type')

        # OPEN FILE AND PREPARE HOLDERS
        with open(path, 'r') as file:
            self.meta['raw'] = file.readlines()
        meta_content = str()
        data_content = list()
        flag = 'meta'

        # HANDLE LINE DEPENDING ON CURRENT FLAG AND CONTENT
        for line in self.meta['raw']:
            if not line.strip():
                pass
            elif flag == 'meta':
                meta_content += line
                if line.strip().split()[0:2] == ['Atom1', 'Atom2']:
                    flag = 'data'
            elif flag == 'data':
                line = line.rstrip('\n').replace('-', ' -').replace(', -', ',-')
                line = line.split()
                line = [line[index] for index in indices]
                data_content.append(line)

        # TRANSFORM DATA_CONTENT TO PANDAS DATAFRAME
        self.data = pd.DataFrame(data_content, columns=labels)

        # TRANSFORM META_CONTENT TO META
        self.meta['comment'] = meta_content
        self.meta['path'] = os.path.abspath(path)

    def concentrate(self):
        indices, unique = list(), list()
        for index, row in self.data.iterrows():
            fingerprint1 = (row['atom1'], row['atom2'], row['den'])
            fingerprint2 = (row['atom2'], row['atom1'], row['den'])
            if (fingerprint1 not in unique) and (fingerprint2 not in unique):
                indices.append(index)
                unique.append(fingerprint1)
        self.data = pd.DataFrame(self.data, index=indices)


if __name__ == '__main__':
    path = '/home/dtchon/git/kesshou/test_data/Gcp_Vcp.dat'
    a = CPFrame()
    a.read(hkl_path=path)
    a.concentrate()
    a.data.to_csv('/home/dtchon/git/kesshou/test_data/Gcp_Vcp.csv')
