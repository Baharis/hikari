class LstFrame:
    def __init__(self):
        pass

    @staticmethod
    def read_r1(path):
        """Read and return the final value of R1 from lst file"""

        # READ THE FILE AND JOIN LINES SEPARATED BY '=' SIGN
        lines = [line.strip() for line in open(path, 'r') if line.strip()]
        for index, line in enumerate(lines):
            while line[-1:] == '=':
                lines[index] = line[:-1] + lines[index+1]
                del(lines[index+1])

        # READ THE LINES, LOOK FOR R1, RETURN THE LAST FOUND
        r1 = ''
        for line in lines:
            if 'R1 =' in line:
                r1 = line.split()[2]

        return r1
