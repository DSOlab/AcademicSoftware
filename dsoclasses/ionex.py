
class Ionex:
    def read_header(self, fn):
        with open(fn, 'r') as fin:
            line = fin.readline()
            self.version = float(line[0:8])
            assert line[20] == 'I'
            self.satsys = line[40:60].strip()
