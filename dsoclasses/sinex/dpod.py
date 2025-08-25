
def dpod_freq_corr(fn: str, site_list=[]):
    result = {}
    # result = { 'ADEA': [{'soln':1, 'freq': 365.25, 'xyz': {'xcos': 2.263, 'xsin': 123., 'ycos': ...}}, 
    #                     {'soln':1, 'freq': 182625, 'xyz': {'xcos': 2.263, 'xsin': 123., 'ycos': ...}}]}
    def append_component(site, freq, soln, xyz, cosamp, sinamp):
        def cmp(xyz): return xyz.lower()+'cos', xyz.lower()+'sin'
        if site not in result:
            result[site] = [{'soln':soln, 'freq': freq, 'xyz': {cmp(xyz)[0]: cosamp, cmp(xyz)[1]: sinamp}}]
            return
        else:
            for j, sf in enumerate(result[site]):
                if sf['soln'] == soln and sf['freq'] == freq:
                    result[site][j]['xyz'][cmp(xyz)[0]] = cosamp
                    result[site][j]['xyz'][cmp(xyz)[1]] = sinamp
                    return
            result[site].append({'soln':soln, 'freq': freq, 'xyz': {cmp(xyz)[0]: cosamp, cmp(xyz)[1]: sinamp}})
        return
    
    cur_freq = -999.

    with open(fn, 'r') as fin:
        line = fin.readline()
        while line:
            if line.startswith('# Frequency'):
                # Frequency  1 : 365.250 days
                l = line.split()
                cur_freq = float(l[4])
                assert(l[5] == 'days')
            elif line.startswith('#'):
                pass
            else:
                ##CODE PT __DOMES__SOLN_XYZ_COSAMP__COSSTD__SINAMP__SINSTD                       
                # ADEA  A 91501S001  1   X   2.263   0.077   0.464   0.078
                l = line.split()
                assert(len(l) == 9)
                if l[0] in site_list or site_list == []: 
                    append_component(l[0], cur_freq, l[3], l[4], float(l[5]), float(l[7]))
            line = fin.readline()
    return result
