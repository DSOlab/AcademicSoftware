#! /usr/bin/python

from dsoclasses.orbits import sp3c

if __name__ == "__main__":
    # validate using a DORIS sp3c file, with only one satellite
    sp3 = sp3c.Sp3('ssaja320.b24212.e24222.DG_.sp3.001')
    print('Number of sats in sp3 file: {:}'.format(len(sp3.sat_ids)))
    print('Satellites: {:}'.format(' '.join(sp3.sat_ids)))
    print('SP3 time scale: {:}'.format(sp3.time_sys))
    print('SP3 Reference frame: {:}'.format(sp3.crd_system))
    print('SP3 Orbit type: {:}'.format(sp3.orb_type))
    print('Agency : {:}'.format(sp3.agency))
    print(sp3.get_satellite('L39'))
    
    # validate using a GNSS sp3c file, with only multiple satellites
    sp3 = sp3c.Sp3('COD0MGXFIN_20241050000_01D_05M_ORB.SP3')
    print('Number of sats in sp3 file: {:}'.format(len(sp3.sat_ids)))
    print('Satellites: {:}'.format(' '.join(sp3.sat_ids)))
    print('SP3 time scale: {:}'.format(sp3.time_sys))
    print('SP3 Reference frame: {:}'.format(sp3.crd_system))
    print('SP3 Orbit type: {:}'.format(sp3.orb_type))
    print('Agency : {:}'.format(sp3.agency))
    print(sp3.get_satellite('G03'))
    print(sp3.get_satellite('R11'))
    print(sp3.get_satellite('E07'))
    print(sp3.get_satellite('C21'))
    print(sp3.get_satellite('J04'))
