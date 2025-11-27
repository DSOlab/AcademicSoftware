# Frequencies in MHz
GPS_L1_FREQ = 1575.42
GPS_L2_FREQ = 1227.60
GPS_L5_FREQ = 1176.45

GAL_E1_FREQ  = 1575.42
GAL_E5a_FREQ = 1176.45
GAL_E5b_FREQ = 1207.140
GAL_E5_FREQ  = 1191.795 # (E5a+E5b)
GAL_E6_FREQ  = 1278.75

GLO_G1_FREQ = 1602. #+k*9/16
GLO_G2_FREQ = 1246. #+k*7/16
GLO_G3_FREQ = 1202.025
def glonass_sat_freqs(slot_nr):
    k = slot_nr
    return GLO_G1_FREQ+(k*(9./16.)), GLO_G2_FREQ+(k*(7./16.)), GLO_G3_FREQ


# Speed of light m/sec
C = 299792458e0

# Rotation Rate of Earth (along +Z axis) in [rad/sec]
OmegaEarth = 7.2921159e-5
