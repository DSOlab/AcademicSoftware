#! /usr/bin/python

from dsoclasses.orbits import interpolator
import sys

if __name__ == "__main__":
    
    sp3int = interpolator.Sp3Interpolator(sys.argv[1], ['G', 'J'])
