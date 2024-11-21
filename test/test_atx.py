#! /usr/bin/python

import sys
from dsoclasses.gnss import atx

a = atx.Atx(sys.argv[1])
p = a.get_noazi("JAVTRIUMPH_3NRA", ["G01", "R02"])
print(p.pcv)
print(p.pco('G01'))
print(p.pco('R02'))
