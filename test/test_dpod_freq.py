#! /usr/bin/python

from dsoclasses.sinex.dpod import dpod_freq_corr
import sys

print(dpod_freq_corr(sys.argv[1], ['DIOB', 'ADEA']))
