#! /usr/bin/python

import sys
from dsoclasses.sinex import sinex

snx = sinex.Sinex(sys.argv[1])
print("Sinex Header Info:")
print("version           {:}".format(snx.version))
print("agency            {:}".format(snx.agency))
print("created at        {:}".format(snx.time_creation))
print("data provider     {:}".format(snx.data_provider))
print("data start        {:}".format(snx.data_start))
print("data stop         {:}".format(snx.data_stop))
print("technique         {:}".format(snx.technique))
print("num estimates     {:}".format(snx.num_estimates))
print("constraint code   {:}".format(snx.constraint_code))
print("solution contents {:}".format(snx.solution_contents))

l = snx.site_id(['DIOA', 'DIOB', 'DYNG'])
print(l)
# l = snx.solution_epochs(['DIOA', 'DIOB', 'DYNG'])
# print(l)
l = snx.solution_estimate(['DIOA', 'DIOB', 'DYNG'])
print(l)
