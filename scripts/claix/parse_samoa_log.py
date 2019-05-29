import sys
import re


filename=sys.argv[1]

file=open(filename, 'r')
#     (r0,t0)  Phase time:                        3089.0989 s
phase_time_pattern = re.compile(".*Phase time: *([0-9]+\.[0-9]+)") 
phase_times = []
#     (r0,t0)  sam(oa): Space filling curves and Adaptive Meshes for Oceanic and Other Applications
#     (r0,t0)  Scenario: SWE
#     (r0,t0)  OpenMP: Yes, with tasks, threads: 23, procs: 48
#     (r0,t0)  MPI: Yes, ranks: 32
#     (r0,t0)  ASAGI: Yes, without NUMA support, mode: 2: no mpi
#     (r0,t0)  Debug Level: 1
#     (r0,t0)  Assertions: Yes
#     (r0,t0)  Precision: Double
#     (r0,t0)  Compiler: Intel
#     (r0,t0)  Sections per thread: 16
#     (r0,t0)  Adaptivity: min depth: 15, max depth: 25, start depth: 0
#     (r0,t0)  SWE: Patches: Yes, order: 7, vectorization: Off
#     (r0,t0)  Load balancing: timed load estimate:  No, split mode: 0, serial:  No, frequency: 100000 , threshold: .010
#     (r0,t0)  Load balancing: for heterogenous hardware (HH):  No, ratio: 1.0/1.0
#     (r0,t0)  Load balancing: cell weight: 1.00, boundary weight: .00
ranks_pattern = re.compile(".*ranks: ([0-9]+)")
sections_pattern = re.compile(".*Sections per thread: ([0-9]+)")
threads_pattern = re.compile(".*, threads: ([0-9]+),.*")
min_max_pattern = re.compile(".*min depth: ([0-9]+), max depth: ([0-9]+)")
lbfreq_pattern = re.compile(".* frequency: ([0-9]+)")

threads=-1
ranks=-1
dmin=-1
dmax=-1
sections=-1
lbfreq=-1

for line in file:
  m=re.match(phase_time_pattern, line)
  if m:
    phase_times.append(float(m.group(1)))
  m=re.match(ranks_pattern, line)
  if m:
    ranks = int(m.group(1))
    #print m.group(1)
  m=re.match(sections_pattern, line)
  if m:
    sections = int(m.group(1))
  m=re.match(threads_pattern, line)
  if m:
    #print line
    threads = int(m.group(1))
  m=re.match(min_max_pattern, line)
  if m:
    dmin = int(m.group(1))
    dmax = int(m.group(2))
  m=re.match(lbfreq_pattern, line)
  if m:
    lbfreq = int(m.group(1))

if(len(phase_times)<3):
  time = -1
else:
  time = phase_times[-1]

chameleon="no"
if 'chameleon' in filename:
  chameleon="yes"

print chameleon,ranks,threads,lbfreq,sections,dmin,dmax,time
