scenario='swe'
swe_dg_order='2'
#swe_patch_order=''
#swe_scenario='splashing_pool'
#swe_scenario='resting_lake'
#swe_scenario='gaussian_curve'
#swe_scenario='parabolic_isle'
#swe_scenario='linear_dam_break'
swe_scenario='radial_dam_break'
#swe_scenario='oscillating_lake'
dg_limiter='all'
#dg_limiter='height'
flux_solver='aug_riemann'
flux_time_averaging=True
#flux_solver='fwave'
exe="samoa_"+swe_dg_order+"_"+dg_limiter
compiler='intel'
asagi='false'
#mpi='mpi'
#mpi='mpi'
target='release'
#target='debug'
data_refinement='sample'
#debug_level=6
chameleon=0
