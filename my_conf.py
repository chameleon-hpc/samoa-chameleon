scenario='swe'
swe_dg_order='1'
#swe_patch_order=''
#swe_scenario='splashing_pool'
#swe_scenario='resting_lake'
#swe_scenario='gaussian_curve'
#swe_scenario='parabolic_isle'
#swe_scenario='linear_dam_break'
#swe_scenario='radial_dam_break'
swe_scenario='oscillating_lake'
dg_limiter='unlimited'
#dg_limiter='height'
flux_solver='aug_riemann'
#flux_solver='fwave'
exe="samoa_"+swe_dg_order+"_"+dg_limiter
compiler='gnu'
asagi='false'
mpi='nompi'
target='release'
#target='debug'
data_refinement='sample'
#debug_level=6

