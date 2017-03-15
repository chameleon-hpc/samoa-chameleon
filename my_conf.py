scenario='swe'
swe_dg_order='2'
#swe_scenario='splashing_pool'
#swe_scenario='resting_lake'
#swe_scenario='gaussian_curve'
#swe_scenario='parabolic_isle'
swe_scenario='linear_dam_break'
dg_limiter='unlimited'
#dg_limiter='height'
flux_solver='aug_riemann'
exe="samoa_"+swe_dg_order+"_"+dg_limiter
compiler='gnu'
asagi='false'
mpi='nompi'
target='debug'
data_refinement='sample'
#debug_level=6

