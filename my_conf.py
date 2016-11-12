scenario='swe'
#swe_patch_order='3'
swe_dg_order='2'
#swe_dg_basis='bernstein_l2'
#swe_dg_basis='bernstein_nodal'
#swe_scenario='parabolic_isle'
#swe_scenario='splashing_pool'
#swe_scenario='gaussian_curve'
#swe_scenario='resting_lake'
#swe_scenario='oscillating_lake'
#swe_scenario='radial_dam_break'
#swe_scenario='linear_beach'
#swe_scenario='linear_dam_break'
#swe_dg_basis='bernstein_nodal'
swe_scenario='convergence'
#flux_solver='aug_riemann'
flux_solver='hlle'
dg_limiter='unlimited'
exe="samoa_"+swe_dg_order+"_"+dg_limiter
compiler='gnu'
asagi='false'
mpi='nompi'
target='release'
data_refinement='sample'


