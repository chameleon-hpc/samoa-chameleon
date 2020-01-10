# Scenario 4.3: Long wave resonance in a paraboloid basin

## Notes

The refinement level `15` was chosen, which results in a triangle leg length of about `88.39` - the same as in [1]: see section 4.3, paragraph two.

The dry/wet tolerance of `10e-2` resp. `10e-8` was directly taken from [1].

The courant number of `0.2` was found by using a estimated average of the values given in [1]. The average computation time step size is about `0.451` seconds, while the fixed time step size in [1] is `2.534` seconds.

The datasets only contain the timesteps at `t = n * 0.25P` from `t = 0` up to `t = 2P` to save space.

## Compilation and Execution

## Evaluation

A z-cross-section is interpolated at `z = 0` and plotted for all relevant variables, and compared to the analytical solution at `t = 1.5P, 1.75P, 2P`. These results can then be compared with [1], see figure 7 and 8.

Contour plots at `t = 2P` are generated for the full domain. These results can then be compared with [1], see figure 9.
