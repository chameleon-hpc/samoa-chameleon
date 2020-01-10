# Scenario 4.4: Oscillatory flow in a parabolic bowl

## Notes

The refinement level `12` was chosen, which results in `64*64*2` elements - the same as in [1]: see section 4.4, paragraph two.

The dry/wet tolerance of `10e-3` was directly taken from [1].

The courant number of `0.3` was found by experimentation. The average computation time step size is about `0.002218` seconds, while the fixed time step size in [1] is `0.004487` seconds.

The datasets only contain the timesteps at `t = n * (P / 1000)` from `t = 0` up to `t = 2P` to save space.

A series of runs are made for refinement levels `10` to `18`, and written only at `t = 2P`.

## Compilation and Execution

## Evaluation

A z-cross-section is interpolated at `z = 0` and plotted for all relevant variables, and compared to the analytical solution at `t = 2P`. These results can then be compared with [1], see figure 11.

Contour plots at `t = 2P` are generated for the full domain. These results can then be compared with [1], see figure 12.

A limiter error series is plotted for different refinement levels. These results can then be compared with [1], see figure 14.

Mass and energy errors are calculated and plotted over time. These results can then be compared with [1], see figure 15.