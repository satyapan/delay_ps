# Making delay power spectra from measurement set

## Import module
```
import sys
sys.path.append('/home/users/satyapan/repos/delay_ps')
from delay_ps import *
```

## Make delay ps generator
```
dpg = delay_ps_gen(ms_file, timeavg=False, data_col="DATA")
```
Arguments:
  - ms_file: Path to ms
  - data_col: Data column in ms, default = "DATA"
  - timeavg (bool): If True then 1 delay ps produced for entire ms. If False, multiple delay ps are produced separated by n_timeavg steps
  - n_timeavg: Needed only when timeavg = False. Number of time steps to average for a single delay ps.
  - stokes: Stokes parameter to use 


## Make delay ps
```
dps = dpg.get_delay_ps()
```

## Plot delay ps
```
dps.plot(plot_w_line=False, flip_w=False)
```
Arguments:
  - ax: Axis on which to plot delay ps. If no axis is given, a new figure is created and the delay ps is plotted on it.
  - plot_w_line (bool): If True, plots the modified horizon line for NCP/SCP phasing
  - flip_w (bool): If True, flips the delay spectra to have all baselines pointing along North (for NCP phasing)
