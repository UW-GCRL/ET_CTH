"""R1 analysis 6: transpiration and GPP responses by PFT (Fig. S6c; Section 3.1).
Relative change of forest-area-weighted transpiration (FCTR) and GPP (FPSN) for the 20% height increase
(Default 1.2) and for GEDI Mean, each relative to CLM Default. Output: an6_phs_diag.json."""
import json
import numpy as np
from common import *

m = load_meta(); f = m['forest']
w = m['wt'] * m['area'][m['jx'], m['ix']]
D = {s: load_scen(s, ['FCTR', 'FPSN', 'QVEGE', 'QVEGT', 'QFLX_EVAP_TOT', 'QFLX_LIQ_GRND', 'QFLX_SNOW_GRND',
                      'FCEV', 'FGEV']) for s in ['CLM Default', 'Default 1.2', 'GEDI Mean']}


def rel(a, v, sel):
    return float(100 * (np.nansum((D[a][v] * w)[sel]) / np.nansum((D['CLM Default'][v] * w)[sel]) - 1))


res = {}
for name, sel in [('All forest', f)] + [(p, f & (m['iv'] == k)) for p, k in PFTS.items()]:
    r = dict(dT_p20=rel('Default 1.2', 'FCTR', sel), dGPP_p20=rel('Default 1.2', 'FPSN', sel),
             dT_mean=rel('GEDI Mean', 'FCTR', sel), dGPP_mean=rel('GEDI Mean', 'FPSN', sel))
    r['ratio_p20'] = r['dGPP_p20'] / r['dT_p20']
    r['ratio_mean'] = r['dGPP_mean'] / r['dT_mean']
    res[name] = r
    print(f"{name:10s} +20%: T {r['dT_p20']:+.2f}%  GPP {r['dGPP_p20']:+.2f}%   GEDI Mean: T {r['dT_mean']:+.2f}%  GPP {r['dGPP_mean']:+.2f}%")
json.dump(res, open('an6_phs_diag.json', 'w'), indent=1)
