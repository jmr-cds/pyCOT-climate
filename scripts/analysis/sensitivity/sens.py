"""Sensitivity probes for pyCOT hyperparameters (read-only use of the library).

Usage: python scripts/analysis/sensitivity/sens.py net1.txt [net2.txt ...]
Produces the table in docs/CODEBASE_ANALYSIS.md, section 6.
"""
import sys, time, json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'src'))
from pyCOT.io.functions import read_txt
from pyCOT.analysis.organizations.organizations import compute_organizations
from pyCOT.analysis.organizations.self_maintenance import check_self_maintenance
from pyCOT.analysis.organizations.closure import is_closed, is_ssm

nets = sys.argv[1:]
out = {}
for path in nets:
    rn = read_txt(path)
    rec = {}
    for mo in (2, 5, 10, 100):
        t = time.time()
        r = compute_organizations(rn, max_so_order=mo, verify_organizations=False)
        rec[f'n_so_max_order_{mo}'] = len(r.semiorganizations)
        rec[f'max_order_found_{mo}'] = max((s.order for s in r.semiorganizations), default=-1)
        rec[f't_{mo}'] = round(time.time() - t, 2)
    r2 = compute_organizations(rn, max_so_order=100, use_vertical_lift=False, verify_organizations=False)
    rec['n_so_no_lift'] = len(r2.semiorganizations)
    # closure / SSM re-check on quotiented+E0 raw masks
    rd = r.rn_data
    bad_closed = bad_ssm = 0
    for s in r.semiorganizations:
        m = s.species_mask | rd.E0_mask
        bad_closed += not is_closed(rd.supp_raw, rd.prod_raw, m)
        bad_ssm += not is_ssm(rd.supp_raw, rd.prod_raw, m)
    rec['n_not_closed'] = bad_closed
    rec['n_not_ssm'] = bad_ssm
    rec['n_ercs'] = len(r.ercs)
    rec['E0_size'] = bin(rd.E0_mask).count('1')
    # epsilon sensitivity of LP verdicts (cap number of SOs for runtime)
    sos = r.semiorganizations[:400]
    for eps in (1e-3, 1e-4, 1e-6, 1e-8, 1e-10):
        n = 0
        for s in sos:
            sp = [x for x in rn.species() if x.name in set(s.species_names)]
            ok, _, _ = check_self_maintenance(sp, rn, epsilon=eps)
            n += bool(ok)
        rec[f'n_org_eps_{eps:g}'] = n
    rec['n_checked_for_eps'] = len(sos)
    # module size with vs without E0 (order-0 only)
    el = [s for s in r.semiorganizations if s.order == 0]
    if el:
        rec['elem_size_with_E0_mean'] = round(sum(len(s.species_names) for s in el) / len(el), 2)
        rec['elem_size_without_E0_mean'] = round(sum(bin(s.species_mask).count('1') for s in el) / len(el), 2)
        rec['n_elementary'] = len(el)
    out[path.split('/')[-1]] = rec
    print(path, json.dumps(rec), flush=True)
