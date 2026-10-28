"""Inspect semi-organizations returned by compute_organizations that are not closed.

Usage: python scripts/analysis/sensitivity/nc.py [network.txt]
Default network: data/biochemical_databases/BiGG/bigg_e_coli_core.txt
See docs/CODEBASE_ANALYSIS.md, section 3 item 1.
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / 'src'))

from pyCOT.io.functions import read_txt
from pyCOT.analysis.organizations.organizations import compute_organizations
from pyCOT.analysis.organizations.closure import is_closed, is_ssm
from pyCOT.analysis.organizations.closure_oracle import closure_oracle

path = sys.argv[1] if len(sys.argv) > 1 else str(REPO / 'data/biochemical_databases/BiGG/bigg_e_coli_core.txt')
rn = read_txt(path)
r = compute_organizations(rn, max_so_order=100, verify_organizations=True)
rd = r.rn_data
nm = rd.species_names
masks = {s.species_mask for s in r.semiorganizations}
print('n_so', len(r.semiorganizations), 'n_org', len(r.organizations))
for s in r.semiorganizations:
    m = s.species_mask | rd.E0_mask
    if not is_closed(rd.supp_raw, rd.prod_raw, m):
        c = closure_oracle(rd.supp_raw, rd.prod_raw, m)
        extra = [nm[j] for j in range(rd.n_species) if (c & ~m) >> j & 1]
        cq = c & ~rd.E0_mask
        print('order', s.order, 'size', bin(m).count('1'), 'closure adds', extra,
              '| closure SSM', is_ssm(rd.supp_raw, rd.prod_raw, c),
              '| closure already listed', cq in masks, '| is_org', s.is_organization)
