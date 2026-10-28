# pyCOT-climate: codebase, hyperparameters, and project analysis

Snapshot: branch `claude/admiring-mccarthy-umrppx`, 2026-10-28. Line numbers refer to this snapshot.

**Epistemic tags used throughout**
- **[measured]**: reproduced by running the code in this snapshot (Python 3, `rustworkx`, `scipy`/HiGHS, installed at the versions current on 2026-10-28). Script in §6.
- **[read]**: verified by reading the code or a committed output file.
- **[reported]**: stated in a project README/report/paper draft. Supporting output not found or not checked.
- **[plausible]**: inferred from code reading, consistent with measurements, but not proven by a dedicated test.

---

## 1. Architecture (src/pyCOT)

### 1.1 Data model and I/O
- **`core/rn_rustworkx.py`, `ReactionNetwork(PyDiGraph)`**: a bipartite directed graph. Species and reaction nodes share one index space. Edges carry `ReactionEdge(coefficient, type ∈ {reactant, product})` (`core/rn_types.py`).
- **Stoichiometry matrix**: `S = products − reactants`. A catalyst therefore has a net-zero entry, so catalysis is visible only in the *support* (bitmask) representation, not in S.
- **`get_reactions_from_species`**: always includes every inflow reaction.
- **`generated_closure`**: seeds with `S ∪ inflow species`. The LP self-maintenance check therefore evaluates `closure(X ∪ inflows)`, not `X` [read].

**Text format** (`io/functions.py:47 read_txt`):
- One reaction per line: `name: a + 2 b => c;`.
- Everything after the first `;` is stored as a comment.
- An empty LHS is an inflow; an empty RHS is an outflow/decay.
- `#` begins a comment.
- Coefficient forms: `2*A`, `2 A`, `2A`, `1.496 3pg_c`.
- Pitfalls:
  - Compact `2A` parsing turns names like `2OG` or `3PG` into coefficient+species [read].
  - A leading-decimal coefficient `.3*B` is mis-parsed. The repo's own test fails on it: `test_add_from_reaction_string_mixed_coefficients` [measured].
  - A file with several `;`-separated reactions on one line (e.g. `data/Examples_tests/LUCA/KEGG_data.txt`) loads only its first reaction via `read_txt` [read].

### 1.2 Organization pipeline (`analysis/organizations/`)
Entry point is `organizations.compute_organizations` (`organizations.py:462`). It runs these stages:
1. **`build_rndata`** (`io_pyCOT.py`): compiles supports and products to integer bitmasks.
   - **E0** = closure of the inflow products, i.e. everything producible from food alone.
   - All later combinatorics run on *quotiented* masks with E0 removed (`supp_q = supp & ~E0`).
   - Stoichiometric coefficients are ignored until the LP.
2. **ERCs** (`erc.py`): for each reaction, `closure(supp_q[r])`. Reactions are grouped by the closure they generate. `min_bases` are the ⊆-minimal supports in each group. E0 is injected as a persistent ERC.
3. **Hierarchy** (`hierarchy.py`): strict-subset order on ERC masks, with Hasse reduction.
4. **Synergies** (`synergy.py`, `compute_synergies_basis_first`): (i, j) → k when a minimal basis of k is covered by `E_i ∪ E_j` but by neither alone. Results are filtered to *maximal*, then *fundamental* (no descendant pair yields the same k).
5. **Complementarities** (`complementarity.py`): `prod_i ∩ req_j ≠ ∅` (or the reverse). Also filtered to *fundamental* (⊆-minimal producer/consumer per species).
6. **SO0 search** (`so_search.compute_elementary_sos`): Mode-1 DFS in ERC-index space.
   - A state is an ERC set, with `req(S) = ∪req \ ∪prod` and `prod(S) = ∪prod`.
   - Moves:
     - (a) add a ⊆-minimal producer of a required species;
     - (b) add a synergy partner whose target produces something required.
   - Canonical-ordering rules and Pareto dominance prune symmetric paths.
   - A state with `req = ∅` is semi-self-maintaining (SSM).
7. **SOi search** (`so_search.compute_so_hierarchy:781`): BFS by "order". From each known SO, Mode-2 seeds are fundamental synergy partners, consumers of produced species, and (optionally) Hasse parents ("vertical lift"); each seed is followed by Mode-1. Order is assigned by containment among the SOs that were *found*.
8. **LP self-maintenance** (`self_maintenance.minimize_sv:34`): a feasibility check, solved with HiGHS.
   - Constraints: `v_r ≥ ε` for every reaction r ∈ R_X, `Σv = |R_X|`, `S·v ≥ −margin`.
   - An SO that passes is an **organization**.

Other modules:
- **Reference oracles**: `closure_oracle.py`, plus `closure.is_closed` and `closure.is_ssm`.
- **`max_semiorg.py`**: iteratively removes ERCs whose requirements cannot be met.
- **`generators.py`**: primitive (non-synergy-target) ERCs.
- **Legacy, parallel implementations**: `analysis/ERC_Hierarchy.py`, `ERC_Synergy_Complementarity.py`, `Persistent_Modules_Generator.py` (brute force over ERC combinations; only *pairs* for non-elementary SOs).

### 1.3 Decomposition theorem (`analysis/decomposition/`)
Implements X = (E ∪ F) ∪ D₁ ∪ … ∪ Dₘ, where:
- E = catalysts (zero rows of S over R_X; `catalysts.py`);
- F = overproducible species (one LP per species, maximizing (S·v)_s with `v ≥ 0`, summed into a joint witness; `overproduction.py:32-87`);
- Dᵢ = fragile circuits (union-find over the remaining species; `circuits.py`).

X is self-maintaining iff every Dᵢ is self-maintaining on its own consuming reactions. `hierarchy.decompose_hierarchy` reuses parent results along the SO lattice. `dependency.py` builds the interface DAG.

### 1.4 Kinetics and simulation
- **Rate laws** (`kinetics/`, registry in `kinetics/__init__.py`): `mak`, `mmk` (first reactant only), `hill` (S = sum of reactants), `cosine`, `saturated`, `threshold_memory`, `inhibition`.
- **ODE** (`simulations/ode.py:44`): `solve_ivp` with LSODA by default. On failure it silently retries at looser tolerances, then falls back to `odeint`. The state is clamped at ≥ 0 both inside the RHS and afterwards.
- **Stochastic**:
  - `simulations/stochastic.gillespie`: direct-method SSA with combinatorial mass-action propensities.
  - `kinetics/stochastic_events.GillespieSSA`: a *second* SSA with deterministic-law propensities (k·xⁿ). The two disagree whenever a stoichiometric coefficient is > 1 [read].
  - `TauLeaping` and `LangevinDynamics` are empty stubs.
  - `BudgetAllocationSimulator`: a discrete-time Poisson budget model, used by XCEPT Article 2.
- **Spatial** (`simulations/spatial.py`): 2-D grid, 5-point Laplacian, periodic boundaries.
- **Metapopulation** (`simulations/metapopulation.py:780`):
  - Each patch runs as an ODE or SSA model.
  - Patches are coupled by Lie–Trotter splitting, with transport `P = expm(Lᵀ dt)`.
  - The coupling step defaults to the output step (`dt_couple = dt_out`), so `n_steps` changes the results through O(dt) splitting error [read].

### 1.5 Process analysis (`analysis/process_analyzer.py`, `process_structure.py`)
**Mode classification** (`classify_process_mode:122`). With s = S·v and an absolute tolerance τ:
- *Stationary*: all |sᵢ| ≤ τ.
- *Overproduction*: s ≥ −τ with some sᵢ > τ.
- *Problem*: s ≤ τ with some sᵢ < −τ.
- *Challenge*: mixed signs.
- The "cognitive domain" is (Stationary or Overproduction) with every vⱼ > τ.

**Multiscale analysis**:
- Windowed sums of v (aggregation).
- Random decompositions of v (disaggregation: coordinate, Dirichlet, or mixed).
- Mode 1: pathway variability.
- Mode 2: temporal scales.
- Mode 3: both, compared by z-score.
- `detect_thresholds` looks for jumps in category statistics across window sizes.

**Cone geometry**: sympy nullspace, plus a dense grid of `grid_resⁿ` points filtered by S·v ≥ −tol. There is a second, different `analyze_cone` in `process_analyzer.py`.

### 1.6 scripts/ and visualizations/
- **`scripts/analysis`**: test-suite-style drivers for each analysis module.
- **`scripts/simulations`**: ODE/SSA/spatial/metapopulation drivers.
- **`scripts/demos`**: stoichiometry toy examples.
- **`scripts/experimental`**: scratch work. `script9_1/9_2_*` are the cleanest examples and the only ones with valid `data/` paths.
- Most scripts reference `Txt/` or `networks/` directories that no longer exist; some import removed modules (`pyCOT.ERC_synergy`, `pyCOT.simulations_old`). Outputs are written relative to the current working directory.
- **Top-level `visualizations/`**: two HTML ERC/organization hierarchy views (`MLM_Afidos`, `Test`). They are produced by `rn_visualize.hierarchy_visualize_html`, called from `projects/COT_Fundamental_Generators_Exploration/gen_v2/run_organizations.py`.

---

## 2. Hyperparameters and their effects

### 2.1 Exhaustiveness of organization enumeration (largest effect on counts)

| Parameter (default) | Location | Effect |
|---|---|---|
| `max_so_order` / `max_order` (10) | `organizations.py:466`, `so_search.py:789,902` | Caps the number of Mode-2 BFS rounds; stops silently. **[measured]** BIOMD0000000943: 94 / 416 / **1224** / 1267 SOs at caps 2 / 5 / **10** / 100. The default loses 43 SOs (3.4%). e_coli_core saturates at cap ≥ 10 (309). Order labels can exceed the cap: labels come from containment, while the cap limits rounds. |
| `use_vertical_lift` (True) | `organizations.py:467`, `so_search.py:790` | Switches between the paper's Conjecture 1 (True) and Conjecture 2 (False). Neither is proven complete. **[measured]** False finds 252 vs 309 SOs on e_coli_core (−18%) and 84 vs 1267 on BIOMD0000000943 (−93%). Conjecture 2 is therefore not a safe approximation. |
| Fundamental-only moves (by design) | `so_search.py:1123-1161`; `organizations.py` module docstring | Latent joins (unions of SOs with no synergy/complementarity edge) and free-species additions are never explored unless `include_spurious=True`. The output is "fundamental" SOs, **not** all closed SSM sets. |
| Canonical-ordering pruning (no flag) | `so_search.py:553,564,583,593` | Symmetry breaking under an unproven assumption that every SSM can be built in increasing ERC-index order. A real risk of missing SOs [plausible]. The code itself documents a previously "silently dropped" class that was later fixed (`so_search.py:807-819`). |
| `include_spurious` (False), `exhaustive_cap` (10), `pairwise_cap` (500), `max_rounds` (20) | `organizations.py:469,159,293-294` | With spurious enumeration on, mixed free-species subsets are skipped above 10 free species, and pairwise latent-join probing is skipped above 500 masks. Both happen silently. |
| `verify_organizations` (True) | `organizations.py:468` | False skips the LP; every SO is then flagged `is_organization=False`. |
| `compute_synergies(level="maximal")` | `synergy.py:478` | The default level returns no *fundamental* synergies. `direct=True` and the staged filter give different fundamental sets, even though their docstrings say the results are identical [read]. |
| Legacy caps: `max_generator_size` (8), `max_organization_size` (5, unused), `max_combo_size` | `Persistent_Modules_Generator.py:563,627,793` | Generator size caps make the legacy enumeration non-exhaustive. Non-elementary SOs are only pairs, and are never LP-checked. |
| `max_chains` (5000) | `decomposition/org_graph.py:88` | Chain enumeration truncates silently. |

### 2.2 Numerical tolerances

| Parameter (default) | Location | Effect |
|---|---|---|
| LP `epsilon` (1e-6), the minimum flux per reaction | `self_maintenance.py:34,112` | **The most consequential numerical setting [measured].** With ε ≥ 1e-6, verdicts are stable (identical at 1e-3, 1e-4, 1e-6). With ε ≤ 1e-8, **every** SSM is declared an organization: e_coli_core 13 → 309/309, BIOMD0000000943 23 → 400/400, Farm.txt 4 → 14/14. Cause: a forced drain of order coefficient × ε falls below HiGHS's primal feasibility tolerance (~1e-7; the module docstring says this at l.59-65), so infeasible systems are reported feasible. The safety margin is only one order of magnitude; networks with fractional coefficients < 0.1 erode it [plausible]. ε too large risks false negatives if a valid flux needs dynamic range > 1/ε. |
| `margin` (1e-9) | `self_maintenance.py:34` | Slack on S·v ≥ −margin. Effectively dominated by the solver tolerance, since `S·v` is not re-checked after solving. |
| Normalization `Σv = n_reactions` | `self_maintenance.py:87` | Sets the scale that ε is compared against (mean flux 1). |
| Exceptions → "not SM" | `self_maintenance.py:107,156` | A solver crash is indistinguishable from a negative verdict. |
| Overproduction `_EPS` (1e-7), `_MARGIN` (1e-9), bounds `v ≥ 0` | `decomposition/overproduction.py:28-51` | `_EPS` thresholds the objective (S·v)_s, not a flux floor. The substantive issue is `v ≥ 0`: reactions in R_X may be switched off, which is the flaw `minimize_sv` was fixed for (`self_maintenance.py:46-58`). F can be over-inclusive, which shrinks the fragile circuits and could make the decomposition accept an X that the whole-set LP rejects [plausible, not measured]. |
| `_TOL` (1e-9) | `catalysts.py`, `decomposition/hierarchy.py` | Zero-row test for catalysts; sign test for Prop 2.19. |
| Process `tol` (1e-3, absolute) | `process_analyzer.py:127,210,365,…` | Not scaled by flux magnitude or window length. Large fluxes or wide windows leave almost no "Stationary"/"balanced" results; small fluxes are mostly "Stationary"/"Incomplete". Mode labels and proportions are scale artefacts unless τ is set relative to ‖S·v‖. |
| Cone `grid_res` (5 or 20), feasibility tol (1e-6) vs classification tol (1e-3) | `process_structure.py:614-616,675` | `grid_resⁿ` points; e.g. 25⁵ ≈ 9.8 M in `script_process_cone_geometry.py`. Coverage and counts depend on the grid. The two tolerances are inconsistent. |
| ODE `rtol`/`atol` (1e-8/1e-10), `method` (LSODA) | `ode.py:51-53` | The fallback loosens tolerances to 1e-6/1e-8 without raising an error. Clamping to ≥ 0 breaks smoothness and conservation. |
| SSA `max_iter` (100 000) | `stochastic.py:174,744` | When reached, the last state is held to `t_end`. A truncated run then looks like a steady state; there is only a verbose summary, no warning. |

### 2.3 Modelling and parsing flags

| Parameter (default) | Location | Effect |
|---|---|---|
| `exact_names` (False) | `rn_rustworkx.py:189,273`; `functions.py:47` | Species differing only in case or whitespace are **merged** (`ATP`/`atp`), which changes network topology. Endosymbiosis uses `exact_names=True`; most other callers do not. |
| Default rate constants U(0.01, 1), initial conditions U(50, 100), etc. | `simulations/core.py:163-191`, `ode.py:114`, `stochastic.py:426` | Unseeded random defaults. Any result produced from defaults is non-reproducible and reflects one random parameterization. |
| `dt_couple` (= `dt_out`) | `metapopulation.py:598-628` | Operator-splitting error is tied to output resolution. |
| `cosine` period | `deterministic_advanced.py:54` | Implements cos(w·t / 2π) rather than the documented cos(w·t); the period is 4π²/w. |
| `threshold_memory` memory/decay | `deterministic_advanced.py:133-265` | Memory state is never passed by the simulators, so μ and decay have no effect. |
| `BudgetAllocationSimulator` `dt` (0.05), `climate` (0.5), `rng_seed` | `stochastic_events.py:456-550` | Large dt pushes per-unit probabilities above 1; they are renormalised with a warning. Degradation amplitude is `1/climate`, not the documented `1−climate`, which is 4× asymmetric at 0.5 and divides by zero at 0. |

### 2.4 Sampling and statistics in process analysis

| Parameter (default) | Location | Effect |
|---|---|---|
| `n_samples` (100), `n_components` (10), `strategy` ('mixed') | `process_analyzer.py:622-625` | Coverage of the decomposition space. `pathway_variance` scales with `n_components`. `pattern_robustness` is computed on the original process every time, so it is always 0 or 1 [read]. |
| Mode 3 `n_decomp_samples` (50), `n_decomp_components` (5), z-cutoff 2 | `process_analyzer.py:830-834,954` | Compares per-step means (∝ agg/window) with component means (∝ agg/n_components). Unless window = n_components, \|z\| is large regardless of the dynamics. |
| `stride` (1) | `process_analyzer.py:736` | Windows overlap, so samples are autocorrelated and dispersion is understated. |
| `threshold_sensitivity` (0.5), `min_confidence` (0.6) | `process_analyzer.py:993-994` | Flags dᵢ > 0.6·mean(d). Because windowed sums grow ~linearly with window size, "thresholds" are largely artefacts. `variance_spike` is a no-op. |
| `seed` arguments | many | Often not propagated to default parameters or initial conditions. Global `np.random.seed` calls (`spatial.py:132`, `process_analyzer.py:1204`, `process_structure.py:295`) have side effects on unrelated code. |

---

## 3. Correctness findings that bear on reported results

1. **Reported SOs are not always closed [measured].**
   - 8 of the 309 SOs on e_coli_core fail `is_closed`. All 8 lack `succoa_c`. Its producer R33 (`akg_c + coa_c + nad_c`) needs three reactants, so the likely mechanism is that ERC-space synergy propagation is pairwise only (`fundamental_graph.py:167-267`) [plausible].
   - The closure of each such set is SSM and is *already* in the list. These are duplicates that inflate SO counts (to 309 from 301 distinct closed sets) and order labels (orders 3–7).
   - They are not mislabelled organizations, because the LP re-closes before checking.
   - Nothing in `compute_organizations` re-checks closure. A `closure == mask` assertion would catch this.
2. **The core SO engine and the `gen_v2` engine disagree [read].** In `projects/COT_Fundamental_Generators_Exploration/outputs/gen_v2_comparison/old_vs_new_80_800.csv`:
   - e_coli_core: 309 (old) vs 334 (new);
   - BIOMD0000000943: 1267 vs 2355.
   - Elementary counts agree.
   - The core engine in this snapshot reproduces the "old" numbers exactly [measured: 309 and 1267 at `max_so_order=100`]. Which count is correct is unresolved. Every project that uses the core hierarchy inherits this uncertainty.
3. **Order labels are relative to what was found.** A missed sub-SO lowers the order of its supersets, and can make a set look elementary.
4. **Label mismatch in process analysis [read].**
   - `classify_process_mode` returns `"Problem Mode"`/`"Challenge Mode"`. `classify_response_to_disturbance` compares against `"Problem"`/`"Challenge"` (`process_analyzer.py:290,296`), so those branches never fire.
   - Several plotting functions also look up `"Cognitive Domain"`/`"Problem"` keys that are never produced, so bars and proportions come out empty or zero.
   - `classify_response_to_disturbance` returns `None` when `verbose=False` on some paths.
5. **Legacy `ERC_Synergy_Complementarity.py` inverts containment** (`get_contained` returns ancestors). It yields minimal rather than maximal synergy targets, and one complementarity type can never fire. Results from this module differ from `organizations/`; it is used by `projects/Generative_Structure_Orgs` [read].
6. **Tests.**
   - The `src/pyCOT/tests` suite: 124 pass, 2 fail (decimal-coefficient parse; `Reaction.__str__` format). `test_load_from_file.py` errors at collection because of a hard-coded `Txt/` path [measured].
   - There are **no unit tests** for `analysis/organizations` or `analysis/decomposition` under `src/`. Project-local tests and oracles exist in `projects/COT_Fundamental_Generators_Exploration/tests`, `oracles/`, and `projects/Decomposition_Theorem/tests`.

---

## 4. How the projects use pyCOT to reach their claims

### 4.1 COT_Fundamental_Generators_Exploration (E. coli comparative paper)

**Pipeline:** `read_txt` → `build_rndata` → `compute_ercs` → `build_hierarchy` → `compute_synergies_basis_first` → `compute_complementarities` → `so_search.compute_elementary_sos` → `check_self_maintenance(ε=1e-6)`.
- Food sets are defined per scenario: aerobic/anaerobic/carbon-starvation core media, trace metals, and native media from BiGG exchanges.
- Input: `data/biochemical_databases/BiGG/bigg_*.txt`.

**Claims:**
- Reproduction of Centler et al. 2006.
- Module-size invariance across reconstructions.
- Trace-metal "fusion".
- More modules under anaerobic conditions.
- EO (self-maintaining) fraction of 29–40% (core) vs 11–15% (genome scale).
- Yeast modules are 100% compartment-hybrid.

**Assessment:**
- **Module size includes E0 [read + measured].**
  - `eso_eo_analysis.py:80` sets `full_mask = mask | E0_mask`.
  - In `outputs/ecoli_epm_comparison/results.csv`, `epm_size_min` equals `e0_size` in every food scenario (13/13, 11/11, 12/12).
  - On e_coli_core, mean elementary-SO size is 15.13 with E0 vs **2.13 without** [measured].
  - Consequences [read, from the same CSV]:
    - **"Size invariance" across reconstructions.** Under the aerobic core medium, iAF1260, iJO1366 and iML1515 all have E0 = 23 and mean size 25.7/25.7/25.6. The non-E0 excess is ≈ 2.6–2.7 in every case. The invariance therefore holds for the non-E0 part too, but it is a ~2.7-species signal, not a ~25-species one.
    - **Trace-metal / full-native "fusion".** iAF1260 mean size goes from 25.7 (aerobic core) to 63.7 (full native). E0 goes from 23 to 61, so the non-E0 excess stays 2.7. The near-tripling is E0 growth. The structural change is the module *count* (182 → 164). The same holds for iJO1366 (E0 71, mean 73.7) and iML1515 (E0 67, mean 69.7).
- **Yeast "100% hybrid" is tautological.** E0 contains both `_c` and `_m` species (co2, h2o, h), so every module spans both compartments regardless of topology [read].
- **Centler reproduction** needed a hand-added free-species step on top of the legacy brute force [read]. So "exact reproduction" depends on the spurious-organization extension, not on the fundamental engine alone.
- **Oracle validation**: `scripts/compare_oracles.py` covers 10 small networks; the ESO comparison columns are blank. ESO search therefore lacks an independent oracle check [read].
- **EO fractions** depend entirely on the ε setting (§2.2). No ε sensitivity analysis is reported.

### 4.2 COT_Endosymbiosis

**Claims:** fusigenic inflows exist; all 4 mergers are 100% hybrid; mean module size rises 54–87%; two mechanisms (capability gain vs re-integration); 1 of 7 interface reactions is load-bearing.

**Pipeline:** the same core calls (`endosymbiosis_analysis.py:66-84`, `exact_names=True`). `merge_networks.py` suffixes symbiont species with `__endo`. Hand-built interfaces and food sets: `real_data/build_*_case.py`. Fusion search is random sampling with seed 0 (`fusigenic_inflow_search.py:126-137`).

**Assessment:**
- **"100% hybrid" and the size increase follow by construction [read].** The merged E0 contains both partners' food, and module size includes E0 (log: merged sizes 26–30 = E0 26 + 2–4; host alone 13–17 = E0 13 + 2–4).
- The "mitochondrial" case uses e_coli_core as both host and symbiont.
- n = 3 hand-built cases plus a toy model, and interface reactions are modelling choices (acknowledged by the authors) [reported].
- The iAF1260 inflow search covers 15 hand-picked candidates.

### 4.3 Decomposition_Theorem

**Claim:** Theorem 2.16 (Veloz & Razeto-Barry 2017): an SO is self-maintaining iff each fragile circuit is.

**Pipeline:** `compute_so_hierarchy` (`SO_MAX_ORDER=100`) → `cot_gen` lattice → `decompose_hierarchy`. `run_decomposition.py:49-72` reassigns `NETWORK` several times; the last assignment wins (`Farm_agro_stages.txt`).

**Assessment:**
- Validated on one network (`tests/validate_decomposition.py`, e_coli_core).
- Inherits the `v ≥ 0` overproduction issue (§2.2) and the core-vs-gen_v2 SO count discrepancy (§3.2).
- No theorem write-up in the repo.

### 4.4 RAF_Comparison

**Purpose:** "Decomposing Autocatalytic Sets through COT", §4. The paper source is not in the repo.

**Pipeline:**
- Heteropolymer network generator: `MAX_LENGTH=3`, `P_CATALYST=0.2`, `CLEAVAGE_PROB=0.5`, `SEED=3` (`generate_heteropolymer_network.py:45-51`).
- Local `raf/` package (`compute_maxRAF`, `cot_bridge.closure`), then `pyCOT.analysis.decomposition`, with `SO_MAX_ORDER=60`.

**Assessment:**
- One random seed and polymer length ≤ 3.
- The cleavage probability was chosen after observing that full reversibility removes the effect (documented in the script).
- e_coli_core has no catalysts, so its maxRAF is empty by definition and that comparison is one-sided.
- The README points to a moved module path.

### 4.5 Generative_Structure_Orgs

**Purpose:** scaling of ERC, synergy and complementarity counts across BioModels and BiGG. No README.

**Pipeline:** `Cleaned/config.py:33-37` (`MAX_REACTIONS=1000`, `MAX_ERCS=900`, `MIN_ERCS=4`). Uses the **legacy** `ERC_Hierarchy` / `ERC_Synergy_Complementarity` with `.pkl` caches.

**Assessment:**
- Synergy and complementarity statistics come from the legacy module with inverted containment (§3.5), so they are not comparable with `organizations/` results.
- The size filter keeps 11/108 BiGG models.
- The dataset label is blank in 428/438 rows of `synergy_stats.csv` [read, agent].
- `scripts/` and `Cleaned/` are near-duplicate trees.

### 4.6 XCEPT_Article_1 (Lake Chad conflict, the climate/conflict application)

**Claims:**
- Transition structure between four socio-political organizations (Tribe, Chief, State, Chief-State), with a dominant 4-cycle TrC → ChC → CSC → StP (`reports/section_cot_lp.tex`).
- README: multi-pronged interventions succeed 70–90% of the time vs 20–40% for single-focus ones [reported; no supporting output found; README references scripts such as `run_analysis.py` that do not exist].

**Pipeline:**
- **The four organizations are hard-coded species sets** (`script_complex_COT_ESMO.py:75-79`), not computed with `compute_organizations` [read]. Characteristic peace/conflict states are hand-chosen.
- Elementary modes (ESMOs) come from randomized LP vertex sampling: `N_RAND_1=2000`, `N_EXPLORE=40`, `N_RAND_3=600`, `PATIENCE=200`, `SEED=42`, `EPS_SIGN=1e-3`, `NEUTRAL_TOL=5.0`.
- Mode counts are used as Markov transition weights.
- ODE/stochastic simulations use hand-set rate constants (`script_basic_model.py:94-106`, `T_END=400`, seed 42; `script_section4.py` `N_SEEDS=10`).

**Assessment:**
- Transition probabilities rest on a **non-exhaustive, single-seed** sample of modes, with `Pr_feasible = 1.0` everywhere and one state pair per edge.
- No sensitivity to the sampling budget or seed is reported.
- Whether the four hard-coded sets are actually organizations of the model was not checked against `compute_organizations` in the repo.
- Rate constants are not calibrated to data. ODE (Michaelis–Menten-like) and stochastic (min-based) kinetics differ.

### 4.7 XCEPT_Article_2 (Adaptive Peace Simulation, Lake Chad)

**Pipeline:** `aps_lib_stochastic.py` hard-codes reactions and rate constants (l.337-394). Its claimed source file `Conflict_Model_Article_2_v1.txt` does not exist. Runs on `BudgetAllocationSimulator` (CLI defaults: boost 5, 15 periods, 10 runs, dt 0.05, seed 42).

**Assessment:**
- No COT computation (closure, organizations, self-maintenance) is performed, despite the self-maintenance language in the paper.
- Rate constants were tuned by eye (`script_stochastic_model_tuning.py`).
- Inherits the `1/climate` degradation bug and the floor/clip conservation bias of the simulator (§2.3).
- The README budget table sums to 120–155%.

### 4.8 AMF_model / AMF (arbuscular mycorrhizal fungi trade)

**Pipeline:** `simulations.ode.simulation` (LSODA) on `data/Ecological_models/AMF.txt`.
- `t_span (0, 600)`, 601 steps; the last 20% is treated as steady state.
- Four scenarios via soil N/P and input rates (`AMF.py:353-385`).
- Monte Carlo ±10%, n = 100, seed 123.

**Assessment:**
- Scenario IV uses "calibrated" Ns = 92.64, Ps = 78.40, far outside the documented ranges (Ns 10–40, Ps 5–20), with no source [read, agent].
- Mass-action kinetics mixes % and g/kg units.
- `x0_dict` is mutated across scenarios.
- The Monte Carlo script points to a missing file.

---

## 5. Data catalogue (`data/`)

| Folder | Files | Provenance | Consumers |
|---|---|---|---|
| `Conflict_models/` | 8 `.txt` | Hand-built conflict/resource toy models. Some lines carry LaTeX labels (`#{eq:…}\\`), which `#`-stripping removes. | Root `scripts/` (stale paths); XCEPT keeps its own copies |
| `Ecological_models/` | `AMF.txt`, `MLM_Afidos/RN_MLM26_V02.txt`, `Marine_Ecosystem/Las_Cruces_251021.txt` | Hand-built | AMF_model, Fundamental Generators (visualizations), root scripts |
| `Examples_tests/` | 123 files (≈107 txt, 6 pkl) | Toys; Centler 2006 transcribed by script; FarmVariants (15); LUCA copies | Fundamental Generators, Decomposition, Generative |
| `LUCA/` | 1 xlsx (`42003_2021_1918_MOESM6_ESM.xlsx`; the file-name prefix 42003 is the *Communications Biology* supplement scheme, 2021 [read from filename; article identity not verified]) | Journal supplement | Unused directly |
| `biochemical_databases/` | 6513 files: `biomodels_all/` (SBML), `biomodels_all_txt/` (txt + ERC `.pkl` caches), `BioMD_*` category folders, `BiGG/` (108 txt + 14 pkl), `biomodels_interesting/`, `Other/` | BioModels SBML → txt; BiGG models (reversible reactions split into `_fwd`/`_rev`, biomass kept) | All biology projects |
| `ai_generated/` | 302 `.txt` (RandomAlife, Ns20, GlobalWarming0/1, EColi) | Generated | Root `scripts/analysis` |

The `.pkl` files are cached ERC computations. They become stale if the ERC code changes, and they are trusted without a version check [plausible].

---

## 6. Reproduction of the [measured] items

Environment: `pip install rustworkx scipy numpy pandas sympy networkx matplotlib pytest`.

1. `python -m pytest -q src/pyCOT/tests --ignore=src/pyCOT/tests/test_load_from_file.py` gives 124 passed, 2 failed.
2. For each of `data/Examples_tests/Farm.txt`, `data/biochemical_databases/BiGG/bigg_e_coli_core.txt` and `data/biochemical_databases/biomodels_all_txt/BIOMD0000000943.txt`:
   - Run `compute_organizations(rn, max_so_order=m, verify_organizations=False)` for m ∈ {2, 5, 10, 100}.
   - Run it again with `use_vertical_lift=False`.
   - Call `check_self_maintenance(species, rn, epsilon=ε)` for ε ∈ {1e-3, 1e-4, 1e-6, 1e-8, 1e-10} on up to 400 SOs.
   - Test `is_closed` / `is_ssm` on `species_mask | E0_mask`.

| Network | SOs at cap 2/5/10/100 | No lift | Not closed | Orgs at ε = 1e-3 / 1e-4 / 1e-6 / 1e-8 / 1e-10 | Elem. size with / without E0 |
|---|---|---|---|---|---|
| Farm.txt | 14/14/14/14 | 14 | 0 | 4/4/4/14/14 (of 14) | 5.33 / 4.33 |
| e_coli_core | 222/300/309/309 | 252 | 8 | 13/13/13/309/309 (of 309) | 15.13 / 2.13 |
| BIOMD0000000943 | 94/416/1224/1267 | 84 | 0 | 23/23/23/400/400 (of first 400) | 7.88 / 1.88 |

---

## 7. Recommendations (highest leverage first)
1. Report module sizes as |X \ E0| (or both). Re-state the size-invariance, fusion and hybridization claims on that basis.
2. Add an ε sensitivity sweep (≥ 1e-6 and ≥ 1e-5) and a post-solve residual check `min(S·v) ≥ −margin`, `min(v) ≥ ε`. Raise instead of returning False on solver exceptions.
3. Assert `closure(mask) == mask` and `is_ssm(mask)` for every SO returned. De-duplicate by closure.
4. Resolve the core vs `gen_v2` SO count discrepancy with an exhaustive oracle on small networks. Raise the default `max_so_order`, or warn when the cap is hit.
5. XCEPT-1: compute the organizations of the model with `compute_organizations` rather than hard-coding them. Report mode-sampling convergence across seeds and budgets.
6. Process analysis: fix the mode-label strings and make τ relative to the flux scale.
7. Fix `read_txt` for multi-reaction lines and the leading-decimal coefficient parse. Make `exact_names=True` the default for curated models.
