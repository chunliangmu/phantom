# AGENTS.md

Phantom: 3D SPH/MHD astrophysics code, written in Fortran 90 (`gfortran`/`ifort`/`ifx`/`aocc`). The whole codebase is Fortran; Python exists only for analysis scripts (`scripts/pyphantom`, docs build).

## Build

- The root `Makefile` is a thin wrapper that delegates everything to `build/Makefile` (`cd build; make ...`). Run make from the repo root.
- Both `SYSTEM` and `SETUP` are **required** on the command line, e.g.:
  `make SYSTEM=gfortran SETUP=star`
  - `SYSTEM` selects a compiler block in `build/Makefile_systems`; this machine always uses `gfortran` (flags in `build/Makefile_defaults_gfortran`).
  - `SETUP` selects a compile-time configuration block in `build/Makefile_setups` (e.g. `star`, `disc`, `sedov`, `test`, `testgrav`, ...). Unknown SETUP = hard error with a hint. Day-to-day work uses `SETUP=star` (sometimes `dustystar` or `radstar`) for binary-star simulations.
- Default target `make` builds the main `phantom` binary. Other binaries, each its own target (all need SETUP): `phantomsetup`, `phantomtest`, `phantommoddump` (alias `moddump`), `phantomanalysis` (alias `analysis`), `phantomevcompare`, ...
- Changing SETUP/SYSTEM/FFLAGS/FPPFLAGS automatically triggers `make clean` (tracked via `build/.make_last*` files). Full rebuilds on config change are normal, not a bug.
- Binaries land in `bin/` (gitignored via `bin/.gitignore` containing `*`); `build/*.o` and `*.mod` are gitignored, as is `build/phantom-version.h` (generated from version numbers at the top of `build/Makefile`).
- `OPENMP=yes` (OpenMP) and `DOUBLEPRECISION=yes` are defaults; `DEBUG=yes` adds debug flags. Extra switches via env var, e.g. `MPI=yes` (or `MPI=openmpi`).

### Config system (important)

- `.F90` sources are C-preprocessed (cpp); `.f90` are not. Preprocessor macros (`-DGRAVITY`, `-DMHD`, ...) come from the SETUP block in `build/Makefile_setups` and are printed by `checkparams` at the start of each build.
- **APR is NOT enabled by `SETUP=star`**: in the committed `Makefile_setups`, `APR=yes` (→ `-DAPR`) is set only by `SETUP=testapr`. To build the APR binary used for the binary-star runs, pass it explicitly: `make SYSTEM=gfortran SETUP=star APR=yes`. `build/.make_lastfppflags` shows what the last build used (`-DGRAVITY -DAPR` for the current run setup).
- Compile-time switches are `#ifdef` blocks in `src/main/config.F90` (module `dim`, e.g. `ISOTHERMAL`, `MAXPTMASS`, `APR`). Do NOT edit `build/config.F90` if it exists — the source of truth is `src/main/config.F90`; `build/` is just the object dir.
- Setup source files live in `src/setup/setup_*.f90`; a new setup must also get a `SETUP=` block with `KNOWN_SETUP=yes` in `build/Makefile_setups` or make will reject it. A few obsolete SETUPs (e.g. `mcfost`, `planets`, `binarydisc`, `warp`) are silently remapped to generic setups with a warning.

## Tests

- Build+run test suite: `make SETUP=test phantomtest && ./bin/phantomtest` (the `make test` target is the shortcut).
- Run a subset by passing an argument: `./bin/phantomtest derivs`, or `gravity`, `dust`, `gr`, `ptmass`, `apr`, `growth`, `nimhd`, `kdtree`, ... (valid names in `src/tests/testsuite.f90`).
- Tests are compiled conditionally by SETUP: e.g. `SETUP=testgrav` + `./bin/phantomtest gravity`; `SETUP=testgr` + `gr ptmass`; `SETUP=testdust` + `dust`; `SETUP=testapr` + `apr derivs`. Full mapping in `.github/workflows/test.yml`.
- The Fortran unit-test suite is `src/tests/`; standalone analysis tests are separate Python/Cython scripts in `scripts/` (e.g. `test_analysis_ce.py`).

## Run-directory workflow

- Simulations are NOT run from the repo root. Use a separate run directory (e.g. `2ia-mini-test/` — an untracked local AGB-star test case, SETUP=star).
- `scripts/writemake.sh` writes a `Makefile` that calls repo make with `RUNDIR="${PWD}"` and copies binaries (`phantom`, `phantom_version`) into the run dir.
- `2ia-mini-test/` is a **temporary** scratch folder for the user's test runs; it may be removed or renamed at any time (as may other untracked run dirs like `apr-diag-run/`). Create your own run folder (and then `writemake.sh star > Makefile` inside it) for any test work. `2ia-mini-test/nosplitmerge/` is a comparison run with APR splitting/merging disabled (same dumps) — useful as an energy-conservation control.
- Run `scripts/writemake.sh` ONLY inside a run folder — never overwrite the `Makefile` in the repo root. `SETUP` is taken from the first argument (e.g. `~/phantom/scripts/writemake.sh star > Makefile`); without an argument no SETUP is written.
- Input files (`star.in`, etc.), output dumps and such run dirs are **never committed**; keep them untracked.
- `phantomsetup` writes an initial dump like `star_00000.tmp` — the `.tmp` suffix means the file is "incomplete" (some column not yet filled in, though all info to derive it is present). The complete `star_00000` appears once the main `phantom` binary runs (`write_initial_dump` in `src/main/initial.F90` strips the suffix and deletes the `.tmp`).
- Runtime data tables (EOS files etc.) are found relative to the `PHANTOM_DIR` env var (else `./`); the repo's `data/` dir must be reachable for runs. In a run dir: `export SYSTEM=gfortran PHANTOM_DIR=/home/coder/phantom; make APR=yes && ./phantom star > run.log 2>&1`.
- **Memory**: this machine has only 8 GB RAM (8 GB swap). APR multiplies the allocation by 8 (`update_max_sizes` in `src/main/config.F90`), and the default `maxp_alloc = 5200000` (`src/main/config.F90`) needs ~42 GB and gets OOM-killed. Always run with `-maxp=100000` (single or double dash both work, e.g. `./phantom -maxp=100000 star`).
- **Timing**: ~3 min per 100 time units on this machine (4-core Celeron); a full 20000-unit run takes hours. For debugging, edit `tmax` (and `dtmax`, the dump interval) in `star.in` down to ~300-600 units.
- The binary **rewrites `star.in` at the end of each run** (`dumpfile`/`logfile` are bumped to the next numbers), so restore `dumpfile = star_00000` (and `logfile`) before restarting from the original dump. The `.ev` file name follows the `logfile` stem (star01.log → star01.ev), and its columns are labelled in the `#` header (time, ekin, etherm, emag, epot, etot, ...).
- Run simulations in the foreground with a long tool timeout — backgrounded (`nohup ... &`) processes die when the shell session ends. For energy work, also capture the log: the periodic summary prints `split since last` / `merged since last` counts and the per-level particle distribution (npart change per summary = splits − merges/2, because `nkilled` counts 2 per merged pair).

## Data analysis with Python

- The `.venv` at the repo root has the packages for analysis — activate it (`source .venv/bin/activate`) or call `.venv/bin/python` directly. **`scripts/pyphantom` is Python 2 only** and will not import under this venv — don't use it.
- **`sarracen` reads dumps directly**: `from sarracen import sarracen; sdf, sdf_sink = sarracen.read_phantom("star_00000")` (returns pandas DataFrames). Run it from the repo root so the venv python is used. DataFrames have per-particle `x,y,z,vx,vy,vz,u,h,poten,apr_level,mass,...` columns.
- Dump-file semantics that are easy to get wrong:
  - `poten(i)` is the **total** potential energy of particle i (the 0.5·m_i factor is already folded in, `force.F90` `epoti = 0.5*pmassi*(...)`). The `.ev` `epot` column is `Σ poten(i)` and is correct — do NOT re-weight by mass.
  - `mass` per particle = `massoftype(1)/2**(apr_level-1)` (APR). `apr_level` is stored in dumps and restored on restart; `init_apr` only rescales `massoftype` by `2**(apr_max-1)` when the dump has no levels.
  - `massoftype` and the header live in the tagged binary header; the arrays are Fortran sequential unformatted records (little-endian on this machine, real*8 by default), with each array written as a [16-char tag record][data record] pair and interleaved block-header records — a naive tag/data pairing will misalign.

## APR debugging notes (state as of this branch)

- The dominant energy-non-conservation bug found so far: **merging does not conserve potential energy** (still unfixed in this branch, HEAD = `ee3668fca`). `combine_two_particles` places the survivor at the pair midpoint; for pairs crossing a spherical APR boundary (the usual case), the midpoint is ~d²/8r closer to the centre, so the system becomes more bound and the energy is silently lost. Measured ~-1e-9 per merge (63% sink-gas + 37% self-gravity potential; |ΔE| ∝ M_sink·m·d²/r³) vs +3.6e-11 per split (children placed tangentially, slightly less bound). This accounts for most of the ~8e-4 relative etot drift over 20000 units in the star test. The existing `delta_ekin` → `vxyzu(4)` adjustment conserves E_kin+E_therm exactly per pair (and is mass-independent), so the fix is to also add the potential change to `vxyzu(4)` of the merged particle (`u -= ΔE_pot/m_new`, O(N) pair-sum) — a test of this reduced the drift ~10×. `src/main/apr.f90` is back to pristine after this investigation.
- Three smaller bugs found in the same audit were fixed by the user (commits `0dcc2e69a` "(apr) parallelization bug fix" and `f4d3ff77d` "(apr) bug fix") and are now in this branch (HEAD = `ee3668fca` == `origin/apr-flexbound9-relax-merge2`): (1) `merge_with_special_tree` now uses `aprmassoftype(igas, apr_level(mergelist(inodeparts(...))))` — it previously indexed `apr_level` with the special tree's `xyzh_merge` indices (1..nmerge) instead of real particle indices (check `src/main/apr.f90:589`); (2) `combine_two_particles` now does `treecache(5,keep) = treecache(5,keep)*2` instead of `- int(1,kind=1)`, which had set the tree-cache mass to ≈ -1 (`src/main/part.F90:1444`); (3) the merge loop's `relaxlist` update uses `!$omp critical` instead of `atomic capture` to avoid a read/write race with the follow-up scan (`src/main/apr.f90:637`). None of the three affected the measured energy drift.
- Still-open latent issues (verified by reading, all currently benign): `adjust_entropy` (utils_apr.f90) is entirely commented out; `integrate_geodesic` (split_dir=2 / GR) has a dt-shrinking loop that can hang when `iexternalforce=0`; `ref_dir=-1` maps the **innermost** region to level 1 = the **largest** particle mass (coarsest) — the centre of the star gets derefined, and merges happen at the inner boundary where |ΔE| ∝ 1/r³ is largest (the APR test suite uses `ref_dir=1`, centre finest; re-check the region↔level↔mass mapping before assuming the configuration does what it seems to).
- `poten` bookkeeping at events (halved on split, summed on merge) is consistent and transient — `poten` is recomputed from scratch in every force call, so diagnostics should ignore the stored values between dumps.
- Instrumentation gotchas: any new local used inside the merge loop's OpenMP region must be added to the `!$omp private(...)` lists or the build fails; keep lines ≤132 chars and declare all variables (gfortran `-Wall`/`-std=f2008` flag truncation and `implicit type`; `-Werror` only with `NOWARN=yes`, which the current build flags in `build/.make_lastfflags` do not have); after editing `apr.f90` rebuild with `make SYSTEM=gfortran SETUP=star APR=yes` (full rebuild on config change is normal).

## Docs

- `docs/` is a Sphinx site (readthedocs, config `.readthedocs.yml`); build requires `docs/requirements.txt`. The developer guide lives in `docs/developer-guide/` (styleguide, testing, setup, staging, bots, fortran, vscode).

## Git & environment

- `upstream/main` has an **old/early** version of adaptive particle refinement (APR) merged (its `src/main/apr.f90` is ~560 lines vs ~670 here). Since then the APR code has had considerable development — the latest versions live in the current working branch (`apr-flexbound9-relax-merge2`) and other APR branches (local `apr`, `apr-optimize`; remotes `bec/apr-*`, `origin/apr-*`). Bec (Bec Nealon, the lead APR developer) keeps her latest work on `bec/apr-conservation-bug2` (most recently updated of the `bec/apr-*` branches). Check `git log --oneline --all --grep=apr` / recent merges before assuming `upstream/main` represents current APR state.
- A main purpose of this working copy is debugging the new APR code (splitting/merging, conservation).
- Stay inside the working directory; avoid files outside it (especially `/tmp` — scratch files there may be wiped mid-session by reboot/hibernation). The `.venv` at the repo root is provided for Python analysis work.

## CI

- PRs must pass: build (all SETUPs × 4 compilers), test, mpi, mcfost, binary, growth, krome, coolingra workflows (`.github/workflows/`). Local approximation: `scripts/testbot.sh` and `scripts/buildbot.sh`.
