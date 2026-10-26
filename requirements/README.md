# Unified mock, truth and fitting environment

New SPH mock rendering, mock truth and five-stage fitting all use the same
Python 3.11 environment with PyAutoLens, PyAutoGalaxy, PyAutoArray, PyAutoFit
and PyAutoConf **2025.1.18.7**, plus lenstronomy **1.13.2**.

Install the resolver-compatible bootstrap, then restore the direct production
versions in the same isolated environment:

```bash
python -m pip install -c requirements/bootstrap.txt '.[mock,lensing]'
python -m pip install --no-deps -r requirements/production.txt
```

`bootstrap.txt` fixes NumPy 1.26.4 and SciPy 1.14.0 so that the published PyAuto
dependency metadata can resolve. `production.txt` fixes the direct versions
used in the fitting environment, including NumPy 2.1.3 and SciPy 1.16.0, and
adds the same pinned lenstronomy version for mock rendering. AutoConf declares
`numpy<=2.0.1`, anesthetic 2.8.14 declares `numpy<2`, and AutoFit/AutoArray
declare `scipy<=1.14.0`; the second command deliberately reconstructs the
observed direct versions. Consequently `pip check` reports these known
upstream metadata conflicts after that step. This is not a claim that upstream
supports that NumPy/SciPy combination. See `docs/INSTALLATION.md` for validation
scope and the numerical comparison against the saved historical mock data.

Nautilus 1.0.4, Dynesty 2.1.4 and Colossus 1.3.1 are explicitly fixed in both
constraint files. They are direct-version snapshots, not cross-platform
transitive lockfiles. The base package requires NumPy below 2.2 because Astropy
6.1 imports APIs removed in newer NumPy versions.

For compatibility, `mock.txt` and `fitting.txt` include `production.txt`, and
`fitting-bootstrap.txt` includes `bootstrap.txt`. Apply the runtime aliases
with `--no-deps` after bootstrapping, as shown above.

`legacy-mock-2024.txt` preserves the old mock environment for historical audits
in a separate environment; new release commands use the unified 2025 stack.
Saved 2024 mock manifests retain their original versions and provenance.

TNG extraction additionally needs the separately installed `illustris_python`
reader described in `docs/INSTALLATION.md`; selection-radius analysis uses the
`selection` extra. Installation, CLI and resource checks are distinct from a
complete particle-to-posterior refit.
