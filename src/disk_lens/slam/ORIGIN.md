# Local SLaM stage provenance

The release inventory records a copy on 25 September 2026 from the historical
fitting project's local SLaM directory, associated with PyAutoLens
`2025.1.18.7`. The original checkout revision was not recorded.

On 26 October 2026 these modules were compared with the official
[PyAutoLens Workspace tag `2025.1.18.7`](https://github.com/PyAutoLabs/autolens_workspace/tree/2025.1.18.7),
resolved from the locally available upstream Git history to commit
[`cb5a94f41273c0df26ad9bc269d25014594ebd3d`](https://github.com/PyAutoLabs/autolens_workspace/tree/cb5a94f41273c0df26ad9bc269d25014594ebd3d).
The tag object is `f0436d7e1bc2a95c251214c7fac1fea4bf923d94`.
This identifies a reproducible comparison baseline; it does not establish
that the historical copy was made from that exact checkout.

## Differences from that baseline

- `source_lp.py`: adds `n_live=800` to `run` and passes it to Nautilus,
  replacing the upstream hard-coded value of 200.
- `source_pix.py`: adds `n_live=800` to `run_1` and `n_live=1000` to `run_2`,
  replacing hard-coded values of 150 and 100 respectively. The second search
  remains `DynestyStatic`. Remaining differences are whitespace/trailing comma.
- `light_lp.py`: adds `n_live=1000` and passes it to Nautilus instead of the
  upstream hard-coded 150.
- `mass_light_dark.py`: adds `n_live=1500` instead of the upstream hard-coded
  250. The local version omits the upstream assignments that overwrite
  `dark.mass_at_200`, `dark.redshift_object` and `dark.redshift_source` inside
  the dark-profile setup. This preserves the locally supplied model choices.
- `slam_util.py`: adds `update_progress_csv_file` and `update_result_json_file`.
  It changes `output_fit_multi_png`, `output_result_to_fits`,
  `output_source_multi_png`, `plot_fit_png_row` and `plot_source_png_row`.
  The other four top-level functions have identical ASTs at this baseline.
- `__init__.py`: imports only the five required modules. Three unused stage
  variants remain in the external historical archive.

The live-point values above describe the module signature defaults, not every
historical production run. The pipeline passes the configured stage values;
saved run manifests remain the authority for historical 1500/2000 runs.

## File identity

SHA-256 hashes of the Python files at this audit:

| File | SHA-256 |
| --- | --- |
| `__init__.py` | `767e4eea0e716b15551d80058f1d34231e5255419fbec0c12e9b95ab12850017` |
| `source_lp.py` | `43db34f99381fe0ee3a8135f3cadabe07aa58b18ab518f0ad445fd8bb085488b` |
| `source_pix.py` | `b202a29897cc0a2d52728aafd3e8eb6b1dc81be296eb03ae0c29482e0a4886f4` |
| `light_lp.py` | `ee8b940de0aac6dae24813e6830a046456f5804ed665697ef820ce0e60628957` |
| `mass_light_dark.py` | `bfa6d3d1d4c71ada02f44913c5f41b6b367f0a9157b3cc47aaee11c1b78e6512` |
| `slam_util.py` | `6b0f3ad6966b6e739e20b51e749cc20089836695824ebe519fe6619904d08f2c` |

The audit's per-file comparisons and unified diffs are retained outside the
distributable source tree. No scientific stage implementation was changed by
the attribution audit.

## Redistribution status

The applicable permission/license for this historical workspace material
still requires confirmation; see [THIRD_PARTY.md](../../../THIRD_PARTY.md).
No upstream copyright notice has been invented or replaced with the license
of a separately installed library.
