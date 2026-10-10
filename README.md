# UDCT v7.12 Preregistration: Jeans Sign Transitions in Hydrus I

**Author:** Won Shik Paik, Independent Researcher, Auckland, New Zealand.

**Protocol prepared:** 10 October 2026 (Auckland). **Links added:** 11 October 2026 (Auckland).

## Links and registration status

- GitHub: https://github.com/wspaik5/UDCT_v7.12_Preregistration_Jeans_Transitions
- Zenodo DOI (reserved): https://doi.org/10.5281/zenodo.23280972
- v7.11 preregistration: https://doi.org/10.5281/zenodo.23244522
- v7.11 results: https://doi.org/10.5281/zenodo.23263099

At package preparation, the v7.12 DOI is reserved in a Zenodo draft; publication has not been confirmed. A DOI reservation is not public preregistration. Publish the complete protocol and source package before executing new galaxy cases, and retain the actual UTC publication timestamp. This README's status describes package preparation; consult the public deposit for its subsequent publication status.

This is a prospective follow-up motivated by already-known v7.11 results. It is not a blind prediction of those earlier results. The package contains no v7.12 galaxy outcomes. Only synthetic validation has been executed during preparation.

## Scope

The fixed design evaluates the signed stationary spherical-average Jeans closure for Hydrus I with an unsoftened tapered cusp. It varies constant anisotropy beta from 0 to 0.5 in steps of 0.05 and screening length across 1, 1.2, 1.5, 2, 2.5 and 3 pc. Five screened laws and three QUMOND controls are prescribed in the PDF.

The fixed matrix contains 198 field evaluations, 2,178 pressure records and 363 primary N/P/I decisions. Pre-specified adaptive refinement examines adjacent raw sign-changing brackets only, with numerical-margin stopping. Consult the PDF for the complete rules; it is the authoritative scientific protocol.

N indicates a resolved negative common band under the stated closure and boundary conditions. P is a positive sampled second-moment candidate, not proof of equilibrium, a positive distribution function or stability. I indicates an inconclusive numerical classification. This study does not establish a universal Hessian threshold or exclude MOND generally. Core softening and local-u extensions are deferred.

## Files

Keep all four Python files together in one directory:

- `UDCT_v7_12_Preregistration_Jeans_Transition.py`: new plan and execution adapter.
- `UDCT_v7_11_Preregistration_Anisotropic_Jeans.py`: frozen signed Jeans solver and decisions.
- `UDCT_v7_10_Postregistration_run.py`: frozen gravity implementation.
- `UDCT_v7_06_Hessian_Term_Sign_Audit_reproduce.py`: frozen antecedent dependency.
- `UDCT_v7_12_Preregistration_Jeans_Transitions_Linked.pdf`: linked protocol.
- `requirements.txt`: pinned NumPy and SciPy dependencies.
- `SHA256SUMS.txt`: checksums of the other seven files, excluding itself and the ZIP.
- `README.md`: this guide.

Uploading only the new Python runner is insufficient to execute the study. Extract the complete ZIP or upload all its files individually to GitHub.

## Preparation and synthetic validation

Preparation was checked with Python 3.12, NumPy 2.3.5 and SciPy 1.17.0.

```bash
python -m pip install -r requirements.txt
python UDCT_v7_12_Preregistration_Jeans_Transition.py
python UDCT_v7_12_Preregistration_Jeans_Transition.py --selftest
```

The default command prints the plan without running galaxy calculations. `--selftest` checks frozen source hashes and synthetic analytic solutions, decisions and bracket behavior. It does not evaluate Hydrus I or the Galactic field. The runner's plan retains its original preparation-time status text; actual publication provenance is separately supplied at execution.

To verify the package on systems with `sha256sum`:

```bash
sha256sum -c SHA256SUMS.txt
```

## Execute only after public registration

After publishing the exact protocol and source package, verify the DOI deposit is publicly accessible. Set the following shell variable to the actual public publication timestamp in timezone-aware UTC, for example in the format `YYYY-MM-DDTHH:MM:SSZ`. Do not use a guessed time or the DOI reservation time.

```bash
REGISTRATION_PUBLISHED_UTC='REPLACE_WITH_ACTUAL_PUBLICATION_TIMESTAMP_UTC'
python UDCT_v7_12_Preregistration_Jeans_Transition.py \
  --run \
  --registration-url https://doi.org/10.5281/zenodo.23280972 \
  --registration-published-utc "$REGISTRATION_PUBLISHED_UTC" \
  --confirm-public-registration \
  --output UDCT_v7_12_results_attempt01.json.gz
```

The placeholder is intentional and must be replaced. The runner rejects an invalid timestamp. Publication is attested by the operator; the runner does not independently retrieve the deposit.

The run writes a compressed incremental journal and a compressed final result. Existing output or attempt files are not overwritten. Preserve failures, partial files and journals; use a new output name for each attempt. Automatic resume is not implemented. Resource interruptions permit only an explicitly partial report.

## Source integrity and deviations

The link update and README do not change the four scientific Python files or the experimental design. The new runner checks the three frozen dependencies against their registered hashes. The PDF records all four Python hashes. The checksum manifest locks the linked PDF, README and requirements as well.

Any correction affecting registered quantities requires a public timestamped amendment before the affected cases are rerun. Preserve both versions, explain deviations, and report unfavorable and unresolved outcomes alongside favorable cases.
