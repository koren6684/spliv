NOT RELEASE READY

# Final release decision for SPLIV 0.2.0

Audit date: 2026-07-31.

The staged package is a locally validated release candidate: 227 tests pass,
the clean pkgdown site builds, and the exact source tarball passes
`R CMD check --as-cran` with 0 errors, 0 warnings, and 1 expected recent-update
NOTE. The release as a whole is not ready because the remaining blockers affect
the identity and reproducibility of the public scientific artifact.

## Blocking items

1. The remote lightweight tag `v0.2.0` already points to baseline commit
   `41993897f83f9388abd2cd13368c2e36ae91aefa`, rather than the current audited,
   uncommitted hardening tree. The maintainer must choose a new immutable
   version/tag strategy; this pass did not move or overwrite the tag.
2. `spliv-reproducibility/renv.lock` resolves `spliv` from that same older
   tagged commit. After the definitive package release exists, regenerate and
   validate the lockfile against its immutable commit.
3. The full 1,000-replicate simulation reference suite has not been rerun after
   the BPE scale correction. The quick deterministic pilot records a real BPE
   eligibility/fit-success change from `0.1/0.1/0.0` to `0.0/0.1/0.0`.
   Therefore the existing full BPE simulation reference table cannot be
   certified until the production run and comparison are completed.
4. The release-candidate changes are intentionally uncommitted and hosted
   multi-platform CI has not run on the definitive commit/tag.
5. R Journal submission still requires a manuscript-level related-software
   comparison, an end-to-end manuscript build, and archival of the exact
   package, lockfile, code, data-acquisition record, and article outputs.

## Completed scientific evidence

- Baseline IV coefficients and IID, HC1, and clustered standard errors match
  equivalent `fixest` models to `1e-9`; UCI endpoints match manual adjusted-
  outcome calculations to `1e-9`.
- BPE scaling is invariant to multiplying the instrument by 100 or 0.01.
- The full Koren workflow completed on checksum-verified external data and the
  compact article outputs match their staged references byte-for-byte.
- The full Lelkes workflow completed on checksum-verified external data; its
  corrected standardized BPE decision remains ineligible.
- Clean-room package and reproducibility tests load the exact built package
  from an installed library, not a neighboring source checkout.

No repository was committed, pushed, tagged, released, or submitted to CRAN
during this audit.
