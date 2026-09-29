# research/external — third-party checkouts (not committed)

This directory holds local clones of the open-source NRW / T-R extraction codes
reviewed in `docs/github_survey.md`. **Its contents are gitignored and are never
committed or redistributed** (see `.gitignore`). Only this README is tracked.

Why:

* Several surveyed repositories have **no licence file**. Without a licence the
  default is "all rights reserved": we may read and run the code locally for
  evaluation, but we may not copy or redistribute it. Nothing from these
  repositories is copied into `src/` or anywhere else in this repository.
* Even for MIT/BSD code, keeping a pinned *pointer* instead of a vendored copy
  keeps the licensing surface of this project clean and the provenance exact.

How to reproduce the checkout:

```bash
scripts/fetch_external.sh            # all repositories, at the pinned commits
scripts/fetch_external.sh scatter-opt permittivitycalc   # a subset
FORCE=1 scripts/fetch_external.sh    # re-fetch even if present
```

The pins (name, URL, full commit hash) are in `research/external_pins.csv`.
The script checks out exactly that commit (shallow fetch by SHA) and verifies
`git rev-parse HEAD` afterwards.

`research/external/.build/` is created by `research/adapters/musho_cpp.py`
(compiled C++ binary); it is also ignored.

Runtime requirements of the adapters (`research/adapters/`):

* Python adapters: the project venv plus `uncertainties lmfit emcee corner
  seaborn tqdm` (permittivitycalc) and CPU `torch` (S-Param-Inversion).
* Octave adapters: GNU Octave ≥ 8 (`apt-get install octave`, ~8 min here).
* `musho_cpp`: `g++` with C++17.
