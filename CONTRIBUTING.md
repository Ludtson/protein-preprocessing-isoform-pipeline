# Contributing

This project is maintained by the Casola Lab (Ecology & Conservation Biology
Program, Texas A&M University) as part of the LINGUA comparative genomics
framework. Contributions are welcome, particularly around species-annotation
edge cases (non-standard GFF dialects, additional source databases) and
performance improvements.

## Getting set up

The pipeline requires bash and Python 3.7 or later. There are no third-party
Python dependencies; everything is standard library.

```bash
git clone https://github.com/Ludtson/protein-preprocessing-isoform-pipeline.git
cd protein-preprocessing-isoform-pipeline
```

## Running the test dataset

The repository includes a small three-species dataset under `example_data/`
with a known-good output under `example_data/test_run/`. This is the primary
way to check that a change hasn't altered pipeline behavior:

```bash
bash scripts/longest_protein_pipeline.sh \
  -g example_data/gff \
  -f example_data/protein \
  -o test_output

diff test_output/master_summary_report.tsv \
  example_data/test_run/master_summary_report.tsv
```

No output from `diff` means the run reproduced the expected result. If a
change intentionally alters the summary report's columns or values, update
the fixture under `example_data/test_run/` in the same commit and explain the
change in the pull request description.

## Code style

- Bash scripts should keep `set -euo pipefail` at the top and fail loudly on
  unexpected conditions rather than silently continuing.
- Python scripts should stick to the standard library unless there's a strong
  reason to add a dependency; this keeps the pipeline easy to deploy on HPC
  clusters without a package manager.
- Match the logging conventions already in use (`[INFO]`, `[WARNING]`,
  `[ERROR]` with a timestamp) rather than introducing a new style.

## Submitting changes

1. Open an issue first for anything beyond a small fix, so the approach can
   be discussed before implementation.
2. Keep pull requests focused on one change. Unrelated cleanup makes review
   harder and should go in a separate PR.
3. Include the reproducibility check output (or an explanation of why the
   fixture needed updating) in the PR description.
