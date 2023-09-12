# Contributing to `scpca-nf`

We welcome contributions to the `scpca-nf` workflow, including usage reports and suggestions from users.

## Issues

Most updates to the `scpca-nf` workflow will begin with an issue that describes a change to be made and the reasoning behind that change.
This provides an opportunity for discussion before implementation of any changes.

## Pull requests and branch structure

The `main` branch holds the current release version of the `scpca-nf` workflow.

New features and other workflow updates that are ready in advance of a new release are found in the `development` branch.

Contributions and updates to the `scpca-nf` repository operate on a pull request model.
Changes will typically be made in a new branch that is created from the `development` branch, followed by a pull request back to the `development` branch.
All pull requests must be reviewed before merging to `development`.
To allow for efficient review, please include in any pull request a concise and clear explanation of the changes you have made and the issues addressed.

When the changes in `development` merit a new release, a pull request will be filed to merge the current version of the `development` branch into `main`, followed by tagging a release on the `main` branch.

## Stub workflows

All Nextflow processes should include a [`stub` block](https://www.nextflow.io/docs/latest/process.html#stub) with a minimal script that can be run quickly to produce files in the expected output locations.
At this stage this is purely used to allow for testing of the main workflow logic rather than the internal logic of each process.

The [`test/stub-run-metadata.tsv`](test/stub-run-metadata.tsv) file is used to define input libraries that will be used for testing.
Each sample ID in the `stub-run-metadata.tsv` file should have a corresponding entry in [`test/stub-sample-metadata.tsv`](test/stub-sample-metadata.tsv) (which can be filled with `NA` values).
Any additions to the overall workflow that will allow processing of a new library type should include adding new example data.
This will involve adding rows to `test/stub-run-metadata.tsv` and `test/stub-sample-metadata.tsv`, along with the appropriate input files (usually empty files with the expected names) for each library in the `test/runs/` directory.
If a new reference type is needed, that should be defined in the [`test/stub-refs.json`](test/stub-refs.json) file.

## Code style

While there is not necessarily an established code style for nextflow code, we try to keep code neat and readable.
Line length should generally be kept under 100 characters, and indentation should be consistent.

For R code, we try to follow [`tidyverse` style conventions](https://style.tidyverse.org), and encourage the use of the [`styler`](https://styler.r-lib.org/) package to ensure that code is formatted consistently.

For python code, we encourage the use of the [`black` code formatter](https://black.readthedocs.io/en/stable/) to ensure consistent formatting.
The `black` package can be installed with `pip install black`, and can be run on a file with `black <filename>`.
Alternatively, if you use [Visual Studio Code](https://code.visualstudio.com), you can install the [`black` extension](https://marketplace.visualstudio.com/items?itemName=ms-python.black-formatter).
