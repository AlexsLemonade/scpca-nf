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

The [`test/stub-metadata.tsv`](test/stub-metadata.tsv) file is used to define input libraries that will be used for testing.
Any additions to the overall workflow that will allow processing of a new library type should be added into `test/stub-metadata.tsv`, along with the appropriate input files (usually empty files with the expected names) for that library type in the `test/runs/` directory.
If a new reference type is needed, that should be defined in the [`test/stub-refs.json`](test/stub-refs.json) file.

