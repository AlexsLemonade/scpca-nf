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

### Updating `nextflow_schema.json`

Any changes that affect the Nextflow configuration files should be reflected in the [`nextflow_schema.json` file](https://nf-co.re/tools#pipeline-schema).
This file can most easily be updated using the [`nf-core/tools` package](https://nf-co.re/tools), which can be installed with `conda install nf-core`.
Then run `nf-core schema build` in the `scpca-nf` directory to update the schema file to match the current config file.
You can also use the web editor that it launches to further customize the schema file.
Note that you may get warnings about any config parameters that include `${projectDir}`, as the build tool wants those to be replaced with absolute paths.
This is not necessary for the schema to be valid, so please keep those paths with the `${projectDir}` variable (enter `n` at the prompt).
### Continuous integration in pull requests

There are several automatic checks performed by GitHub Actions in all pull requests filed to `main` or `development`:

- [Check Nextflow config](.github/workflows/nextflow-config-check.yaml): This workflow ensures that there are no syntax errors in the Nextflow configuration files and that the `nextflow_schema.json` file is up to date. This check is required to pass before pull requests can be merged.
- [Check Nextflow stub](.github/workflows/nextflow-stub-check.yaml): This workflow ensures that the [stub workflow](#stub-workflows) runs without errors. This check is required to pass before pull requests can be merged.
- [Spell check R Markdown and Markdown files](.github/workflows/spell-check.yml): This workflow ensures there are no spelling errors in R Markdown and Markdown files. This check is not required to pass before pull requests can be merged.

There is also one additional `pre-commit ci` workflow which runs all [pre-commit hooks as described in this section](#pre-commit-hooks), except for the spell check pre-commit hook.
Although highly recommended, it is not required that this workflow passes before pull requests can be merged.

## Stub workflows

All Nextflow processes should include a [`stub` block](https://www.nextflow.io/docs/latest/process.html#stub) with a minimal script that can be run quickly to produce files in the expected output locations.
At this stage this is purely used to allow for testing of the main workflow logic rather than the internal logic of each process.

The [`test/stub-run-metadata.tsv`](test/stub-run-metadata.tsv) file is used to define input libraries that will be used for testing.
Each `scpca_sample_id` value in the `stub-run-metadata.tsv` file should have a corresponding entry in [`test/stub-sample-metadata.tsv`](test/stub-sample-metadata.tsv) (which can be filled with `NA` values).
Any additions to the overall workflow that will allow processing of a new library type should include adding new example data.
This will involve adding rows to `test/stub-run-metadata.tsv` and `test/stub-sample-metadata.tsv`, along with the appropriate input files (usually empty files with the expected names) for each library in the `test/runs/` directory.
If a new reference type is needed, that should be defined in the [`test/stub-refs.json`](test/stub-refs.json) file.

## Code style

### Nextflow

While there is not necessarily an established code style for Nextflow code, we try to keep code neat and readable.
Line length should generally be kept under 100 characters, and indentation should be set at 2 spaces per tab.
This can be set in [Visual Studio Code](https://code.visualstudio.com) with the following entries in `settings.json`:

```json
"editor.rulers": [100],
"[nextflow]": {
    "editor.tabSize": 2,
}
```

#### A note on variables in Nextflow/Groovy

Variables in Groovy are set as global by default, which can have some unexpected consequences.
To avoid this, any variables declared in functions or closures (such as `.map{}` statements) should be prefixed with `def`, which defines them as locally-scoped variables.

```groovy
// don't do this:
bad.map{
  my_var = it;
  my_var + 1
}

// do this instead:
good.map{
  def my_var = it;
  my_var + 1
}
```

### R and R Markdown

For R code, we try to follow [`tidyverse` style conventions](https://style.tidyverse.org), and encourage the use of the [`styler`](https://styler.r-lib.org/) package to ensure that code is formatted consistently.

### Python

For python code, we encourage the use of the [`black` code formatter](https://black.readthedocs.io/en/stable/) to ensure consistent formatting.
The `black` package can be installed with `pip install black`, and can be run on a file with `black <filename>`.
If you use Visual Studio Code, you can install the [`black` extension](https://marketplace.visualstudio.com/items?itemName=ms-python.black-formatter).

As an alternative to `black`, you can also use [`ruff`](https://docs.astral.sh/ruff/) for formatting and linting; it follows the same code style conventions.
`ruff` is also available as a [VS Code extension](https://marketplace.visualstudio.com/items?itemName=charliermarsh.ruff).

## Pre-commit hooks

For convenience, we have included a set of [pre-commit hooks](https://pre-commit.com/) that can be used to automatically format code according to the above specifications, as well as to spellcheck and check for other common errors.

To use these hooks, install the `pre-commit` package according to your favorite method (`pip install pre-commit` or `conda install pre-commit`), then run `pre-commit install` in the `scpca-nf` directory.
This will install the hooks in the `.git/hooks` directory, and they will be run automatically when you commit changes.
If any of the hooks fail, the commit will be aborted, and you will need to fix the errors and re-commit.

Notably, the spellcheck hook will report spelling errors, but will also add any words it finds to the dictionary file.
This is convenient for many cases (where the word is real but unknown), but be sure to remove truly misspelled words from the dictionary file before committing, or they will not be caught in the future!
