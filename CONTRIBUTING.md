# Contributing to `scpca-nf`

We welcome contributions to the `scpca-nf` workflow, including usage reports and suggestions from users.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Issues](#issues)
- [Pull requests and branch structure](#pull-requests-and-branch-structure)
  - [Updating `nextflow_schema.json`](#updating-nextflow_schemajson)
  - [Continuous integration in pull requests](#continuous-integration-in-pull-requests)
- [Stub workflows](#stub-workflows)
- [Code style](#code-style)
  - [Nextflow](#nextflow)
    - [Using groovy variables](#using-groovy-variables)
    - [Declaring groovy variables](#declaring-groovy-variables)
    - [Channel transformations](#channel-transformations)
    - [Explicit variables in closures](#explicit-variables-in-closures)
    - [Spacing](#spacing)
  - [R and R Markdown](#r-and-r-markdown)
  - [Python](#python)
- [Environment management with `pixi`](#environment-management-with-pixi)
- [Pre-commit hooks](#pre-commit-hooks)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

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
This file can most easily be updated using the [`nf-core/tools` package](https://nf-co.re/tools), which is included in the [`pixi` environment](#environment-management-with-pixi) (or it can be installed separately).
Run `nf-core pipeline schema build` in the `scpca-nf` directory to update the schema file to match the current config file.
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

We additionally adopt the Nextflow style conventions described below.

#### Using groovy variables

We prefer to refer to variables directly as `variable` and not `${variable}` unless string interpolation is used.
In addition, when indexing a value from a map/dictionary, use `.` and not `[""]`.
Consider for example a variable `meta.unique_id`:

```groovy
// do this:
meta.unique_id

// this is an acceptable scenario with string interpolation:
outfile = "${meta.unique_id}_output.txt"

// don't do this:
meta["unique_id"]
${meta.unique_id}
```


#### Declaring groovy variables 

Variables in Groovy are set as global by default, which can have some unexpected consequences.
To avoid this, any variables declared in functions or closures (such as `.map{ }` statements) should be prefixed with `def`, which defines them as locally-scoped variables.

```groovy
// do this:
good.map{
  def my_var = it;
  my_var + 1
}

// don't do this:
bad.map{
  my_var = it;
  my_var + 1
}
```

Using `def` also improves code readability in closures:

```groovy
// do this:
// note that we preface unused variables with underscores
good.map{ meta, thing1, _thing2 -> 
  // define the long item first for clarity in the returned list
  def long_file = file(meta.long_file_name_in_meta, checkIfExists: true)
  [meta, thing1, long_file]
}


// don't do this:
bad.map{ meta, thing1, _thing2 -> 
  [meta, thing1, file(meta.long_file_name_in_meta, checkIfExists: true)]
}
```

Note that we do not use `def` in script blocks in Nextflow processes.

#### Channel transformations

We return lists, not tuples, from channel transformations.
In multi-line channel transformations, the `->` should appear at the _end_ of a line.

```groovy
// do this:
good.map{ meta, thing1, _thing2 -> 
  def long_file = file(meta.long_file_name_in_meta, checkIfExists: true)
  [meta, thing1, long_file]
}


// don't do this:
bad.map{ meta, thing1 -> [meta, thing1, file(meta.long_file_name_in_meta, checkIfExists: true)] }

// also don't do this:
bad.map{ meta, thing1 -> tuple(
    meta, 
    thing1, 
    file(meta.long_file_name_in_meta, checkIfExists: true)
  ) 
}
```



#### Explicit variables in closures

We use explicit variables in closures for all `.map{ }` statements and all other multi-line closures, but not for small one-line closures.

```groovy
// do this:
good.map{ it -> 
  def id = it[0].unique_id
  [it, id]
}

// don't do this:
bad.map{ 
  def id = it[0].unique_id
  [it, id]
}
```


In `.map{ }` closures, the code base follows the convention that `it[0]` will refer to the `meta` map/dictionary that carries library/sample metadata through the workflow.
If specifying just `it` as the explicit closure variable would require additional indexing during the channel transformation, we prefer to be even more explicit:

```groovy
// do this:
good.map{ meta, _unfiltered, _filtered, processed -> 
  [meta, processed]
}

// don't do this:
bad.map{ it -> 
  [it[0], it[3]]
}
```

#### Spacing

We generally adopt these spacing conventions:

```groovy
// For closures, use spaces surrounding braces {}:
good.map{ it -> [it.unique_id, it] } // do this
bad.map{it -> [it.unique_id, it]} // not this

// Parentheses () and brackets [] do not need spaces, e.g.:
good.join(other_thing) // do this
bad.join( other_thing ) // not this

// if/else spacing as follows:
if (condition) { // do this
if(condition){ // not this, or any other spacing variant
```


### R and R Markdown

For R code, we try to follow [`tidyverse` style conventions](https://style.tidyverse.org), and encourage the use of the [`styler`](https://styler.r-lib.org/) package to ensure that code is formatted consistently.


### Python

For python code, we encourage the use of the [`black` code formatter](https://black.readthedocs.io/en/stable/) to ensure consistent formatting.
The `black` package can be installed with `pip install black`, and can be run on a file with `black <filename>`.
If you use Visual Studio Code, you can install the [`black` extension](https://marketplace.visualstudio.com/items?itemName=ms-python.black-formatter).

As an alternative to `black`, you can also use [`ruff`](https://docs.astral.sh/ruff/) for formatting and linting; it follows the same code style conventions.
`ruff` is also available as a [VS Code extension](https://marketplace.visualstudio.com/items?itemName=charliermarsh.ruff).


## Environment management with `pixi`

To ease installation of dependencies and help with reproducibility, we have included a [pixi](https://pixi.sh) configuration file (`pixi.toml`) and lockfile (`pixi.lock`) in the root of the repository.
`pixi` is a package manager that takes advantage of the conda ecosystem, but faster and with some nice extra features like automatically keeping the environment in sync with the lockfile and the ability to create multiple environments for a project.
The `pixi.toml` file contains a list of the high-level dependencies needed to run the workflow, and the `pixi.lock` file contains the exact versions of those dependencies that were used when the lockfile was created.

The included default environment includes `nextflow`, `nf-core`, `pre-commit`, and their dependencies.

To use this, you will need to have installed the `pixi` package manager, which is most easily done with the following command ([other installation methods](https://pixi.sh/latest/advanced/installation/)):

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

To use the default environment, you can then either prefix your command with `pixi run` (e.g. `pixi run nextflow`), or you can launch a shell with the environment activated by running `pixi shell`.
Either method will automatically update the environment to match the package versions defined the `pixi.toml` and/or `pixi.lock` files.


## Pre-commit hooks

For convenience, we have included a set of [pre-commit hooks](https://pre-commit.com/) that can be used to automatically format code according to the above specifications, as well as to spellcheck and check for other common errors.

Once you have [installed `pixi`](#environment-management-with-pixi), you can run `pixi run pre-commit install` from the `scpca-nf` directory to use these hooks.
Alternatively, install the `pre-commit` package according to your favorite method (`pip`, `conda`, etc.) and then run `pre-commit install` in the `scpca-nf` directory.
This will install the hooks in the `.git/hooks` directory, and they will be run automatically when you commit changes.
If any of the hooks fail, the commit will be aborted, and you will need to fix the errors and re-commit.

Notably, the spellcheck hook will report spelling errors, but will also add any words it finds to the dictionary file.
This is convenient for many cases (where the word is real but unknown), but be sure to remove truly misspelled words from the dictionary file before committing, or they will not be caught in the future!
