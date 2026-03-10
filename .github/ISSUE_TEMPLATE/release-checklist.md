---
name: Release checklist
about: Prepare for a new release version of scpca-nf
title: Prepare for scpca-nf release `vX.X.X`
labels: release
assignees: ''

---

## Steps for a new release of `scpca-nf`

### Preparing for the release

- [ ] Are all of the issues planned for this release resolved? If there are any issues that are unresolved, mark this issue as blocked by those on ZenHub.
- [ ] Update code and documentation with the latest version number in the `development` branch:
  - [ ] [`nextflow.config`](https://github.com/AlexsLemonade/scpca-nf/blob/main/nextflow.config)
  - [ ] [`internal-instructions.md`](https://github.com/AlexsLemonade/scpca-nf/blob/main/internal-instructions.md)
  - [ ] [`external-instructions.md`](https://github.com/AlexsLemonade/scpca-nf/blob/main/external-instructions.md)
  - [ ] [`sb_nextflow_schema.yaml`](https://github.com/AlexsLemonade/scpca-nf/blob/main/cavatica/sb_nextflow_schema.yaml) and [`sb_doc.md`](https://github.com/AlexsLemonade/scpca-nf/blob/main/cavatica/sb_doc.md)
- [ ] Update the tag for `scpca-nf` containers in `config/containers.config` and `nextflow-schema.json`
  _Note:_ Any processes that use those containers will not work until the new tag has been created.
- [ ] Check that the `nextflow-schema.json` is up to date with `pixi run nf-core pipelines schema build`
- [ ] Check that [`sb_nextflow_schema.yaml`](https://github.com/AlexsLemonade/scpca-nf/blob/main/cavatica/sb_nextflow_schema.yaml) entries are up to date (sadly, no automated tool for this)
- [ ] Test that the workflow is in good working order with `nextflow run alexslemonade/scpca-nf -latest -r development`
- [ ] File a PR from the `development` branch to the `main` branch. This should include all of the changes that will be associated with the next release.
- [ ] (Optional) Generate new example `scpca-nf` output files.
If updating the example output is not necessary for this release, check these boxes off for free.
  - [ ] [Re-process the example data](https://github.com/AlexsLemonade/scpca-nf/blob/main/internal-instructions.md#processing-example-data) through the `scpca-nf` workflow and ensure it looks correct.
  - [ ] Compress the example output in `scpca_out` to create `scpca_out.zip`, as described in the instructions, and ensure the file is set to public and read-only.


### Creating a release
- [ ] On the [releases page](https://github.com/AlexsLemonade/scpca-nf/releases), choose `Draft a new release`.
- [ ] In `Choose a tag`, type a new release number using semantic versioning (`vX.X.X`) (you did update the title of this issue to match, right?), then click `Create a new tag: vX.X.X on publish`.
- [ ] Write a description of the major changes in this release. You may want to start with the auto-generated release notes to save time.
- [ ] Optional: If not all issues have been addressed, save a draft to return to later.
- [ ] Publish the release!
