---
name: Release checklist
about: Prepare for a new release version of scpca-nf
title: Prepare for scpca-nf release vX.X.X
labels: release
assignees: ''

---

## Steps for a new release of `scpca-nf`

### Preparing for the release

- [ ] Are all of the issues planned for this release resolved? If there are any issues that are unresolved, mark this issue as blocked by those on ZenHub.
- [ ] Update code and documentation with the latest version number in the `development` branch:
  - [ ] [nextflow.config](https://github.com/AlexsLemonade/scpca-nf/blob/main/nextflow.config)
  - [ ] [internal-instructions.md](https://github.com/AlexsLemonade/scpca-nf/blob/main/internal-instructions.md)
  - [ ] [external-instructions.md](https://github.com/AlexsLemonade/scpca-nf/blob/main/external-instructions.md)
- [ ] Test that the workflow is in good working order with `nextflow run alexslemonade/scpca-nf -latest -r development`
- [ ] File a PR from the `development` branch to the `main` branch. This should include all of the changes that will be associated with the next release.

### Creating a release
- [ ] On the [releases page](https://github.com/AlexsLemonade/scpca-nf/releases), choose `Draft a new release`.
- [ ] In `Choose a tag`, type a new release number using semantic versioning (vX.X.X) (you did update the title of this issue to match, right?), then click `Create a new tag: vX.X.X on publish`.
- [ ] Write a description of the major changes in this release. You may want to start with the auto-generated release notes to save time.
- [ ] Optional: If not all issues have been addressed, save a draft to return to later.
- [ ] Publish the release!
