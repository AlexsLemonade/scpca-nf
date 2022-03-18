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
- [ ] Update code and documentation with the latest version number:
  - [ ] [nextflow.config](https://github.com/AlexsLemonade/scpca-nf/blob/main/nextflow.config) 
  - [ ] [README.md](https://github.com/AlexsLemonade/scpca-nf/blob/main/README.md)
  - [ ] [external-data-instructions.md](https://github.com/AlexsLemonade/scpca-nf/blob/main/external-data-instructions.md)

### Creating a release
- [ ] On the [releases page](https://github.com/AlexsLemonade/scpca-nf/releases), choose `Draft a new release`.
- [ ] In `Choose a tag`, type a new release number using semantic versioning (vX.X.X) (you did update the title of this issue to match, right?), then click `Create a new tag: vX.X.X on publish`.
- [ ] Write a description of the major changes in this release. You may want to start with the Auto-generated release notes to save time.
- [ ] Save a draft to return to later if testing isn't done yet, otherwise:
- [ ] Publish the release!
