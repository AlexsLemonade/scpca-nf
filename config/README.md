This folder contains additional configuration files for running the `scpca-nf` workflow, as described below.
To use a configuration file, it must be referenced in the main configuration file (i.e., [`nextflow.config`](../nextflow.config)) using `includeConfig`.


- `containers.config`: Docker containers used by `scpca-nf` workflow processes
- `process_base.config`: default resource settings for `scpca-nf` workflow processes
- `profile_awsbatch.config`: legacy AWS Batch configuration for Data Lab users
- `profile_awsbatch_auto.config`: AWS Batch configuration with auto-scaling disks for Data Lab users
- `profile_ccdl.config`: file locations for Data Lab users only, as described in [`internal-instructions.md`](../internal-instructions.md#running-scpca-nf-as-a-data-lab-staff-member)
- `profile_stub.config`: file locations used for testing `scpca-nf` workflow logic, as described in [`internal-instructions.md`](../internal-instructions.md#testing-the-workflow)
- `profile_example.config`: file locations for re-generating example `scpca-nf` workflow output, as described in [`internal-instructions.md`](../internal-instructions.md#processing-example-data)
