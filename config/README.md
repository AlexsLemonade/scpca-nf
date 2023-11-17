This folder contains additional configuration files for running the `scpca-nf` workflow.
To use a configuration file, it must be referenced in the main configuration file (i.e., [`nextflow.config`](../nextflow.config)) using `includeConfig`.


- `containers.config` specifies Docker containers used by `scpca-nf` workflow processes
- `process_base.config`: default resource settings for `scpca-nf` workflow processes
- `profile_awsbatch.config`: legacy AWS Batch configuration for Data Lab users
- `profile_awsbatch_auto.config`: AWS Batch configuration with auto-scaling disks for Data Lab users
- `profile_ccdl.config`: file locations used internally by the Data Lab, as described in [`internal-instructions.md`](../internal-instructions.md#running-scpca-nf-as-a-data-lab-staff-member)
- `profile_stub.config` is used for testing the `scpca-nf` workflow logic, as described in [`internal-instructions.md`](../internal-instructions.md#testing-the-workflow)
- `profile_example.config` is used to re-generate example `scpca-nf` workflow output, as described in [`internal-instructions.md`](../internal-instructions.md#processing-example-data)
