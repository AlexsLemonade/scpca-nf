This folder contains additional configuration files for running the `scpca-nf` workflow.
To use a configuration file, it must be referenced in the main configuration file (i.e., [`nextflow.config`](../nextflow.config)) using `includeConfig`.


- `containers.config` specifies Docker containers to use when running the `scpca-nf` workflow
- `process_base.config` specifies resource settings for the `scpca-nf` workflow
- `profile_awsbatch.config` specifies AWS resources for the `scpca-nf` workflow
- `profile_awsbatch_auto.config` specifies AWS resources for the `scpca-nf` workflow
- `profile_ccdl.config` is used internally by the Data Lab for the `scpca-nf` workflow, as described in [`internal-instructions.md`](../internal-instructions.md#running-scpca-nf-as-a-data-lab-staff-member)
- `profile_stub.config` is used for testing the `scpca-nf` workflow logic, as described in [`internal-instructions.md`](../internal-instructions.md#testing-the-workflow)
- `profile_example.config` is used to re-generate example `scpca-nf` workflow output, as described in [`internal-instructions.md`](../internal-instructions.md#processing-example-data)
