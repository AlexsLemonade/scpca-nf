This folder contains additional configuration files for running the `scpca-nf` workflow.
To use a configuration file, it must be referenced in the main configuration file (i.e., [`nextflow.config`](../nextflow.config)) using `includeConfig`.

| Config file | Description |
|-----------|-------------|---|
| `containers.config` | Specifies Docker containers to use when running `scpca-nf` | Yes |
| `process_base.config` | Specifies resource settings for the `scpca-nf` workflow | Yes |
| `profile_awsbatch.config`| Specifies AWS resources for the `scpca-nf` workflow. | Yes |
| `profile_awsbatch_auto.config`| Specifies AWS resources for the `scpca-nf` workflow. | Yes |
| `profile_ccdl.config` | Used internally by the Data Lab for the `scpca-nf` workflow, as described in [`internal-instructions.md`](../internal-instructions.md#running-scpca-nf-as-a-data-lab-staff-member). | No |
| `profile_stub.config` | Used for testing the `scpca-nf` workflow logic, as described in [`internal-instructions.md`](../internal-instructions.md#testing-the-workflow). |  No |
| `profile_example.config` | Used to re-generate example `scpca-nf` workflow output, as described in [`internal-instructions.md`](../internal-instructions.md#processing-example-data). | No |
