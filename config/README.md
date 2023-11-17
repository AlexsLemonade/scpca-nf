This folder contains additional configuration files for running the `scpca-nf` workflow.
Some of these files are called from `nextflow.config`, and others can be optionally specified at run time.


| Config file | Description | Included from `nextflow.config`? |
|-----------|-------------|---|
| `containers.config` | Specifies Docker containers to use when running `scpca-nf` | Yes |
| `process_base.config` | Specifies resource settings for the `scpca-nf` workflow | Yes |
| `profile_awsbatch.config`| Specifies AWS resources for the `scpca-nf` workflow. | Yes |
| `profile_awsbatch_auto.config`| Specifies AWS resources for the `scpca-nf` workflow. | Yes |
| `profile_ccdl.config` | Used internally by the Data Lab for the `scpca-nf` workflow, as described in [`internal-instructions.md`](../internal-instructions.md#running-scpca-nf-as-a-data-lab-staff-member). | No |
| `profile_stub.config` | Used for testing the `scpca-nf` workflow logic, as described in [`internal-instructions.md`](../internal-instructions.md#testing-the-workflow). |  No |
| `profile_example.config` | Used to re-generate example `scpca-nf` workflow output, as described in [`internal-instructions.md`](../internal-instructions.md#processing-example-data). | No |
