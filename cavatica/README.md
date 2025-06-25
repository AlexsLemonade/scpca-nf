# Cavatica port documentation

This document provides instructions for developers on how to port the `scpca-nf` workflow to Cavatica, a cloud-based platform for running bioinformatics workflows.

Instructions for Cavatica users can be found in the [Cavatica documentation](sb_doc.md) file.
Additional details about required input files and other options can also be found in the [external instructions](../external-instructions.md).

## Instructions for developers

### Software installation

The software used for porting to Cavatica is `sbpack`, available at https://github.com/rabix/sbpack

To ease installation and use, we have created a [`pixi`](https://pixi.sh) environment named `cavatica` that includes `sbpack` and its dependencies.
This is defined in the `pixi.toml` file in the root of the repository.
To use `sbpack`, within this environment, prefix your command with `pixi run -e cavatica` (e.g. `pixi run -e cavatica sbpack_nf -h`).
Alternatively, you can launch a shell with the environment activated by running `pixi shell -e cavatica`, and then run `sbpack` command without the prefix.

### Creating and updating the `sb_nextflow_schema.yaml` file

The `sb_nextflow_schema.yaml` file is a configuration file that defines the Nextflow schema for the Cavatica port.
This file can be generated or updated with the following command:

```bash
./cavatica-build.sh
```

This will use the existing file as a starting point in order to preserve any customizations.
If you want to start fresh, you can delete the existing `sb_nextflow_schema.yaml` file and run the command again.

### Customizing the `sb_nextflow_schema.yaml` file

The `sb_nextflow_schema.yaml` file created by the `sbpack_nf` command will not automatically include all of the customizations that we need for the Cavatica port.
Below is a (partial) list of some of the changes that were made to prepare the file for use on Cavatica:

- Removed many parameters that are not expected to be used on Cavatica
  - This includes most parameters related to rerunning with cached results, as we do not currently expect rerunning to be a common use case on Cavatica, and we are not certain that we will be able to support reusing `checkpoint` files.
  - Docker image parameters are also removed, so that the default values will always be used.
  - In the `nextflow_schema.json` file, we use the `"hidden": true` property to hide these parameters, but this option does not appear to be supported by Cavatica.
- Modified required parameters to remove the `null` type. This seems to be the way to designate a parameter as required in Cavatica.
  - For example, the `run_metafile` and `sample_metafile` parameters are required, so they are changed to:

  ```yaml
  type: File
  ```

  instead of the original:

  ```yaml
  type:
    - "null"
    - File
  ```

- Updated the list of profiles to include only the profiles to be used on Cavatica: `cavatica` and `stub` for now.
- Added a customized `outputs` block to fill in the correct path for the output directory from the `--outdir` parameter:

```yaml
outputs:
  - id: output_dir
    label: Workflow output directory
    outputBinding:
      glob: $(inputs.outdir)
    type: Directory
```

- Add default values to options using the `default:` property (which is distinct form the `sbg:toolDefaultValue:` property that is only used for documentation):
  - set default value for the `--outdir` to `scpca_out` to avoid errors if this is left blank
  - Set the default `-profile` to `cavatica`



#### Pushing to Cavatica

To push the workflow to Cavatica, use the following command:

```bash
op run -- ./cavatica-build.sh push
```


### API token

To use the Cavatica API, you will need to create an account at https://cavatica.sbgenomics.com and obtain an authentication token.
You can do this by going to https://cavatica.sbgenomics.com/developer/token and copying the authentication token from there.
We recommend storing this token in 1Password or a similar password manager, along with the API URL, which is `https://cavatica.sbgenomics.com/api/v2/`.

Then you can set the environment variables `SB_AUTH_TOKEN` and `SB_API_ENDPOINT` to allow access to these values.
Again, we recommend doing this via a password manager that will only expose these values as needed, e.g. with the `op run` command from the [1Password CLI](https://developer.1password.com/docs/cli/).
