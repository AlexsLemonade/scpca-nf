# Cavatica port documentation

## Prerequisites

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
cavatica-build.sh
```

This will use the existing file as a starting point in order to preserve any customizations.
If you want to start fresh, you can delete the existing `sb_nextflow_schema.yaml` file and run the command again.

#### Pushing to Cavatica

To push the workflow to Cavatica, use the following command:

```bash
op run -- cavatica-build.sh push
```


### API token

To use the Cavatica API, you will need to create an account at https://cavatica.sbgenomics.com and obtain an authentication token.
You can do this by going to https://cavatica.sbgenomics.com/developer/token and copying the authentication token from there.
We recommend storing this token in 1Password or a similar password manager, along with the API URL, which is `https://cavatica.sbgenomics.com/api/v2/`.

Then you can set the environment variables `SB_AUTH_TOKEN` and `SB_API_ENDPOINT` to allow access to these values.
Again, we recommend doing this via a password manager that will only expose these values as needed, e.g. with the `op run` command from the [1Password CLI](https://developer.1password.com/docs/cli/).
