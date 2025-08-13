#!/usr/bin/env python3

import argparse
import csv
import itertools
import pathlib
import textwrap


def write_project_configs(
    project: str, project_dir: pathlib.Path, project_metadata: dict[str, str]
) -> None:
    project_dir.mkdir(parents=True, exist_ok=True)
    # cellbrowser.conf
    short_label = project
    title = project_metadata.get("project_title", project)
    abstract = project_metadata.get("abstract", "No abstract available.")
    project_cb_file = project_dir / "cellbrowser.conf"
    project_cb_file.write_text(f"shortLabel='{short_label}'\n")
    # desc.conf
    project_desc_file = project_dir / "desc.conf"
    project_desc_contents = f"{title=}\n{abstract=}\nhideDownload=True\n"
    if project.startswith("SCPCP"):
        project_desc_contents += f"other_url='https://scpca.alexslemonade.org/projects/{project} ScPCA project page for {project}'\n"
    project_desc_file.write_text(project_desc_contents)


def write_library_configs(
    library_id: str, sample_id: str, library_dir: pathlib.Path, label_field: str
) -> None:
    library_dir.mkdir(parents=True, exist_ok=True)
    # cellbrowser.conf
    library_cb_file = library_dir / "cellbrowser.conf"
    library_cb_contents = f"""\
        name='{library_id}'
        shortLabel='{sample_id} {library_id}'
        exprMatrix='matrix.mtx.gz'
        meta='meta.tsv'
        geneIdType='auto'
        defColorField='{label_field}'
        labelField='{label_field}'
        enumFields=['{label_field}']
        coords=[
            {"file": "umap_coords.tsv", "shortLabel": "UMAP"}
        ]
        markers = [{"file": "markers.tsv", "shortLabel":"Cluster Markers"}]
        quickGenesFile="quickGenes.tsv"
    """
    library_cb_file.write_text(textwrap.dedent(library_cb_contents))
    # desc.conf
    library_desc_file = library_dir / "desc.conf"
    library_desc_contents = f"""\
        title='{sample_id} {library_id}'
        hideDownload=True
    """
    library_desc_file.write_text(textwrap.dedent(library_desc_contents))


def main():
    parser = argparse.ArgumentParser(
        description="Create UCSC Cell Browser config files"
    )
    parser.add_argument(
        "--outdir", "-o", help="Path to the output directory.", type=pathlib.Path
    )
    parser.add_argument(
        "--conf_type",
        type=str,
        choices=["project", "library"],
        help="Type of configuration files to write ('project' or 'library').",
    )
    parser.add_argument(
        "--ids", required=True, help="Comma separated list of project or library ids."
    )
    parser.add_argument("--sample_ids", help="Comma separated list of sample ids.")
    parser.add_argument("--project-metadata", help="Path to the project metadata file.")
    parser.add_argument(
        "--label-field",
        help="Name of the display label (grouping) field.",
        default="cluster",
    )
    args = parser.parse_args()

    # Arg consistency checks
    if args.conf_type == "project" and not args.project_metadata:
        parser.error(
            "Project metadata file is required to generate project config files."
        )

    args.outdir.mkdir(parents=True, exist_ok=True)
    ids = args.ids.split(",") if args.ids else []
    sample_ids = args.sample_ids.split(",") if args.sample_ids else []

    if args.conf_type == "project":
        # Read project metadata and index by project ID
        project_info = {}
        with open(args.project_metadata, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                project_info[row["scpca_project_id"]] = row

        # create directories and config files for each project
        for project in ids:
            if project in project_info:
                project_metadata = project_info[project]
            else:
                Warning(f"Warning: Project {project} not found in metadata.")
                project_metadata = {}
            project_dir = args.outdir / project
            write_project_configs(project, project_dir, project_metadata)

    elif args.conf_type == "library":
        # create directories and config files for each library
        for library, sample in itertools.zip_longest(ids, sample_ids):
            library_dir = args.outdir / library
            sample = sample if sample else ""
            write_library_configs(library, sample, library_dir, args.label_field)


if __name__ == "__main__":
    main()
