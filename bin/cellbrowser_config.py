#!/usr/bin/env python3

import argparse
import csv
import pathlib
import textwrap


def write_project_config(
    project: str, project_dir: pathlib.Path, project_metadata: dict[str, str]
) -> None:
    project_dir.mkdir(parents=True, exist_ok=True)
    short_label = project
    title = project_metadata.get("project_title", project)
    abstract = project_metadata.get("abstract", "No abstract available.")
    project_cb_file = project_dir / "cellbrowser.config"
    project_cb_file.write_text(f"shortLabel='{short_label}'\n")

    project_desc_file = project_dir / "desc.config"
    project_desc_contents = f"{title=}\n{abstract=}\n"
    if project.startswith("SCPCP"):
        project_desc_contents += (
            f"other_url='https://scpca.alexslemonade.org/projects/{project}'\n"
        )
    project_desc_file.write_text(project_desc_contents)


def main():
    parser = argparse.ArgumentParser(
        description="Create UCSC Cell Browser config files"
    )
    parser.add_argument(
        "--output-dir",
        help="Path to the output directory.",
        type=pathlib.Path,
        default=pathlib.Path("."),
    )
    parser.add_argument("--projects", help="Comma separated list of project ids.")
    parser.add_argument("--project-metadata", help="Path to the project metadata file.")
    args = parser.parse_args()

    # Create the site metadata files
    site_cb_file = args.output_dir / "cellbrowser.config"
    site_cb_file.write_text("shortLabel='ScPCA Portal'\n")

    site_desc_file = args.output_dir / "desc.config"
    site_desc_contents = """\
        title='Single-cell Pediatric Cancer Atlas Portal'
        other_url='https://https://scpca.alexslemonade.org'
    """
    site_desc_file.write_text(textwrap.dedent(site_desc_contents))

    # Read project metadata and index by project ID
    project_info = {}
    with open(args.project_metadata, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            project_info[row["scpca_project_id"]] = row

    # create directories and config files for each project
    projects = args.projects.split(",") if args.projects else []
    for project in projects:
        if project in project_info:
            project_metadata = project_info[project]
        else:
            Warning(f"Warning: Project {project} not found in metadata.")
            project_metadata = {}
        project_dir = args.output_dir / project
        write_project_config(project, project_dir, project_metadata)


if __name__ == "__main__":
    main()
