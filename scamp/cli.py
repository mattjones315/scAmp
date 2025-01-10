"""
Command-line tools for scAmp.
"""

from __future__ import annotations

import os
import pathlib
import pickle
from typing import Annotated, Union

import typer

from scamp import io
from scamp.mixins import CLIError
from scamp import predict

scamp_app = typer.Typer(help="Tools for single-cell analysis of ecDNA.")

CopyNumberFile = Annotated[
    str, typer.Option(help="File path to tab-delimited file of copy numbers")
]
FragDirArg = Annotated[
    str,
    typer.Argument(help="Path to directory containing ATAC fragment files."),
]
ModelFileArg = Annotated[str, typer.Argument(help="Path to model file.")]
OutputDirArg = Annotated[str, typer.Argument(help="Directory of output files.")]
WhitelistFileArg = Annotated[
    str, typer.Argument(help="File path to cellBC whitelist.")
]


@scamp_app.command(name="atac-cnv", help="Quantify single-cell copy-numbers.")
def quantify_copy_numbers(
    fragment_directory: FragDirArg,
    output_directory: OutputDirArg,
    whitelist_file: WhitelistFileArg,
    window_size: Annotated[
        int, typer.Option(help="Base pair width for genomic windows")
    ] = 3000000,
    step_size: Annotated[
        int,
        typer.Option(help="Step size from previous genomic window."),
    ] = 1000000,
    n_neighbors: Annotated[
        int,
        typer.Option(
            help="Number of genomic windows to compare against for normalization"
        ),
    ] = 200,
    reference_genome_name: Annotated[
        str,
        typer.Option(help="Reference genome name, to pair with a blacklist."),
    ] = "hg38",
):

    binned_copy_number_script = (
        f"{os.path.dirname(__file__)}/utilities/scATAC_CNV.R"
    )
    gene_aggregation_script = (
        f"{os.path.dirname(__file__)}/utilities/aggregate_gene_copy_numbers.R"
    )

    # compute copy-numbers in genomic windows
    os.system(
        f"Rscript {binned_copy_number_script} "
        f"{fragment_directory} {window_size} "
        f"{step_size} {n_neighbors} "
        f"{whitelist_file} {output_directory} "
        f"{reference_genome_name} {os.path.dirname(__file__)}"
    )

    # bin by gene
    os.system(
        f"Rscript {gene_aggregation_script} "
        f"{output_directory} {output_directory} {reference_genome_name}"
    )


@scamp_app.command(name="predict", help="Predict ecDNA status.")
def predict_ecdna(
    output_dir: OutputDirArg,
    model_file: ModelFileArg,
    copy_numbers_file: CopyNumberFile,
    anndata: Annotated[
        str, typer.Option("File path to Anndata with copy-number data")
    ] = None,
    mode: Annotated[
        str, typer.Option("Mode: (currently only offering `copynumber`)")
    ] = "copynumber",
    min_copy_number: Annotated[
        float, typer.Option(help="Minimum copy-number to consider.")
    ] = 2.0,
    max_percentile: Annotated[
        float, typer.Option(help="Maximum percentile to cap copy-numbers.")
    ] = 99.0,
) -> None:

    if (copy_numbers_file is None) and (anndata_file is None):
        raise CLIError("Specify one of copy numbers file anndata file.")

    if mode == "copynumber":
        preds = predict.predict_ecdna_from_copy_number(
            copy_numbers_file,
            model_file,
            min_copy_number,
            max_percentile,
        )
