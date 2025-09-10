# biopq/cli.py
import sys
from pathlib import Path
from typing import Annotated

import typer
from loguru import logger

from biopq.vcf_converter import convert_vcf_to_parquet

app = typer.Typer(
    name="biopq",
    help="A modern tool to convert genomic files to Parquet.",
    add_completion=False,
    pretty_exceptions_show_locals=False,
)


@app.command()
def vcf(
    vcf_file: Annotated[
        Path, typer.Option("--vcf", "-i", exists=True, help="Input VCF/BCF file.")
    ],
    output_path: Annotated[
        Path, typer.Option("--out", "-o", help="Output Parquet directory.")
    ],
    verbose: Annotated[
        bool, typer.Option("--verbose", "-v", help="Enable detailed (DEBUG) logging.")
    ] = False,
):
    """
    Convert a VCF/BCF file to a partitioned Parquet dataset.
    """
    log_level = "DEBUG" if verbose else "INFO"
    logger.remove()
    logger.add(sys.stderr, level=log_level)

    try:
        logger.debug(f"Full input path: {vcf_file.resolve()}")
        logger.debug(f"Full output path: {output_path.resolve()}")

        convert_vcf_to_parquet(vcf_file, output_path)

        print(f"\n✅ Success! Parquet dataset created at: {output_path}")

    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")
        print("\n❌ Error: Conversion failed. See the log for details.")
        raise typer.Exit(code=1)


@app.command()
def gtf():
    """
    (Coming soon) Convert a GTF/GFF file to Parquet.
    """
    typer.echo("GTF converter not yet implemented.")


if __name__ == "__main__":
    app()
