# biopq/vcf_converter.py
from __future__ import annotations

import itertools
from pathlib import Path
from typing import Iterable, Iterator

import polars as pl
import pyarrow as pa
import pyarrow.dataset as ds
from cyvcf2 import VCF
from loguru import logger
from pydantic import BaseModel, ValidationError, create_model
from tqdm import tqdm


def build_schema_from_header(vcf_header: VCF.Header) -> type[BaseModel]:
    """
    Dynamically builds a Pydantic model from a VCF header's INFO fields.

    Args:
        vcf_header (VCF.Header): The VCF header object from cyvcf2.

    Returns:
        type[BaseModel]: A dynamically created Pydantic model class for INFO field validation.
    """
    type_map = {
        "Integer": (int | None, None),
        "Float": (float | None, None),
        "String": (str | None, None),
        "Flag": (bool, False),  # Flags are present (True) or absent (False)
    }

    fields = {}
    for header_line in vcf_header.get_lines("INFO"):
        info_id = header_line.get("ID")
        info_type = header_line.get("Type")
        python_type, default_value = type_map.get(info_type, (str | None, None))
        fields[info_id] = (python_type, default_value)

    VcfInfoModel: type[BaseModel] = create_model("VcfInfoModel", **fields)
    return VcfInfoModel


def _record_generator(
    vcf_iterator: Iterable[VCF.Variant], VcfInfoModel: type[BaseModel]
) -> Iterable[tuple]:
    """
    Generator that validates VCF records and yields flat tuples.

    Args:
        vcf_iterator (Iterator[VCF.Variant]): An iterator over VCF records (can be wrapped with tqdm).
        VcfInfoModel (type[BaseModel]): The dynamic Pydantic model for validating INFO fields.

    Yields:
        Iterable[tuple]: A stream of flat tuples, each representing one allele of a variant.
    """
    for rec in vcf_iterator:
        pos1 = rec.POS
        start = pos1 - 1
        end_val = rec.INFO.get("END")
        try:
            end = int(end_val)
        except (ValueError, TypeError, SystemError):
            end = start + len(rec.REF)

        try:
            validated_info = VcfInfoModel(**rec.INFO)
            info_values = validated_info.model_dump().values()
        except ValidationError as e:
            logger.warning(
                f"Schema validation failed at {rec.CHROM}:{rec.POS}, skipping record. Error: {e}"
            )
            continue

        for alt in rec.ALT or []:
            yield (rec.CHROM, start, end, pos1, rec.REF, alt, *info_values)


def _record_batch_generator(
    record_iterator: Iterable[tuple], schema: list[str], chunk_size: int
) -> Iterator[pa.RecordBatch]:
    """
    Generator that converts a stream of tuples into Arrow RecordBatches.

    Args:
        record_iterator (Iterable[tuple]): A stream of flat tuples from _record_generator.
        schema (list[str]): The column names for the DataFrame.
        chunk_size (int): The number of records to include in each batch.

    Yields:
        Iterator[pa.RecordBatch]: A stream of Arrow RecordBatches.
    """
    while True:
        chunk = list(itertools.islice(record_iterator, chunk_size))
        if not chunk:
            break

        df_chunk = pl.DataFrame(chunk, schema=schema, orient="row")
        table_chunk = df_chunk.to_arrow()
        yield from table_chunk.to_batches()


def convert_vcf_to_parquet(
    vcf_file: Path,
    output_path: Path,
    chunk_size: int = 500_000,
    regions: str | None = None,
) -> None:
    """
    Converts a VCF/BCF file to a partitioned Parquet dataset.

    This function orchestrates a streaming, chunked, and schema-validated pipeline
    to process large VCF files with low and constant memory usage.

    Args:
        vcf_file (Path): Path to the input VCF/BCF file.
        output_path (Path): Path for the output Parquet dataset directory.
        chunk_size (int): Number of VCF records to process per memory chunk.
        regions (str, optional): A specific genomic region to process (e.g., "chr1:1-10000").
    """
    logger.info(f"Starting conversion for {vcf_file.name}")
    vcf = VCF(str(vcf_file), lazy=True)

    # Build the dynamic validation schema from header
    logger.info("Building dynamic schema from VCF header...")
    VcfInfoModel = build_schema_from_header(vcf.header_iter)
    info_keys = list(VcfInfoModel.model_fields.keys())
    final_schema = ["chrom", "start", "end", "pos1", "ref", "alt"] + info_keys

    # Set up progress bar (tqdm)
    total_variants = None
    if vcf.index_fn:
        try:
            total_variants = len(vcf)
            logger.info(f"Found index. Total variants to process: {total_variants:,}")
        except Exception:
            logger.warning("Index found but failed to read record count.")
    else:
        logger.warning(
            "No index found; progress bar will not show percentage complete."
        )

    vcf_iterator = vcf(regions) if regions else vcf
    pbar = tqdm(
        vcf_iterator, total=total_variants, desc="Processing variants", unit=" variants"
    )

    # Chain processing functions together
    record_iterator = _record_generator(pbar, VcfInfoModel)
    record_batches = _record_batch_generator(record_iterator, final_schema, chunk_size)

    # Prime generator to get schema from first batch
    try:
        first_batch = next(record_batches)
    except StopIteration:
        logger.warning(f"No valid records found in {vcf_file.name}")
        return

    # Write streaming data to a partitioned Parquet dataset
    final_iterator = itertools.chain([first_batch], record_batches)
    logger.info(f"Writing partitioned Parquet dataset to {output_path}...")
    ds.write_dataset(
        data=final_iterator,
        base_dir=output_path,
        schema=first_batch.schema,
        format="parquet",
        partitioning=["chrom"],
        existing_data_behavior="overwrite_or_ignore",
    )
    logger.info("✅ Conversion complete!")
