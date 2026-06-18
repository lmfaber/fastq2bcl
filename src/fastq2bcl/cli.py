"""
CLI for fastq2bcl app

``[options.entry_points]`` section in ``setup.cfg``::

    console_scripts =
        fastq2bcl = fastq2bcl.cli:run

Then run ``pip install .`` (or ``pip install -e .`` for editable mode)
which will install the command ``fastq2bcl`` inside your current environment.

References:
    - https://setuptools.pypa.io/en/latest/userguide/entry_point.html
    - https://pip.pypa.io/en/stable/reference/pip_install
"""

import signal
import argparse
import csv
import logging
import sys
import os
import re
import textwrap

from pathlib import Path
from rich import print, pretty
from fastq2bcl import __version__
from fastq2bcl.parser import parse_seqdesc_fields
from fastq2bcl.reader import (
    read_first_record,
    get_mask_from_files,
    iter_fastq_records,
    count_fastq_records,
)
from fastq2bcl.writer import (
    write_run_info_xml,
    write_filter,
    write_control,
    write_sample_sheet,
    write_lane_bcls_and_locs,
)

__author__ = "Davide Rambaldi"
__copyright__ = "Davide Rambaldi"
__license__ = "MIT"

_logger = logging.getLogger(__name__)

def _is_lane_token(value):
    if value == None:
        return False
    parts = re.split(r"[;,]", value)
    numeric_parts = [part.strip() for part in parts if part.strip()]
    return bool(numeric_parts) and all(part.isdigit() for part in numeric_parts)


def _parse_lanes(values):
    lanes = []
    for value in values:
        for part in re.split(r"[;,]", value):
            part = part.strip()
            if part:
                lanes.append(int(part))
    if not lanes:
        raise ValueError("Input samplesheet row does not define any lanes")
    if len(lanes) != len(set(lanes)):
        raise ValueError("Input samplesheet row assigns the same lane more than once")
    invalid_lanes = [lane for lane in lanes if lane < 1 or lane > 8]
    if invalid_lanes:
        raise ValueError(
            "Input samplesheet lanes must be between 1 and 8; invalid lanes: "
            + ",".join(str(lane) for lane in sorted(set(invalid_lanes)))
        )
    return lanes


def _optional_path(value, base_dir):
    if value == None or value == "":
        return None
    path = Path(value)
    if not path.is_absolute():
        path = base_dir / path
    return path


def _sample_id_from_path(path):
    name = Path(path).name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return re.sub(r"([._-])R1([._-].*)?$", "", name)


def parse_input_samplesheet(input_samplesheet):
    """
    Parse lane,r1,r2,i1,i2[,project][,sample] input samplesheet rows.

    Lane lists can be semicolon-separated in one cell or comma-separated across
    leading numeric CSV cells, e.g. ``1;2;3;4,5,6,7,file_R1.fastq.gz,...``.
    Additional columns are ignored.
    """
    samplesheet = Path(input_samplesheet)
    base_dir = samplesheet.parent
    samples = []
    with open(samplesheet, newline="") as f_in:
        reader = csv.reader(f_in)
        header = next(reader, None)
        normalized_header = [h.strip() for h in header] if header != None else []
        if normalized_header[:5] != [
            "lane",
            "r1",
            "r2",
            "i1",
            "i2",
        ]:
            raise ValueError("Input samplesheet header must start with: lane,r1,r2,i1,i2")
        value_header = normalized_header[1:]

        for row_number, row in enumerate(reader, start=2):
            if not row or all(cell.strip() == "" for cell in row):
                continue
            first_fastq_idx = None
            for idx, value in enumerate(row):
                if not _is_lane_token(value.strip()):
                    first_fastq_idx = idx
                    break
            if first_fastq_idx == None:
                raise ValueError(f"Input samplesheet row {row_number} has no R1 path")

            lane_values = row[:first_fastq_idx]
            values = [cell.strip() for cell in row[first_fastq_idx:]]
            row_data = {name: "" for name in value_header}
            row_data.update(dict(zip(value_header, values)))
            lanes = _parse_lanes(lane_values)

            r1 = _optional_path(row_data["r1"], base_dir)
            if r1 == None:
                raise ValueError(f"Input samplesheet row {row_number} has no R1 path")
            sample_id = row_data.get("sample") or _sample_id_from_path(r1)
            samples.append(
                {
                    "lanes": lanes,
                    "r1": r1,
                    "r2": _optional_path(row_data["r2"], base_dir),
                    "i1": _optional_path(row_data["i1"], base_dir),
                    "i2": _optional_path(row_data["i2"], base_dir),
                    "sample_id": sample_id,
                    "sample_name": sample_id,
                    "project": row_data.get("project", ""),
                }
            )

    if not samples:
        raise ValueError("Input samplesheet does not contain any samples")
    return samples


def _samples_by_lane(samples):
    lane_samples = {}
    for sample in samples:
        for lane in sample["lanes"]:
            lane_samples.setdefault(lane, []).append(sample)
    return lane_samples


def _sample_record_count(sample, exclude_umi, exclude_index):
    cache_key = "_record_count"
    if cache_key not in sample:
        sample[cache_key] = count_fastq_records(
            sample["r1"],
            sample["r2"],
            sample["i1"],
            sample["i2"],
            exclude_umi,
            exclude_index,
        )
    return sample[cache_key]


def _sample_lane_bounds(sample, lane, exclude_umi, exclude_index):
    lanes = sample["lanes"]
    if len(lanes) == 1:
        return 0, _sample_record_count(sample, exclude_umi, exclude_index)

    lane_index = lanes.index(lane)
    record_count = _sample_record_count(sample, exclude_umi, exclude_index)
    chunk_size, extra_records = divmod(record_count, len(lanes))
    start = (lane_index * chunk_size) + min(lane_index, extra_records)
    end = start + chunk_size + (1 if lane_index < extra_records else 0)
    return start, end


def _iter_sample_lane_records(sample, lane, exclude_umi, exclude_index):
    start, end = _sample_lane_bounds(sample, lane, exclude_umi, exclude_index)
    if start == end:
        return

    for record_index, record in enumerate(
        iter_fastq_records(
            sample["r1"],
            sample["r2"],
            sample["i1"],
            sample["i2"],
            exclude_umi,
            exclude_index,
        )
    ):
        if record_index >= end:
            break
        if record_index >= start:
            yield record


def _iter_lane_records(lane, lane_samples, exclude_umi, exclude_index):
    for sample in lane_samples:
        yield from _iter_sample_lane_records(sample, lane, exclude_umi, exclude_index)


def _count_lane_records(lane, lane_samples, exclude_umi, exclude_index):
    total = 0
    for sample in lane_samples:
        start, end = _sample_lane_bounds(sample, lane, exclude_umi, exclude_index)
        total += end - start
    return total


def _sample_index_key(sample):
    index_parts = []
    try:
        seqdesc_fields = parse_seqdesc_fields(read_first_record(sample["r1"]).description)
    except StopIteration:
        return None
    if seqdesc_fields["index"] != "1":
        index_parts.append(seqdesc_fields["index"])
    if sample["i1"] != None:
        index_parts.append(str(read_first_record(sample["i1"]).seq))
    if sample["i2"] != None:
        index_parts.append(str(read_first_record(sample["i2"]).seq))
    if not index_parts:
        return None
    return tuple(index_parts)


def _validate_shared_lane_indexes(samples_by_lane, exclude_index):
    for lane, lane_samples in samples_by_lane.items():
        if len(lane_samples) < 2:
            continue
        if exclude_index:
            raise ValueError(
                f"Lane {lane} has multiple samples; --exclude-index would remove "
                "the indexes required to demultiplex them"
            )

        seen_indexes = {}
        for sample in lane_samples:
            index_key = _sample_index_key(sample)
            if index_key == None:
                if _sample_record_count(sample, False, exclude_index) == 0:
                    continue
                raise ValueError(
                    f"Lane {lane} has multiple samples, but sample " f"{sample['sample_id']} does not define an index"
                )
            if index_key in seen_indexes:
                raise ValueError(
                    f"Lane {lane} has multiple samples with duplicate index "
                    f"{'/'.join(index_key)}: {seen_indexes[index_key]} and "
                    f"{sample['sample_id']}"
                )
            seen_indexes[index_key] = sample["sample_id"]


def _samples_from_fastq_args(r1, r2, i1, i2):
    return [
        {
            "lanes": [1],
            "r1": Path(r1),
            "r2": Path(r2) if r2 != None else None,
            "i1": Path(i1) if i1 != None else None,
            "i2": Path(i2) if i2 != None else None,
            "sample_id": _sample_id_from_path(r1),
            "sample_name": _sample_id_from_path(r1),
            "project": "",
        }
    ]


def fastq2bcl(
    outdir,
    r1,
    r2=None,
    i1=None,
    i2=None,
    mask_string=None,
    exclude_umi=False,
    exclude_index=False,
    threads=1,
    input_samplesheet=None,
    sample_sheet_format="bcl2fastq",
):
    """fastq2bcl function call

    :param outdir: output directory to create run flowcell fake dir
    :param r1: R1 fastq.gz
    :param r2: R2 fastq.gz
    :param i1: I1 fastq.gz
    :param i2: I2 fastq.gz

    Content of returned tuple:

    rundir: final absolute path od created rundir
    run_id: generated mock run_id
    seq_fields: fields parsed and validated from first R1 record
    mask_string: mask used to generate RunInfo.xml

    :rtype: tuple
    """

    # First validate outdir
    outdir = Path(outdir).absolute()
    assert outdir.is_dir()
    assert os.access(outdir, os.W_OK)

    _logger.info(f"Output directory: {outdir}")

    if input_samplesheet:
        samples = parse_input_samplesheet(input_samplesheet)
    else:
        if r1 == None:
            raise ValueError("Either r1 or input_samplesheet must be provided")
        samples = _samples_from_fastq_args(r1, r2, i1, i2)

    first_sample = samples[0]
    first_r1 = first_sample["r1"]

    # Validate R1 and extract first read
    assert first_r1.is_file()
    first_record = read_first_record(first_r1)
    seqdesc_fields = parse_seqdesc_fields(first_record.description)
    _logger.info(f"first record seq length: {len(first_record.seq)}")
    _logger.info(f"first record sequence: {str(first_record.seq)}")
    _logger.info(f"first record seqdesc fields: {seqdesc_fields}")

    print(f"[green]First sequence length[/green]: {len(first_record.seq)}")
    print(f"[green]First sequence[/green]:")
    print(textwrap.fill(str(first_record.seq), 50))

    # SEE LSO docs/flow.drawio.png
    # CHECK INDEX
    if seqdesc_fields["index"] != "1":
        if exclude_index:
            print(
                f"[red]Founded INDEX sequence in first record[/red] {seqdesc_fields['index']}",
                f"[red]index sequences will NOT be included in the cycles[/red]",
            )
        else:
            print(
                f"[green]Founded INDEX sequence in first record[/green]: {seqdesc_fields['index']}",
                f"[green]index sequences will be included in the cycles[/green]",
            )

    # CHECK UMI
    if seqdesc_fields["UMI"] != None:
        if exclude_umi:
            print(
                f"[red]Founded UMI sequence in first record[/red] {seqdesc_fields['UMI']}",
                f"[red]umi sequences will NOT be included in the cycles[/red]",
            )
        else:
            print(
                f"[green]Founded UMI sequence in first record[/green]: {seqdesc_fields['UMI']}",
                f"[green]umi sequences will be included in the cycles[/green]",
            )

    # RUNDIR
    run_id = mock_run_id(seqdesc_fields)
    rundir = Path.joinpath(outdir, run_id)

    print(f"[green]RUNDIR[/green]: {rundir}")

    if not mask_string:
        # get cycles string from files
        mask_string = get_mask_from_files(
            first_sample["r1"],
            first_sample["r2"],
            first_sample["i1"],
            first_sample["i2"],
            exclude_umi,
            exclude_index,
        )
        _logger.info(f"mask string from files: {mask_string}")

    print(f"[green]MASK[/green]: {mask_string}")

    # SET MASK FROM STRING
    mask = set_mask(mask_string)
    cycles = sum([int(m["cycles"]) for m in mask])
    samples_by_lane = _samples_by_lane(samples)
    _validate_shared_lane_indexes(samples_by_lane, exclude_index)
    lane_count = max(samples_by_lane)

    # WRITE RUN INFO
    _logger.info(f"Writing RunInfo.mxl to dir: {rundir}")
    run_info = write_run_info_xml(
        rundir,
        run_id,
        seqdesc_fields["run_number"],
        seqdesc_fields["flowcell_id"],
        seqdesc_fields["instrument"],
        mask,
        lane_count,
    )

    print(f"[green]RunInfo.xml:[/green]:\n", run_info)
    write_sample_sheet(rundir, samples, mask, sample_sheet_format, exclude_index)

    for lane, lane_samples in sorted(samples_by_lane.items()):
        cluster_count = _count_lane_records(lane, lane_samples, exclude_umi, exclude_index)

        # WRITE FILTER
        print(f"[bold magenta]Writing filter file for lane {lane}[/bold magenta]")
        _logger.info(f"Writing filter file to dir: {rundir} with cluster count: {cluster_count}")
        write_filter(rundir, cluster_count, lane)

        # WRITE CONTROL
        print(f"[bold magenta]Writing control file for lane {lane}[/bold magenta]")
        _logger.info(f"Writing control file to dir: {rundir} with cluster count: {cluster_count}")
        write_control(rundir, cluster_count, lane)

        # WRITE LOCATIONS, BCL AND STATS
        print(f"[bold magenta]Writing location and cycle files for lane {lane}[/bold magenta]")
        _logger.info(f"Writing {cluster_count} records to dir: {rundir}")
        write_lane_bcls_and_locs(
            rundir,
            cluster_count,
            cycles,
            _iter_lane_records(lane, lane_samples, exclude_umi, exclude_index),
            lane,
        )

    return run_id, rundir, seqdesc_fields, mask_string


def mock_run_id(fields):
    """
    Mock the run directory id and Path
    """
    run_id = "100101_" + fields["instrument"] + "_" + fields["run_number"].zfill(4) + "_" + fields["flowcell_id"]
    return run_id


def set_mask(mask_string):
    if mask_string:
        mask = []
        regexp_mask = r"([0-9]+[NY])([0-9]+[NY])?([0-9]+[NY])?([0-9]+[NY])?"
        m = re.match(regexp_mask, mask_string)
        if not m:
            raise ValueError(f"Incorrect mask parse: {mask_string}")
        reads = [g for g in m.groups() if g != None]
        for g_idx in range(len(reads)):
            read = re.match(r"([0-9]+)([YN])", reads[g_idx])
            mask.append(
                {
                    "cycles": read.groups()[0],
                    "index": read.groups()[1],
                    "id": str(g_idx + 1),
                }
            )
        return mask
    else:
        raise ValueError(f"Incorrect mask string: {mask_string}")


# ---- CLI ----
# The functions defined in this section are wrappers around the main Python
# API allowing them to be called directly from the terminal as a CLI
# executable/script.


def parse_args(args):
    """Parse command line parameters

    Args:
      args (List[str]): command line parameters as list of strings
          (for example  ``["--help"]``).

    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """
    parser = argparse.ArgumentParser(
        description="Convert fastq.gz reads and metadata in a bcl2fastq-able run directory"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"fastq2bcl {__version__}",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="loglevel",
        help="set loglevel to INFO",
        action="store_const",
        const=logging.INFO,
    )
    parser.add_argument(
        "-vv",
        "--very-verbose",
        dest="loglevel",
        help="set loglevel to DEBUG",
        action="store_const",
        const=logging.DEBUG,
    )
    parser.add_argument("-m", "--mask", dest="mask", help="define mask in format 110N10Y10Y110N")

    parser.add_argument(
        "-r1",
        "--read-1",
        dest="r1",
        help="fastq.gz with R1 reads",
        metavar="R1",
    )

    parser.add_argument(
        "--input_samplesheet",
        dest="input_samplesheet",
        help="CSV with lane,r1,r2,i1,i2 rows assigning FASTQ files to lanes",
    )

    parser.add_argument(
        "--sample-sheet-format",
        dest="sample_sheet_format",
        choices=["bcl2fastq", "bcl-convert"],
        default="bcl2fastq",
        help="Generated SampleSheet.csv format. default: bcl2fastq",
    )

    parser.add_argument(
        "-r2",
        "--read-2",
        dest="r2",
        help="fastq.gz with R2 reads (optional)",
        metavar="R2",
    )

    parser.add_argument(
        "-i1",
        "--index-1",
        dest="i1",
        help="fastq.gz with I1 reads (optional)",
        metavar="I1",
    )

    parser.add_argument(
        "-i2",
        "--index-2",
        dest="i2",
        help="fastq.gz with I2 reads (optional)",
        metavar="I2",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        help="Set the output directory for mocked run. default: cwd",
        default=os.getcwd(),
    )

    parser.add_argument(
        "--exclude-umi",
        dest="exclude_umi",
        help="Do not write UMI from the R1 and R2 fastq reads to the cycles",
        action="store_true",
    )

    parser.add_argument(
        "--exclude-index",
        dest="exclude_index",
        help="Do not write Index from the R1 and R2 fastq reads to the cycles",
        action="store_true",
    )

    parser.add_argument(
        "-T",
        "--threads",
        help="Number of threads to use to write bcls. Default 1",
        type=int,
        default=1,
        dest="threads",
    )

    parsed_args = parser.parse_args(args)
    if parsed_args.input_samplesheet == None and parsed_args.r1 == None:
        parser.error("Either --input_samplesheet or -r1/--read-1 is required")
    return parsed_args


def setup_logging(loglevel):
    """Setup basic logging

    Args:
      loglevel (int): minimum loglevel for emitting messages
    """
    logformat = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s"
    logging.basicConfig(level=loglevel, stream=sys.stdout, format=logformat, datefmt="%Y-%m-%d %H:%M:%S")


def main(args):
    """Wrapper allowing :func:`fastq2bcl` to be called with string arguments in a CLI fashion

    Instead of returning the value from :func:`fastq2bcl`, it prints the result to the
    ``stdout`` in a nicely formatted message.

    Args:
      args (List[str]): command line parameters as list of strings
          (for example  ``["--verbose", "42"]``).
    """
    pretty.install()

    args = parse_args(args)
    setup_logging(args.loglevel)
    _logger.info("Starting application...")
    _logger.info(f"User defined mask: {args.mask}")
    _logger.info(f"Input files: R1={args.r1} R2={args.r2} I1={args.i1} I2={args.i2}")

    print("[bold green]fastq2bcl[/bold green]")
    print("Args:", args)

    # call fastq2bcl
    run_id, rundir, seqdesc_fields, mask_string = fastq2bcl(
        args.outdir,
        args.r1,
        args.r2,
        args.i1,
        args.i2,
        args.mask,
        args.exclude_umi,
        args.exclude_index,
        args.threads,
        args.input_samplesheet,
        args.sample_sheet_format,
    )

    _logger.info("Script ends here")


def run():
    """Calls :func:`main` passing the CLI arguments extracted from :obj:`sys.argv`

    This function can be used as entry point to create console scripts with setuptools.
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
