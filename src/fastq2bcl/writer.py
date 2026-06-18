import logging
import struct
import csv
from contextlib import ExitStack
from pathlib import Path

from fastq2bcl.parser import parse_seqdesc_fields
from fastq2bcl.reader import read_first_record

_logger = logging.getLogger(__name__)

BASE_TO_BITS = {ord("A"): 0, ord("C"): 1, ord("G"): 2, ord("T"): 3}
BASE_TO_BITS.update({ord("a"): 0, ord("c"): 1, ord("g"): 2, ord("t"): 3})
WRITE_CHUNK_SIZE = 1024 * 1024


def format_lane(lane):
    return f"L{int(lane):03d}"


def format_tile(lane):
    return f"s_{int(lane)}_1101"


def write_run_info_xml(rundir, run_id, run_number, flowcell_id, instrument, mask, lane_count=1):
    """
    Write RunInfo.xml
    """

    runinfo = generate_run_info_xml(run_id, run_number, flowcell_id, instrument, mask, lane_count)
    _logger.info(f"RunInfo.xml:\n{runinfo}")

    # Create directory and write file
    xmlout = Path.joinpath(rundir, "RunInfo.xml")
    xmlout.parent.mkdir(exist_ok=True, parents=True)
    with open(xmlout, "wt") as f_out:
        f_out.write(runinfo)

    return runinfo


def generate_run_info_xml(run_id, run_number, flowcell_id, instrument, mask, lane_count=1):
    """
    Generate a valid Runinfo xml file.
    """

    # check mask and write mask
    xml_mask = ""
    for m in mask:
        xml_mask += f"""<Read NumCycles="{m['cycles']}" Number="{m['id']}" IsIndexedRead="{m['index']}" />"""

    xml = f"""<?xml version="1.0"?>
<RunInfo xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Version="2">
    <Run Id="{run_id}" Number="{run_number}">
        <Flowcell>{flowcell_id}</Flowcell>
        <Instrument>{instrument}</Instrument>
        <Date>100101</Date>
        <Reads>
            { xml_mask }
        </Reads>
        <FlowcellLayout LaneCount="{lane_count}" SurfaceCount="1" SwathCount="1" TileCount="1" />
    </Run>
</RunInfo>
"""
    return xml


def write_filter(rundir, cluster_count, lane=1):
    """
    Write filter
    """
    path = rundir / f"Data/Intensities/BaseCalls/{format_lane(lane)}/{format_tile(lane)}.filter"
    path.parent.mkdir(exist_ok=True, parents=True)
    with open(path, "wb") as f_out:
        f_out.write(bytes([0, 0, 0, 0]))
        f_out.write(bytes([3, 0, 0, 0]))
        f_out.write(struct.pack("<I", cluster_count))
        remaining = cluster_count
        chunk = b"\x01" * min(WRITE_CHUNK_SIZE, max(remaining, 1))
        while remaining > 0:
            write_size = min(remaining, len(chunk))
            f_out.write(chunk[:write_size])
            remaining -= write_size


def write_control(rundir, cluster_count, lane=1):
    """
    Write control file
    """
    path = rundir / f"Data/Intensities/BaseCalls/{format_lane(lane)}/{format_tile(lane)}.control"
    path.parent.mkdir(exist_ok=True, parents=True)
    with open(path, "wb") as f_out:
        f_out.write(bytes([0, 0, 0, 0]))  # "Zero value (for backwards compatibility)"
        f_out.write(bytes([2, 0, 0, 0]))  # "Format version number"
        f_out.write(struct.pack("<I", cluster_count))  # "Number of clusters"
        remaining = cluster_count * 2
        chunk = b"\x00" * min(WRITE_CHUNK_SIZE, max(remaining, 1))
        while remaining > 0:
            write_size = min(remaining, len(chunk))
            f_out.write(chunk[:write_size])
            remaining -= write_size


def write_locs(outdir, positions, lane=1):
    """
    Write locations.

    Args:
        Positions (List(tuple)): is a List of tuple with x and y values
    """
    # From mkdata.sh of bcl2fastq

    # printf '0: 010000000000803f' | xxd -r -g0 > "$locs_filename"
    # printf '0: %.8x' $clusters_count | sed -E 's/0: (..)(..)(..)(..)/0: \4\3\2\1/' | xxd -r -g0 >> "$locs_filename"

    # So with 1 cluster count should be
    # 01 00 00 00    00 00 80 3f
    # 01 00 00 00    CDCC8C3F 9A99993F

    # Source of this is bcl2fastq/src/cxx/lib/data

    # struct Record
    # {
    #     /// \brief X-coordinate.
    #     float x_;
    #     /// \brief y-coordinate.
    #     float y_;
    # }
    path = Path(outdir) / f"Data/Intensities/{format_lane(lane)}/{format_tile(lane)}.locs"
    path.parent.mkdir(exist_ok=True, parents=True)
    with open(path, "wb") as f_out:
        f_out.write(bytes([1, 0, 0, 0, 0, 0, 0x80, 0x3F]))
        f_out.write(struct.pack("<I", 0))
        positions_count = 0
        for position in positions:
            f_out.write(encode_loc_bytes(position[0], position[1]))
            positions_count += 1
        f_out.seek(8)
        f_out.write(struct.pack("<I", positions_count))
    return positions_count


def encode_loc_bytes(x_pos, y_pos):
    """
    Encode x and y positon.
    FIXME this is not the correct formul according to the bcl2fastq source code.
    """
    x_bytes = struct.pack("<f", (int(x_pos) - 1000) / 10)
    y_bytes = struct.pack("<f", (int(y_pos) - 1000) / 10)
    return x_bytes + y_bytes


def encode_cluster_byte(base, qual):
    """
    Encode cluster byte.
    Bits 0-1 are the bases, respectively [A, C, G, T]
    for [0, 1, 2, 3]:
    bits 2-7 are shifted by two bits and contain the quality score.
    All bits 0 in a byte is reserved for no-call.
    """
    encoded_base = BASE_TO_BITS.get(ord(base))
    if encoded_base == None:
        return bytes([0])  # no call
    return bytes([(qual << 2) | encoded_base])


def encode_cluster_value(base, qual):
    encoded_base = BASE_TO_BITS.get(ord(base))
    if encoded_base == None:
        return 0
    return (qual << 2) | encoded_base


def init_bcl_and_write_cluster_counts(cycledir, cluster_count, filename="s_1_1101.bcl"):
    """
    Create bcl file and write cluster count
    """
    with open(cycledir / filename, "wb") as f_out:
        f_out.write(struct.pack("<I", cluster_count))


def write_cycle(context, progress, task_id):
    """
    Write a cycle file with a thread. with progress, task_id and exit event
    context: tuple with (cycle, cluster_count, outdir, data)
    data: tuple (base, quality) for a cluster
    """
    if len(context) == 4:
        cycle, cluster_count, outdir, data = context
        lane = 1
    else:
        cycle, cluster_count, outdir, data, lane = context
    cycledir = get_cycle_dir(outdir, cycle, lane)
    filename = cycledir / f"{format_tile(lane)}.bcl"
    _logger.info(f"Writing {cluster_count} clusters for cycle: {cycle+1} to dir {cycledir}")

    init_bcl_and_write_cluster_counts(cycledir, cluster_count, f"{format_tile(lane)}.bcl")

    # write data
    sequences_written = 0
    with open(filename, "ab") as f_out:
        for base, quality in data:
            _logger.debug(f"Appending seq: {base}")
            f_out.write(encode_cluster_byte(base, quality))
            sequences_written += 1
            progress[task_id] = {
                "progress": sequences_written,
                "total": cluster_count,
            }

    # write stats
    write_stat_file(cycledir / f"{format_tile(lane)}.stats")


def write_bcl_and_stats(cycle, cluster_count, outdir, sequences, lane=1):
    """
    Single process mode to write bcls
    """
    cycledir = get_cycle_dir(outdir, cycle, lane)
    filename = cycledir / f"{format_tile(lane)}.bcl"
    init_bcl_and_write_cluster_counts(cycledir, cluster_count, f"{format_tile(lane)}.bcl")
    # write data
    with open(filename, "ab") as f_out:
        for basecalls, qualscores in sequences:
            if cycle >= len(basecalls):
                _logger.info(f"Sequence is shorter than expected, adding N")
                f_out.write(encode_cluster_byte("N", 0))
            else:
                f_out.write(encode_cluster_byte(basecalls[cycle], qualscores[cycle]))
                _logger.debug(
                    f"Appending basecall: {basecalls[cycle]} to bcl for cycle {cycle+1} lenght sequence {len(basecalls)}"
                )

    # write stats
    write_stat_file(cycledir / f"{format_tile(lane)}.stats")


def write_lane_bcls_and_locs(outdir, cluster_count, cycles, records, lane=1):
    """
    Stream lane records once and write locs plus all cycle BCL files.
    """
    tile = format_tile(lane)
    cycle_files = []
    for cycle in range(cycles):
        cycledir = get_cycle_dir(outdir, cycle, lane)
        init_bcl_and_write_cluster_counts(cycledir, cluster_count, f"{tile}.bcl")
        cycle_files.append(cycledir / f"{tile}.bcl")

    locs_path = Path(outdir) / f"Data/Intensities/{format_lane(lane)}/{tile}.locs"
    locs_path.parent.mkdir(exist_ok=True, parents=True)

    buffers = [bytearray() for _cycle in range(cycles)]
    positions_count = 0

    with ExitStack() as stack:
        locs_out = stack.enter_context(open(locs_path, "wb"))
        locs_out.write(bytes([1, 0, 0, 0, 0, 0, 0x80, 0x3F]))
        locs_out.write(struct.pack("<I", 0))

        bcl_outputs = [stack.enter_context(open(filename, "ab")) for filename in cycle_files]

        for (basecalls, qualscores), position in records:
            locs_out.write(encode_loc_bytes(position[0], position[1]))
            positions_count += 1

            read_length = len(basecalls)
            for cycle in range(cycles):
                if cycle >= read_length:
                    buffers[cycle].append(0)
                else:
                    buffers[cycle].append(encode_cluster_value(basecalls[cycle], qualscores[cycle]))

                if len(buffers[cycle]) >= WRITE_CHUNK_SIZE:
                    bcl_outputs[cycle].write(buffers[cycle])
                    buffers[cycle].clear()

        for cycle, buffer in enumerate(buffers):
            if buffer:
                bcl_outputs[cycle].write(buffer)

        locs_out.seek(8)
        locs_out.write(struct.pack("<I", positions_count))

    for filename in cycle_files:
        write_stat_file(filename.with_suffix(".stats"))

    if positions_count != cluster_count:
        raise ValueError(f"Location count {positions_count} does not match cluster count {cluster_count}")

    return positions_count


def write_stat_file(filename):
    with open(filename, "wb") as f_out:
        # can I get away with this?
        f_out.write(bytes([0] * 108))


def append_data_to_bcl(base, quality, filename):
    bcl_byte = encode_cluster_byte(base, quality)
    with open(filename, "ab") as f_out:
        f_out.write(bcl_byte)


def get_cycle_dir(outdir, cycle, lane=1):
    cycledir = outdir / f"Data/Intensities/BaseCalls/{format_lane(lane)}/C{cycle+1}.1"
    cycledir.mkdir(exist_ok=True, parents=True)
    return cycledir


def _read_labels(mask):
    labels = []
    read_count = 0
    index_count = 0
    for read in mask:
        if read["index"] == "Y":
            index_count += 1
            labels.append((f"Index{index_count}Cycles", read["cycles"]))
        else:
            read_count += 1
            labels.append((f"Read{read_count}Cycles", read["cycles"]))
    return labels


def _override_cycles(mask):
    return ";".join(f"{'I' if read['index'] == 'Y' else 'Y'}{read['cycles']}" for read in mask)


def _sample_indexes(sample, exclude_index=False):
    indexes = []
    try:
        seqdesc_fields = parse_seqdesc_fields(read_first_record(sample["r1"]).description)
    except StopIteration:
        return "", ""
    if seqdesc_fields["index"] != "1" and not exclude_index:
        indexes.extend(seqdesc_fields["index"].split("+"))
    if sample["i1"] != None:
        indexes.append(str(read_first_record(sample["i1"]).seq))
    if sample["i2"] != None:
        indexes.append(str(read_first_record(sample["i2"]).seq))
    indexes += ["", ""]
    return indexes[0], indexes[1]


def write_sample_sheet(rundir, samples, mask, sample_sheet_format="bcl2fastq", exclude_index=False):
    if sample_sheet_format == "bcl-convert":
        write_bcl_convert_sample_sheet(rundir, samples, mask, exclude_index)
    elif sample_sheet_format == "bcl2fastq":
        write_bcl2fastq_sample_sheet(rundir, samples, mask)
    else:
        raise ValueError(f"Unsupported sample sheet format: {sample_sheet_format}")


def _is_undetermined_sample(sample):
    return any("undetermined" in str(sample.get(key, "")).lower() for key in ("project", "sample_id", "sample_name"))


def write_bcl2fastq_sample_sheet(rundir, samples, mask):
    """
    Write a minimal bcl2fastq SampleSheet.csv for the generated lanes.
    """
    reads = [m["cycles"] for m in mask if m["index"] == "N"]
    path = Path(rundir) / "SampleSheet.csv"
    path.parent.mkdir(exist_ok=True, parents=True)
    with open(path, "wt") as f_out:
        f_out.write("[Header]\n\n")
        f_out.write("[Reads]\n")
        for read in reads:
            f_out.write(f"{read}\n")
        f_out.write("\n[Settings]\n\n")
        f_out.write("[Data]\n")
        f_out.write("Lane,Sample_ID,Sample_Name,Description,Sample_Project\n")
        for sample in samples:
            if _is_undetermined_sample(sample):
                continue
            for lane in sample["lanes"]:
                f_out.write(f"{lane},{sample['sample_id']},{sample['sample_name']},,\n")


def write_bcl_convert_sample_sheet(rundir, samples, mask, exclude_index=False):
    """
    Write a v2 SampleSheet.csv suitable for bcl-convert.
    """
    path = Path(rundir) / "SampleSheet.csv"
    path.parent.mkdir(exist_ok=True, parents=True)
    override_cycles = _override_cycles(mask)
    read_labels = _read_labels(mask)
    has_index1 = any(label == "Index1Cycles" for label, _cycles in read_labels)
    has_index2 = any(label == "Index2Cycles" for label, _cycles in read_labels)

    with open(path, "wt", newline="") as f_out:
        writer = csv.writer(f_out, lineterminator="\n")
        writer.writerow(["[Header]"])
        writer.writerow(["FileFormatVersion", "2"])
        writer.writerow([])

        writer.writerow(["[Reads]"])
        for label in ("Read1Cycles", "Read2Cycles", "Index1Cycles", "Index2Cycles"):
            for read_label, cycles in read_labels:
                if read_label == label:
                    writer.writerow([read_label, cycles])
        writer.writerow([])

        writer.writerow(["[BCLConvert_Settings]"])
        writer.writerow(["FastqCompressionFormat", "gzip"])
        writer.writerow([])

        writer.writerow(["[BCLConvert_Data]"])
        data_header = ["Lane", "Sample_ID"]
        if has_index1:
            data_header.append("Index")
        if has_index2:
            data_header.append("Index2")
        data_header += ["Sample_Project", "OverrideCycles"]
        if has_index1:
            data_header.append("BarcodeMismatchesIndex1")
        if has_index2:
            data_header.append("BarcodeMismatchesIndex2")
        writer.writerow(data_header)

        for sample in samples:
            if _is_undetermined_sample(sample):
                continue
            index, index2 = _sample_indexes(sample, exclude_index)
            for lane in sample["lanes"]:
                data_row = [lane, sample["sample_id"]]
                if has_index1:
                    data_row.append(index)
                if has_index2:
                    data_row.append(index2)
                data_row += [sample.get("project", ""), override_cycles]
                if has_index1:
                    data_row.append("1")
                if has_index2:
                    data_row.append("1")
                writer.writerow(data_row)
