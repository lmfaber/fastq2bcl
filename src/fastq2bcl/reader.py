import logging
import gzip
from fastq2bcl.parser import parse_seqdesc_fields
from Bio import SeqIO
from rich import print

_logger = logging.getLogger(__name__)


def iter_raw_fastq_records(fastq_file):
    """
    Stream FASTQ records as ``(description, id, sequence, qualities)``.

    ``qualities`` are returned as Phred integers for compatibility with the
    writer API, but this avoids BioPython record allocation in the hot loop.
    """
    with gzip.open(fastq_file, "rt") as fastq_fh:
        while True:
            header = fastq_fh.readline()
            if header == "":
                break
            sequence = fastq_fh.readline().rstrip("\n\r")
            plus = fastq_fh.readline()
            quality = fastq_fh.readline().rstrip("\n\r")
            if not sequence or not plus or not quality:
                raise ValueError(f"Malformed FASTQ record in {fastq_file}")
            if not header.startswith("@"):
                raise ValueError(f"Malformed FASTQ header in {fastq_file}: {header.rstrip()}")
            description = header[1:].rstrip("\n\r")
            record_id = description.split(None, 1)[0]
            yield description, record_id, sequence, [ord(char) - 33 for char in quality]


def read_first_record(fastq_file):
    """
    Validate fastq.gz r1 file and extract first read
    """
    _logger.info(f"Opening gz file {fastq_file}")
    with gzip.open(fastq_file, "rt") as fastq_fh:
        return next(SeqIO.parse(fastq_fh, "fastq"))


def get_file_handlers(r1, r2, i1, i2):
    """
    Return list of FH
    """
    files_fh = [gzip.open(r1, "rt")]
    if not i1 == None:
        files_fh.append(gzip.open(i1, "rt"))
    if not i2 == None:
        files_fh.append(gzip.open(i2, "rt"))
    if not r2 == None:
        files_fh.append(gzip.open(r2, "rt"))

    return files_fh


def get_fastq_paths(r1, r2, i1, i2):
    files = [r1]
    if i1 != None:
        files.append(i1)
    if i2 != None:
        files.append(i2)
    if r2 != None:
        files.append(r2)
    return files


def get_mask_from_files(r1, r2, i1, i2, exclude_umi, exclude_index):
    """
    Build a mask string using seq length. In case of index and/or UMI in R1 sequence description, write this length to the Index mask
    """
    record_1 = read_first_record(r1)
    seq_fields = parse_seqdesc_fields(record_1.description)
    index_1_bases = 0

    # Write R1 mask
    mask = f"{len(record_1.seq)}N"

    # check errors on index for R1
    if seq_fields["index"] != "1" and not exclude_index:
        if i1 != None or i2 != None:
            raise ValueError("Usage of index from sequence desc and I1 and I2 files at the same time is not supported")
        # continue and write to index I1 length TODO I2 for double index
        print(f"LENGTH INDEX {seq_fields['index']}")
        index_1_bases += len(seq_fields["index"])

    # check errors on UMI for R1
    if seq_fields["UMI"] != None and not exclude_umi:
        if i1 != None or i2 != None:
            raise ValueError("Usage of UMI from sequence desc and I1 and I2 files at the same time is not supported")
        # continue and write to index I1
        index_1_bases += len(seq_fields["UMI"])

    # Write index I1 based on UMI and index
    if index_1_bases > 0:
        mask += f"{index_1_bases}Y"

    # Write indexes
    if i1 != None:
        index_1 = read_first_record(i1)
        mask += f"{len(index_1.seq)}Y"

    if i2 != None:
        index_2 = read_first_record(i2)
        mask += f"{len(index_2.seq)}Y"

    # Write R2 record
    if r2 != None:
        record_2 = read_first_record(r2)
        # finally add R2 to mask
        mask += f"{len(record_2.seq)}N"

    return mask


def iter_fastq_records(r1, r2, i1, i2, exclude_umi, exclude_index):
    """
    Stream fastq files R1-R2 with I1 and I2 and yield only the data we need.

    Each yielded item is ``((sequence, quality), (x_pos, y_pos))``.
    """
    seq_iterators = [iter_raw_fastq_records(path) for path in get_fastq_paths(r1, r2, i1, i2)]

    for r1_description, record_id, record_seq, record_qual in seq_iterators[0]:
        try:
            opt_data = [next(iterator) for iterator in seq_iterators[1:]]
        except StopIteration:
            raise ValueError("FASTQ files do not contain the same number of records")

        record_fields = parse_seqdesc_fields(r1_description)

        if not exclude_index and record_fields["index"] != "1":
            _logger.info(f"Reading index field: {record_fields['index']}")
            record_seq += record_fields["index"]
            record_qual += [40] * len(record_fields["index"])

        # the UMI is quality MAX (40)
        if not exclude_umi and record_fields["UMI"] != None:
            _logger.info(f"Reading umi field: {record_fields['UMI']}")
            record_seq += record_fields["UMI"]
            record_qual += [40] * len(record_fields["UMI"])

        for _description, opt_record_id, opt_seq, opt_qual in opt_data:
            if opt_record_id != record_id:
                raise ValueError(f"Seq ID mismatch for record {opt_record_id} R1 is {record_id}")
            record_seq += opt_seq
            record_qual += opt_qual

        yield (record_seq, record_qual), (
            record_fields["x_pos"],
            record_fields["y_pos"],
        )

    for iterator in seq_iterators[1:]:
        try:
            next(iterator)
        except StopIteration:
            continue
        raise ValueError("FASTQ files do not contain the same number of records")


def iter_fastq_cycle_data(cycle, r1, r2, i1, i2, exclude_umi, exclude_index):
    """
    Stream ``(base, quality)`` values for a single BCL cycle.
    """
    for (basecalls, qualscores), _position in iter_fastq_records(r1, r2, i1, i2, exclude_umi, exclude_index):
        if cycle >= len(basecalls):
            yield ("N", 0)
        else:
            yield (basecalls[cycle], qualscores[cycle])


def count_fastq_records(r1, r2, i1, i2, exclude_umi, exclude_index):
    """
    Count records without retaining FASTQ data in memory.
    """
    return sum(1 for _record in iter_fastq_records(r1, r2, i1, i2, exclude_umi, exclude_index))


def read_fastq_files(r1, r2, i1, i2, exclude_umi, exclude_index):
    """
    Compatibility helper that materializes FASTQ data in memory.

    Production code should use ``iter_fastq_records`` instead.
    """
    # return a list of tuple with seq, qual
    # and a list of tuple for pos with x and y
    # SINGLE R1
    # sequences = [('AAAA',1111)]
    # positions = [(1,1)]
    #
    # in case of multiple files R1-R2:
    # PAIR R1-R2
    # sequences = [('AAAABBBB',11111111)]
    # positions = [(1,1)]
    #
    # I need way to handle multiple files and merge them in a single with exitstack
    # Ref https://docs.python.org/3/library/contextlib.html#contextlib.ExitStack

    # output Lists
    sequences = []
    positions = []

    for sequence, position in iter_fastq_records(r1, r2, i1, i2, exclude_umi, exclude_index):
        positions.append(position)
        sequences.append(sequence)

    return (sequences, positions)
