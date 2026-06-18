import subprocess
import shutil
import pandas as pd
from pathlib import Path

import hashlib
import gzip
import glob


def md5sum(filename):
    with open(filename, "rb") as f:
        return hashlib.file_digest(f, "md5").hexdigest()


def compare_fastq_files(expected_files, observed_dir):
    print("\nFASTQ md5 comparison")
    print(f"{'file':<35} {'expected':<32} {'observed':<32} status")
    print("-" * 110)

    failures = []
    for expected_file in expected_files:
        observed_file = observed_dir / expected_file.name
        if not observed_file.exists():
            failures.append((expected_file.name, "missing"))
            print(f"{expected_file.name:<35} {'-':<32} {'missing':<32} FAIL")
            continue

        expected_md5 = md5sum(expected_file)
        observed_md5 = md5sum(observed_file)
        status = "OK" if expected_md5 == observed_md5 else "FAIL"
        print(f"{expected_file.name:<35} {expected_md5:<32} {observed_md5:<32} {status}")
        if status != "OK":
            failures.append((expected_file.name, "md5 mismatch"))

    if failures:
        print(f"\nFAILED: {len(failures)} FASTQ file(s) differ")
        for filename, reason in failures:
            print(f"  - {filename}: {reason}")
        raise SystemExit(1)

    print(f"\nOK: {len(expected_files)} FASTQ file(s) match")


def read_fastq_records(filename):
    records = []
    with gzip.open(filename, "rt") as f_in:
        while True:
            header = f_in.readline()
            if header == "":
                break
            sequence = f_in.readline()
            plus = f_in.readline()
            quality = f_in.readline()
            if not sequence or not plus or not quality:
                raise ValueError(f"Malformed FASTQ record in {filename}")
            records.append((header, sequence, plus, quality))
    return records


def write_sorted_fastq(input_files, output_file):
    records = []
    for input_file in input_files:
        records.extend(read_fastq_records(input_file))
    records.sort(key=lambda record: record[0])

    with open(output_file, "wb") as f_out:
        with gzip.GzipFile(filename="", mode="wb", fileobj=f_out, mtime=0, compresslevel=9) as gz_out:
            for record in records:
                gz_out.write("".join(record).encode())
    with open(output_file, "r+b") as f_out:
        f_out.seek(9)
        f_out.write(b"\x03")


def normalize_bclconvert_fastqs(input_pattern, output_file):
    input_files = sorted(glob.glob(input_pattern))
    write_sorted_fastq(input_files, output_file)


def check_requirements():
    if shutil.which("fastq2bcl") is None:
        print("no executable found for command 'fastq2bcl'")
        exit(1)
    if shutil.which("bcl-convert") is None:
        print("no executable found for command 'bcl-convert'")
        exit(1)


check_requirements()

tmp_test_dir: Path = Path("tmp_testing")
if tmp_test_dir.exists():
    shutil.rmtree(tmp_test_dir)
tmp_test_dir.mkdir(exist_ok=True, parents=True)
test_data_dir: Path = Path("data/test")
for test_case in [
    "10_input_samplesheet",
    "11_no_undetermined_specified",
    "12_separate_undetermined",
]:
    current_test_data_dir: Path = Path(tmp_test_dir, test_case)

    ### Run fastq2bcl ###
    run_dir = Path(current_test_data_dir)
    run_dir.mkdir(exist_ok=True, parents=True)
    input_samplesheet = Path(test_data_dir, test_case, "input_samplesheet.csv")
    subprocess.run(
        [
            "fastq2bcl",
            "--input_samplesheet",
            str(input_samplesheet),
            "--sample-sheet-format",
            "bcl-convert",
            "--outdir",
            str(run_dir),
        ],
        check=True,
        # stdout=subprocess.DEVNULL,
        # stderr=subprocess.STDOUT
    )

    ### Run bcl convert ###
    bclconvert_output = Path(current_test_data_dir, "bclconvert")

    bclconvert_output.mkdir(exist_ok=True, parents=True)
    subprocess.run(
        [
            "bcl-convert",
            "--force",
            "--output-directory",
            str(bclconvert_output),
            "--bcl-input-directory",
            f"{str(run_dir)}/100101_run_0001_ABCD",
            "--sample-sheet",
            f"{str(run_dir)}/100101_run_0001_ABCD/SampleSheet.csv",
        ],
        check=True,
        # stdout=subprocess.DEVNULL,
        # stderr=subprocess.STDOUT
    )

    ### Rezip data ###
    input_samplesheet = pd.read_csv(input_samplesheet)
    input_samplesheet["lane"] = input_samplesheet["lane"].str.split(";")

    bclconvert_data = Path(current_test_data_dir, "bclconvert_data")
    bclconvert_data.mkdir(exist_ok=True, parents=True)

    # Rezip Undetermined files
    for lane in input_samplesheet.explode("lane")["lane"].unique():
        input_R1_files = f"{str(bclconvert_output)}/Undetermined_*_L00{lane}_R1_*fastq.gz"
        output_R1_fastq = f"{str(bclconvert_data)}/Undetermined_L00{lane}_R1.fastq.gz"

        input_R2_files = f"{str(bclconvert_output)}/Undetermined_*_L00{lane}_R2_*fastq.gz"
        output_R2_fastq = f"{str(bclconvert_data)}/Undetermined_L00{lane}_R2.fastq.gz"

        for input_files, output_file in [(input_R1_files, output_R1_fastq), (input_R2_files, output_R2_fastq)]:
            normalize_bclconvert_fastqs(input_files, output_file)

    for row in input_samplesheet.itertuples():
        if "undetermined" in str(row.project).lower() or "undetermined" in str(row.sample).lower():
            continue

        sample = row.sample
        input_R1_files = [f"{str(bclconvert_output)}/{sample}_*_L00{lane}_R1_*fastq.gz" for lane in row.lane]
        output_R1_fastq = f"{str(bclconvert_data)}/{sample}_R1.fastq.gz"
        input_R2_files = [f"{str(bclconvert_output)}/{sample}_*_L00{lane}_R2_*fastq.gz" for lane in row.lane]
        output_R2_fastq = f"{str(bclconvert_data)}/{sample}_R2.fastq.gz"

        for input_files, output_file in [(input_R1_files, output_R1_fastq), (input_R2_files, output_R2_fastq)]:
            matched_files = []
            for input_pattern in input_files:
                matched_files.extend(sorted(glob.glob(input_pattern)))
            write_sorted_fastq(matched_files, output_file)

    expected_output = Path(test_data_dir, test_case, "expected_output")
    expected_files = sorted(expected_output.glob("*.fastq.gz"))
    compare_fastq_files(expected_files, bclconvert_data)
