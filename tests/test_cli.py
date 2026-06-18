import struct
import pytest
import unittest.mock

from fastq2bcl.cli import (
    main,
    mock_run_id,
    fastq2bcl,
    set_mask,
    run,
    parse_input_samplesheet,
)

__author__ = "Davide Rambaldi"
__copyright__ = "Davide Rambaldi"
__license__ = "MIT"

test_fields = {
    "instrument": "M11111",
    "run_number": "222",
    "flowcell_id": "000000000-K9H97",
    "lane": "1",
    "tile": "1101",
    "x_pos": "19304",
    "y_pos": "1328",
    "UMI": "AAACGGG",
    "read": "1",
    "is_filtered": "N",
    "control_number": "0",
    "index": "1",
}

expected_mask_110N = [{"cycles": "110", "index": "N", "id": "1"}]
expected_mask_110N10Y10Y110N = [
    {"cycles": "110", "index": "N", "id": "1"},
    {"cycles": "10", "index": "Y", "id": "2"},
    {"cycles": "10", "index": "Y", "id": "3"},
    {"cycles": "110", "index": "N", "id": "4"},
]


def filter_cluster_count(rundir, lane):
    with open(rundir / f"Data/Intensities/BaseCalls/L00{lane}/s_{lane}_1101.filter", "rb") as f_in:
        f_in.seek(8)
        return struct.unpack("<I", f_in.read(4))[0]


def test_run():
    """run CLI test"""
    with pytest.raises(SystemExit):
        run()


def test_main_usage(capsys, tmpdir):
    """CLI Tests"""
    # capsys is a pytest fixture that allows asserts against stdout/stderr
    # https://docs.pytest.org/en/stable/capture.html
    main(["-o", str(tmpdir), "-r1", "data/test/01_single/test_single.fastq.gz"])
    captured = capsys.readouterr()
    assert "100101_M11111_0222_000000000-K9H97" in captured.out


def test_parse_input_samplesheet_with_comma_lane_list(tmpdir):
    samplesheet = tmpdir / "input.csv"
    samplesheet.write(
        "lane,r1,r2,i1,i2,project,sample,extra\n"
        "1;2;3;4,5,6,7,8,data/test/07_pair/R1.fastq.gz,data/test/07_pair/R2.fastq.gz,,,ProjectA,SampleA,ignored\n"
        "8,data/test/01_single/test_single.fastq.gz,,,\n"
    )

    samples = parse_input_samplesheet(str(samplesheet))

    assert samples[0]["lanes"] == [1, 2, 3, 4, 5, 6, 7, 8]
    assert samples[0]["r1"].name == "R1.fastq.gz"
    assert samples[0]["r2"].name == "R2.fastq.gz"
    assert samples[0]["sample_id"] == "SampleA"
    assert samples[0]["sample_name"] == "SampleA"
    assert samples[0]["project"] == "ProjectA"
    assert samples[1]["lanes"] == [8]
    assert samples[1]["sample_id"] == "test_single"
    assert samples[1]["project"] == ""


def test_parse_input_samplesheet_rejects_lane_above_8(tmpdir):
    samplesheet = tmpdir / "input.csv"
    samplesheet.write("lane,r1,r2,i1,i2\n" "9,data/test/01_single/test_single.fastq.gz,,,\n")

    with pytest.raises(ValueError, match="lanes must be between 1 and 8"):
        parse_input_samplesheet(str(samplesheet))


def test_fastq2bcl_with_input_samplesheet(tmpdir):
    run_id, rundir, seqdesc_fields, mask_string = fastq2bcl(
        str(tmpdir),
        None,
        input_samplesheet="data/test/10_input_samplesheet/input_samplesheet.csv",
    )

    assert run_id == "100101_run_0001_ABCD"
    assert mask_string == "309N8Y309N"
    assert 'LaneCount="8"' in (rundir / "RunInfo.xml").read_text()
    assert (rundir / "Data/Intensities/BaseCalls/L001/s_1_1101.filter").is_file()
    assert (rundir / "Data/Intensities/BaseCalls/L008/s_8_1101.filter").is_file()
    assert (rundir / "Data/Intensities/L008/s_8_1101.locs").is_file()
    samplesheet = (rundir / "SampleSheet.csv").read_text()
    assert "1,sample1,sample1,," in samplesheet
    assert "2,sample2,sample2,," in samplesheet
    assert "8,sample2,sample2,," in samplesheet
    assert "8,sample3,sample3,," in samplesheet
    assert "undetermined" not in samplesheet.lower()
    assert "empty_R1" not in samplesheet

    with open(rundir / "Data/Intensities/BaseCalls/L001/s_1_1101.filter", "rb") as f_in:
        f_in.seek(8)
        assert struct.unpack("<I", f_in.read(4))[0] == 1
    with open(rundir / "Data/Intensities/BaseCalls/L002/s_2_1101.filter", "rb") as f_in:
        f_in.seek(8)
        assert struct.unpack("<I", f_in.read(4))[0] == 2
    with open(rundir / "Data/Intensities/BaseCalls/L008/s_8_1101.filter", "rb") as f_in:
        f_in.seek(8)
        assert struct.unpack("<I", f_in.read(4))[0] == 2


def test_fastq2bcl_with_bcl_convert_sample_sheet(tmpdir):
    run_id, rundir, seqdesc_fields, mask_string = fastq2bcl(
        str(tmpdir),
        None,
        input_samplesheet="data/test/10_input_samplesheet/input_samplesheet.csv",
        sample_sheet_format="bcl-convert",
    )

    samplesheet = (rundir / "SampleSheet.csv").read_text()

    assert run_id == "100101_run_0001_ABCD"
    assert mask_string == "309N8Y309N"
    assert "[BCLConvert_Settings]\nFastqCompressionFormat,gzip" in samplesheet
    assert ("Lane,Sample_ID,Index,Sample_Project,OverrideCycles,BarcodeMismatchesIndex1") in samplesheet
    assert "Read1Cycles,309\nRead2Cycles,309\nIndex1Cycles,8" in samplesheet
    assert "1,sample1,CGCGCGCG,project1,Y309;I8;Y309,1" in samplesheet
    assert "2,sample2,AACCACTA,project2,Y309;I8;Y309,1" in samplesheet
    assert "8,sample2,AACCACTA,project2,Y309;I8;Y309,1" in samplesheet
    assert "8,sample3,ATATATAT,project3,Y309;I8;Y309,1" in samplesheet
    assert "1,sample2,AACCACTA,project2,Y309;I8;Y309,1" not in samplesheet
    assert "undetermined" not in samplesheet.lower()
    assert "empty_R1" not in samplesheet


def test_fastq2bcl_with_unspecified_undetermined_reads(tmpdir):
    run_id, rundir, seqdesc_fields, mask_string = fastq2bcl(
        str(tmpdir),
        None,
        input_samplesheet="data/test/11_no_undetermined_specified/input_samplesheet.csv",
        sample_sheet_format="bcl-convert",
    )

    samplesheet = (rundir / "SampleSheet.csv").read_text()

    assert run_id == "100101_run_0001_ABCD"
    assert mask_string == "309N8Y309N"
    assert "1,sample1,CGCGCGCG,project1,Y309;I8;Y309,1" in samplesheet
    assert "2,sample2,AACCACTA,project2,Y309;I8;Y309,1" in samplesheet
    assert "8,sample3,ATATATAT,project3,Y309;I8;Y309,1" in samplesheet
    assert "undetermined" not in samplesheet.lower()
    assert filter_cluster_count(rundir, 1) == 2
    assert filter_cluster_count(rundir, 2) == 2
    assert filter_cluster_count(rundir, 8) == 3


def test_fastq2bcl_with_separate_undetermined_reads(tmpdir):
    run_id, rundir, seqdesc_fields, mask_string = fastq2bcl(
        str(tmpdir),
        None,
        input_samplesheet="data/test/12_separate_undetermined/input_samplesheet.csv",
        sample_sheet_format="bcl-convert",
    )

    samplesheet = (rundir / "SampleSheet.csv").read_text()

    assert run_id == "100101_run_0001_ABCD"
    assert mask_string == "309N8Y309N"
    assert "1,sample1,CGCGCGCG,project1,Y309;I8;Y309,1" in samplesheet
    assert "2,sample2,AACCACTA,project2,Y309;I8;Y309,1" in samplesheet
    assert "8,sample3,ATATATAT,project3,Y309;I8;Y309,1" in samplesheet
    assert "undetermined" not in samplesheet.lower()
    assert filter_cluster_count(rundir, 1) == 2
    assert filter_cluster_count(rundir, 2) == 2
    assert filter_cluster_count(rundir, 8) == 3


def test_fastq2bcl_with_dual_index_bcl_convert_sample_sheet(tmpdir):
    run_id, rundir, seqdesc_fields, mask_string = fastq2bcl(
        str(tmpdir),
        "data/test/05_multi_pair_double_index/R1.fastq.gz",
        "data/test/05_multi_pair_double_index/R2.fastq.gz",
        "data/test/05_multi_pair_double_index/RIndex1.fastq.gz",
        "data/test/05_multi_pair_double_index/RIndex2.fastq.gz",
        exclude_index=True,
        sample_sheet_format="bcl-convert",
    )

    samplesheet = (rundir / "SampleSheet.csv").read_text()

    assert mask_string == "296N8Y8Y309N"
    assert (
        "Lane,Sample_ID,Index,Index2,Sample_Project,OverrideCycles," "BarcodeMismatchesIndex1,BarcodeMismatchesIndex2"
    ) in samplesheet
    assert "Read1Cycles,296\nRead2Cycles,309\nIndex1Cycles,8\nIndex2Cycles,8" in samplesheet
    assert "1,R1,AACCACTA,AACCACTA,,Y296;I8;I8;Y309,1,1" in samplesheet


def test_multithread_usage(capsys, tmpdir):
    """CLI Tests"""
    # capsys is a pytest fixture that allows asserts against stdout/stderr
    # https://docs.pytest.org/en/stable/capture.html
    main(
        [
            "-o",
            str(tmpdir),
            "-r1",
            "data/test/01_single/test_single.fastq.gz",
            "-T",
            "16",
        ]
    )
    captured = capsys.readouterr()
    assert "100101_M11111_0222_000000000-K9H97" in captured.out


def test_multithread_usage_and_different_length(capsys, tmpdir):
    """CLI Test"""
    # capsys is a pytest fixture that allows asserts against stdout/stderr
    # https://docs.pytest.org/en/stable/capture.html
    main(
        [
            "-o",
            str(tmpdir),
            "-r1",
            "data/test/09_multi_pair_different_indexes/R1.fastq.gz",
            "-r2",
            "data/test/09_multi_pair_different_indexes/R2.fastq.gz",
            "-i1",
            "data/test/09_multi_pair_different_indexes/RIndex1.fastq.gz",
            "-i2",
            "data/test/09_multi_pair_different_indexes/RIndex2.fastq.gz",
            "-T",
            "16",
            "--exclude-index",
        ]
    )
    captured = capsys.readouterr()
    assert "100101_run_0001_ABCD" in captured.out


def test_mock_run_id():
    """Mock run id Tests"""
    assert mock_run_id(test_fields) == "100101_M11111_0222_000000000-K9H97"


def test_fastq2bcl(tmpdir):
    """Fastq2bcl main function Tests"""
    run_id, rundir, seqdesc_fields, mask_string = fastq2bcl(str(tmpdir), "data/test/01_single/test_single.fastq.gz")
    assert seqdesc_fields["flowcell_id"] == "000000000-K9H97"
    assert run_id == "100101_M11111_0222_000000000-K9H97"
    assert mask_string == "110N"


def test_fastq2bcl_with_mask(tmpdir):
    """Fastq2bcl main function Tests"""
    run_id, rundir, seqdesc_fields, mask_string = fastq2bcl(
        str(tmpdir),
        "data/test/07_pair/R1.fastq.gz",
        "data/test/07_pair/R2.fastq.gz",
        mask_string="309N309N",
    )
    assert seqdesc_fields["flowcell_id"] == "ABCD"
    assert run_id == "100101_run_0001_ABCD"
    assert mask_string == "309N309N"


def test_set_mask():
    """Test mask generation"""
    assert set_mask("110N") == expected_mask_110N
    assert set_mask("110N10Y10Y110N") == expected_mask_110N10Y10Y110N
    with pytest.raises(ValueError):
        set_mask(None)
    with pytest.raises(ValueError):
        set_mask("100")


def test_fastq2bcl_with_umi(tmpdir):
    """Fastq2bcl main function Tests with UMI"""
    run_id, rundir, seqdesc_fields, mask_string = fastq2bcl(
        str(tmpdir), "data/test/03_single_with_umi/single_with_umi.fastq.gz"
    )
    assert seqdesc_fields["UMI"] == "ACGTAGTAC"


def test_fastq2bcl_with_exclude_umi(tmpdir):
    """Fastq2bcl main function Tests with exclude-umi"""
    run_id, rundir, seqdesc_fields, mask_string = fastq2bcl(
        str(tmpdir),
        "data/test/03_single_with_umi/single_with_umi.fastq.gz",
        exclude_umi=True,
    )
    assert seqdesc_fields["UMI"] == "ACGTAGTAC"


def test_fastq2bcl_with_index(tmpdir):
    """Fastq2bcl main function Tests with INDEX"""
    run_id, rundir, seqdesc_fields, mask_string = fastq2bcl(
        outdir=str(tmpdir),
        r1="data/test/02_single_with_index/single.R1.fastq.gz",
    )
    assert seqdesc_fields["index"] == "AACCACTA"


def test_fastq2bcl_with_exclude_index(tmpdir):
    """Fastq2bcl main function Tests with exclude-index"""
    run_id, rundir, seqdesc_fields, mask_string = fastq2bcl(
        outdir=str(tmpdir),
        r1="data/test/02_single_with_index/single.R1.fastq.gz",
        exclude_index=True,
    )
    assert seqdesc_fields["index"] == "AACCACTA"
