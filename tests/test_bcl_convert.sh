if ! command -v fastq2bcl >/dev/null 2>&1; then
    echo "fastq2bcl could not be found"
    exit 1
fi
rm -r 100101_run_0001_ABCD
fastq2bcl --input_samplesheet data/test/10_input_samplesheet/input_samplesheet.csv --sample-sheet-format bcl-convert

if ! command -v bcl-convert >/dev/null 2>&1; then
    echo "bcl-convert could not be found"
    exit 1
fi
rm -r test10
bcl-convert --force --output-directory test10 --bcl-input-directory 100101_run_0001_ABCD
--bcl-sampleproject-subdirectories true --sample-sheet 100101_run_0001_ABCD/SampleSheet.csv
if [ "$?" -gt 1 ]; then
    echo "bcl-convert failed"
    exit 1
fi
