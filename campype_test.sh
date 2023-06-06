wget -nc --directory-prefix test --input-file test/test_files.txt

cat test/test_config.py > campype_config.py
cat test/input_files_test.csv > input_files.csv

bash -i campype

# Remove the downloaded files
rm test/*.fastq.gz
