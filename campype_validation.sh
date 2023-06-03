wget -nc --directory-prefix campype_validation --input-file campype_validation/validation_files.txt

cat campype_validation/validation_config.py > campype_config.py
cat campype_validation/input_files_validation.csv > input_files.csv

bash -i campype