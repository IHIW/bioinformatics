# filename: Run_GenerateReferences.sh
# date: 2020-11-17
# version: 1.0
# author: Ben Matern <B.M.Matern@umcutrecht.nl>

# description:
# 	This Script runs the Reference Sequence Generator. This chooses valid reference sequences to be used in HML documents for the 18th Workshop.

# My biopython etc. packages are installed in a virtual environment. Change or delete this if your venv is different.
source ../../venv/bin/activate

# Iterate a list of IPD-IMGT/HLA releases, separated by spaces. Or just chose a single release.
#releases="3.43.0 3.42.0 3.41.0 3.40.0 3.39.0 3.38.0 3.37.0 3.36.0 3.35.0 3.34.0 3.33.0 3.32.0 3.31.0 3.30.0 3.29.0 3.28.0 3.27.0 3.26.0"
releases="3.43.0 3.42.0"
version="1.3"
for release in $releases; do
    echo $release
    python GenerateReferences.py --release=$release --output="../../reference_alleles/"$release"_catalog" --version=$version --supplementary
done

# Deactivating my virtual environment
deactivate
