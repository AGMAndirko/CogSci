#!/usr/bin/env bash

# Take the Peyregne et al data and create a pseudo-vcf
sed 's/,/\t/g' 'Peyregne 2017.vcf' | awk '{print $1, $2, "-", $4, $5}' | tail -n +2 | awk '$5=! "NA" {print $0}' | sed 's/ /\t/g' > formatted_peyregne.vcf
 
# from then on, follow ExPecto's intructions, detailed here: https://github.com/FunctionLab/ExPecto
git clone https://github.com/FunctionLab/ExPecto.git
mv formatted_peyregne.vcf Expecto/
cd ExPecto
sh download_resources.sh; tar xf resources_20190807.tar.gz

# Download BEDOPS (https://bedops.readthedocs.io/en/latest/) and run:
closest-features --delim '\t' --closest --dist <(awk '{printf $1"\t"$2-1"\t"$2"\n"}' formatted_peyregne.vcf|sed s/chr//g|sed s/^/chr/g|sort-bed - ) ./resources/geneanno.pc.sorted.bed > formatter_peyregne.bed.sorted.bed.closestgene

# Run ExPecto after satifying whatever dependencies you might miss (rhis might take a long timetime):
python chromatin.py formatted_peyregne.vcf 
python predict.py --coorFile formatted_peyregne.vcf  --geneFile formatted_peyregne.vcf.bed.sorted.bed.closestgene --snpEffectFilePattern formatted_peyregne.vcf.shift_SHIFT.diff.h5 --modelList ./resources/modellist --output output.csv
