## Developmental correspondence of juvenile stages across the locust, harlequin ladybird, and diamondback moth

Hang Zhou1†*, Runguo Shu1†, Chaowei Zhang1, Yiqi Xiao1, Dong Jing1, Jiejing Tang1, Zixiong Cao2, Xi Chen1, Yang Mei1, Fei Li1*

Code version:1.0

The code in "Code S1.md" is primarily used for standardization of expression levels and clustering analysis. The expression TPM matrix of different insects is available in Data S4.

### Trim_galore
'
trim_galore -j 20 -q 20 --phred33 --stringency 3 --length 25 -e 0.1 --paired $dir/cmp/01raw_data/$fq1 $dir/cmp/01raw_data/$fq2 --gzip -o out_dir
'
