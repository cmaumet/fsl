# fsl
FSL extension for voxel-wise thresholding

#### Make
```
FSLDIRUP=/usr/local/fsl_MODIFIED
cd ${FSLDIRUP}/src/cluster
make
```

#### Examples of usage for the updated `cluster` command
##### Voxel-wise corrected
<pre>
${FSLDIR}/bin/cluster -i thresh_zstat1 -c stats/cope1 -t 3.0902 -d 0.700797 --volume=38344 -r 6.60989 
--othresh=thresh_zstat1_2 -o cluster_mask_zstat1 --connectivity=26 --olmax=lmax_zstat1.txt --scalarname=Z 
<b>--voxthresh</b>

Cluster Index   Voxels  Z-MAX   Z-MAX X (vox)   Z-MAX Y (vox)   Z-MAX Z (vox)   
Z-COG X (vox)   Z-COG Y (vox)   Z-COG Z (vox)   COPE-MAX    COPE-MAX X (vox)    
COPE-MAX Y (vox)    COPE-MAX Z (vox)    COPE-MEAN
20  27  7.47    19  40  18  18.8    40.3    17.9    420 18  42  19  243
19  13  7.49    16  25  20  15.9    25  19.9    222 16  25  20  159
18  11  7.26    49  34  20  48.8    34.5    19.4    194 50  34  20  153
17  10  6.76    21  10  11  21.7    10.2    10.6    323 23  9   10  226
16  10  6.1 23  44  13  23.4    43.8    13.6    235 23  44  13  186
15  8   6.03    36  37  23  35.9    36.8    22.5    222 35  37  23  167
14  6   5.86    25  9   14  24.2    9.49    13.7    323 24  9   14  236
13  6   6.52    44  33  22  44  33.2    22  200 44  33  22  166
12  6   6.3 23  35  22  22.9    35  21.5    243 23  35  22  175
11  5   6.19    15  20  22  16  20.2    21.2    241 15  21  22  181
10  3   5.58    20  16  20  20.3    16  20.3    187 20  16  20  170
9   2   5.18    21  14  9   21.5    14  9   198 22  14  9   197
8   2   5.09    15  24  13  15.5    24  13  157 15  24  13  148
7   2   5.43    16  42  16  16  42.5    16  238 16  42  16  198
6   2   5.93    52  37  17  52  37.5    17  141 52  37  17  139
5   2   5.86    13  23  18  13  22.5    18  182 13  22  18  170
4   2   5.32    31  39  22  31.5    39  22  213 32  39  22  203
3   2   5.35    17  19  24  17.5    18.5    23.5    409 18  18  23  268
2   1   5.99    13  23  22  13  23  22  117 13  23  22  117
1   1   5.26    21  41  23  21  41  23  221 21  41  23  221

</pre>

##### Voxel-wise uncorrected
<pre>
${FSLDIR}/bin/cluster -i thresh_zstat1 -c stats/cope1 -t 3.0902 -d 0.700797 --volume=38344 -r 6.60989 
--othresh=thresh_zstat1_2 -o cluster_mask_zstat1 --connectivity=26 --olmax=lmax_zstat1.txt --scalarname=Z 
<b>--voxuncthresh</b>

Cluster Index   Voxels  Z-MAX   Z-MAX X (vox)   Z-MAX Y (vox)   Z-MAX Z (vox)   
Z-COG X (vox)   Z-COG Y (vox)   Z-COG Z (vox)   COPE-MAX    COPE-MAX X (vox)    
COPE-MAX Y (vox)    COPE-MAX Z (vox)    COPE-MEAN
147 795 7.47    19  40  18  24.7    41.7    18.4    632 28  56  17  168
146 243 7.49    16  25  20  16.4    21.6    18.8    428 17  19  23  135
145 161 4.56    42  43  14  39.3    44.2    15.1    479 33  36  11  199
144 152 7.26    49  34  20  47.7    35  19.7    218 52  38  16  117
143 66  6.76    21  10  11  22.3    10.9    11.3    323 23  9   10  171
142 65  5.04    46  23  20  45.3    20.7    20  242 41  16  21  131
141 42  4.95    49  46  17  48.1    45.7    17.9    357 47  47  17  184
140 27  4.03    36  7   6   35.2    8.44    3.52    200 31  8   4   136
139 24  4.41    27  18  17  28.1    14.6    17.6    263 28  13  18  162
138 17  4.77    47  21  15  47.7    20.9    14.3    200 48  19  14  125
137 16  4.71    39  8   12  40.6    8.91    11.2    224 39  8   12  142
136 15  4.67    28  13  0   25.4    14  0   415 25  14  0   226
135 15  4.74    40  14  0   38.2    13.8    0   591 40  14  0   298
134 14  4.48    47  16  10  47.1    15.5    9.86    204 46  13  10  151
133 13  3.91    44  14  8   43.6    16.2    6.95    533 43  16  7   295
[...]

</pre>

##### Cluster-wise corrected (same as before)
```
${FSLDIR}/bin/cluster -i thresh_zstat1 -c stats/cope1 -t 2.3 -d 0.700797 --volume=38344 -r 6.60989 
--othresh=thresh_zstat1_2 -o cluster_mask_zstat1 --connectivity=26  --olmax=lmax_zstat1.txt --scalarname=Z -p 0.05

Cluster Index   Voxels  P   -log10(P)   Z-MAX   Z-MAX X (vox)   Z-MAX Y (vox)  Z-MAX Z (vox)   
Z-COG X (vox)   Z-COG Y (voxZ-COG Z (vox)   COPE-MAX    COPE-MAX X (vox)   COPE-MAX Y (vox)    
COPE-MAX Z (vox)    COPE-MEAN
7   3130    0   68.9    7.47    19  40  18  32.5    39.9    18.1    632 28  56  17  139
6   521 6.02e-20    19.2    7.49    16  25  20  16  22.4    17.5    428 17  19  23  110
5   167 1.7e-08 7.77    4.71    39  8   12  43.9    14.9    9.45    533 43  16  7   146
4   99  1.58e-05    4.8 6.76    21  10  11  23.1    10.5    11.5    323 23  9   10  143
3   43  0.017   1.77    4.03    36  7   6   35.5    8.35    3.75    200 31  8   4   128
2   41  0.0227  1.64    4.41    27  18  17  27.9    14.6    17.6    263 28  13  18  131
1   40  0.0263  1.58    4.77    47  21  15  48  20.8    14  200 48  19  14  102
```
