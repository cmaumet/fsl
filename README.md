# fsl
FSL extension for voxel-wise thresholding

##### How to run make
```
$ cd /usr/local/fsl_MODIFIED/src/cluster
$ make
```

##### Voxel-wise corrected
<pre>
$ /usr/local/fsl_MODIFIED/src/cluster/cluster -i thresh_zstat1 -c stats/cope1 -t 2.3 -d 0.388188 --volume=45089 --othresh=thresh_zstat1_2 -o cluster_mask_zstat1 --connectivity=26 --olmax=lmax_zstat1.txt --scalarname=Z -p 0.05 --voxthresh
Cluster Index   Voxels  P   -log10(P)   Z-MAX   Z-MAX X (vox)   Z-MAX Y (vox)   Z-MAX Z (vox)   Z-COG X (vox)   Z-COG Y (vox)   Z-COG Z (vox)   COPE-MAX    COPE-MAX X (vox)    COPE-MAX Y (vox)    COPE-MAX Z (vox)    COPE-MEAN
2   1217    0.000829    3.08    5.77    41  16  13  34  14.7    14  1.37e+03    40  16  13  269
1   518 0.00487 2.31    5.44    41  41  25  43.5    40.4    24  864 47  38  21  337
</pre>

##### Cluster-wise corrected (same as before)
Original FSL call
```
$ /usr/local/fsl/bin/cluster -i thresh_zstat1 -c stats/cope1 -t 2.3 -d 0.388188 --volume=45089 --othresh=thresh_zstat1_2 -o cluster_mask_zstat1 --connectivity=26  --olmax=lmax_zstat1.txt --scalarname=Z -p 0.05
Cluster Index   Voxels  P   -log10(P)   Z-MAX   Z-MAX X (vox)   Z-MAX Y (vox)   Z-MAX Z (vox)   Z-COG X (vox)   Z-COG Y (vox)   Z-COG Z (vox)   COPE-MAX    COPE-MAX X (vox)    COPE-MAX Y (vox)    COPE-MAX Z (vox)    COPE-MEAN
4   1217    3.61e-24    23.4    5.77    41  16  13  34  14.7    14  1.37e+03    40  16  13  269
3   518 4.57e-13    12.3    5.44    41  41  25  43.5    40.4    24  864 47  38  21  337
2   114 0.000717    3.14    4.6 47  18  18  47.1    19.2    19.7    474 48  17  19  241
1   79  0.00995 2   4.59    33  38  31  32.3    39.2    30.9    370 32  38  31  228
```

Updated FSL call
```
Camilles-Computer:fmri_test.feat $ /usr/local/fsl_MODIFIED/src/cluster/cluster -i thresh_zstat1 -c stats/cope1 -t 2.3 -d 0.388188 --volume=45089 --othresh=thresh_zstat1_2 -o cluster_mask_zstat1 --connectivity=26  --olmax=lmax_zstat1.txt --scalarname=Z -p 0.05
Cluster Index   Voxels  P   -log10(P)   Z-MAX   Z-MAX X (vox)   Z-MAX Y (vox)   Z-MAX Z (vox)   Z-COG X (vox)   Z-COG Y (vox)   Z-COG Z (vox)   COPE-MAX    COPE-MAX X (vox)    COPE-MAX Y (vox)    COPE-MAX Z (vox)    COPE-MEAN
4   1217    3.61e-24    23.4    5.77    41  16  13  34  14.7    14  1.37e+03    40  16  13  269
3   518 4.57e-13    12.3    5.44    41  41  25  43.5    40.4    24  864 47  38  21  337
2   114 0.000717    3.14    4.6 47  18  18  47.1    19.2    19.7    474 48  17  19  241
1   79  0.00995 2   4.59    33  38  31  32.3    39.2    30.9    370 32  38  31  228
```