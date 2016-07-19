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
${FSLDIRUP}/src/cluster/cluster -i thresh_zstat1 -c stats/cope1 -t 2.3 -d 0.388188 --volume=45089 --othresh=thresh_zstat1_2 -o cluster_mask_zstat1 --connectivity=26 --olmax=lmax_zstat1.txt --scalarname=Z -p 0.05 --voxthresh

Cluster Index   Voxels  P   -log10(P)   Z-MAX   Z-MAX X (vox)   Z-MAX Y (vox)   Z-MAX Z (vox)   Z-COG X (vox)   Z-COG Y (vox)   Z-COG Z (vox)   COPE-MAX    COPE-MAX X (vox)    COPE-MAX Y (vox)    COPE-MAX Z (vox)    COPE-MEAN
2   1217    0.000829    3.08    5.77    41  16  13  34  14.7    14  1.37e+03    40  16  13  269
1   518 0.00487 2.31    5.44    41  41  25  43.5    40.4    24  864 47  38  21  337
</pre>

##### Voxel-wise uncorrected
<pre>
${FSLDIRUP}/src/cluster/cluster -i thresh_zstat1 -c stats/cope1 -t 2.3 -d 0.388188 --volume=45089 --othresh=thresh_zstat1_2 -o cluster_mask_zstat1 --connectivity=26 --olmax=lmax_zstat1.txt --scalarname=Z -p 0.05 <b>--voxuncthresh</b>

Cluster Index	Voxels	P	-log10(P)	Z-MAX	Z-MAX X (vox)	Z-MAX Y (vox)	Z-MAX Z (vox)	Z-COG X (vox)	Z-COG Y (vox)	Z-COG Z (vox)	COPE-MAX	COPE-MAX X (vox)	COPE-MAX Y (vox)	COPE-MAX Z (vox)	COPE-MEAN
120	1217	3.89e-09	8.41	5.77	41	16	13	34	14.7	14	1.37e+03	40	16	13	269
119	518	2.73e-08	7.56	5.44	41	41	25	43.5	40.4	24	864	47	38	21	337
118	114	2.12e-06	5.67	4.6	47	18	18	47.1	19.2	19.7	474	48	17	19	241
117	79	2.22e-06	5.65	4.59	33	38	31	32.3	39.2	30.9	370	32	38	31	228
116	59	3.08e-06	5.51	4.52	17	40	18	17	39.5	19.2	675	15	38	19	334
115	55	9.72e-06	5.01	4.27	25	17	4	26.6	15.5	4.68	299	25	17	4	200
114	51	9.08e-05	4.04	3.74	38	45	13	38	42.1	13.6	452	37	42	15	245
[...]
</pre>

##### Cluster-wise corrected (same as before)
```
${FSLDIRUP}/src/cluster/cluster -i thresh_zstat1 -c stats/cope1 -t 2.3 -d 0.388188 --volume=45089 --othresh=thresh_zstat1_2 -o cluster_mask_zstat1 --connectivity=26  --olmax=lmax_zstat1.txt --scalarname=Z -p 0.05

Cluster Index   Voxels  P   -log10(P)   Z-MAX   Z-MAX X (vox)   Z-MAX Y (vox)   Z-MAX Z (vox)   Z-COG X (vox)   Z-COG Y (vox)   Z-COG Z (vox)   COPE-MAX    COPE-MAX X (vox)    COPE-MAX Y (vox)    COPE-MAX Z (vox)    COPE-MEAN
4   1217    3.61e-24    23.4    5.77    41  16  13  34  14.7    14  1.37e+03    40  16  13  269
3   518 4.57e-13    12.3    5.44    41  41  25  43.5    40.4    24  864 47  38  21  337
2   114 0.000717    3.14    4.6 47  18  18  47.1    19.2    19.7    474 48  17  19  241
1   79  0.00995 2   4.59    33  38  31  32.3    39.2    30.9    370 32  38  31  228
```
