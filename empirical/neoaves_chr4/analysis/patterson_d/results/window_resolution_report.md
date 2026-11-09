# Patterson's D window-resolution analysis

Input: `/Users/ytabatabaee/Desktop/msrc-paper/empirical/neoaves_chr4/analysis/patterson_d/results/patterson_d_locus_results.tsv`; locus-level counts were reused without reprocessing nucleotide alignments.

Low-information flag: fewer than 50 informative sites or fewer than 5 loci. Empty windows are retained with blank Patterson_D values.

Window tracks are visualization summaries. They are distinct from the existing 500-kb, 1-Mb, and 2-Mb genomic jackknife blocks, which were not changed.

## Resolution summary

resolution	design	filter	total_windows	defined_D	empty_windows	low_information_windows	low_information_fraction	min_informative_sites	median_informative_sites	max_informative_sites
100kb	Columbea_alt_outgroup	all	902	862	40	254	0.281596	0	98	239
100kb	Columbea_alt_outgroup	transversions_only	902	859	43	901	0.998891	0	15	51
100kb	Columbea_primary	all	902	862	40	246	0.272727	0	102	242
100kb	Columbea_primary	transversions_only	902	858	44	901	0.998891	0	15	53
100kb	N61_alt_outgroup	all	902	862	40	260	0.288248	0	113	240
100kb	N61_alt_outgroup	transversions_only	902	857	45	901	0.998891	0	17	53
100kb	N61_primary	all	902	862	40	253	0.280488	0	116	271
100kb	N61_primary	transversions_only	902	859	43	900	0.997783	0	18	53
100kb	N62_alt_outgroup	all	902	859	43	282	0.312639	0	107	266
100kb	N62_alt_outgroup	transversions_only	902	855	47	900	0.997783	0	17	56
100kb	N62_primary	all	902	859	43	277	0.307095	0	110	253
100kb	N62_primary	transversions_only	902	852	50	900	0.997783	0	17	58
250kb	Columbea_alt_outgroup	all	361	361	0	14	0.038781	8	243	405
250kb	Columbea_alt_outgroup	transversions_only	361	361	0	284	0.786704	1	37	85
250kb	Columbea_primary	all	361	361	0	14	0.038781	8	243	420
250kb	Columbea_primary	transversions_only	361	361	0	281	0.778393	1	38	86
250kb	N61_alt_outgroup	all	361	360	1	16	0.044321	0	268	493
250kb	N61_alt_outgroup	transversions_only	361	359	2	225	0.623269	0	43	101
250kb	N61_primary	all	361	360	1	16	0.044321	0	278	525
250kb	N61_primary	transversions_only	361	359	2	219	0.606648	0	44	103
250kb	N62_alt_outgroup	all	361	360	1	16	0.044321	0	261	497
250kb	N62_alt_outgroup	transversions_only	361	360	1	242	0.670360	0	42	89
250kb	N62_primary	all	361	360	1	16	0.044321	0	267	490
250kb	N62_primary	transversions_only	361	360	1	233	0.645429	0	43	86
500kb	Columbea_alt_outgroup	all	181	181	0	0	0.000000	124	474	720
500kb	Columbea_alt_outgroup	transversions_only	181	181	0	31	0.171271	13	72	151
500kb	Columbea_primary	all	181	181	0	0	0.000000	121	486	735
500kb	Columbea_primary	transversions_only	181	181	0	27	0.149171	13	76	149
500kb	N61_alt_outgroup	all	181	181	0	0	0.000000	120	539	853
500kb	N61_alt_outgroup	transversions_only	181	181	0	21	0.116022	14	84	168
500kb	N61_primary	all	181	181	0	0	0.000000	117	548	883
500kb	N61_primary	transversions_only	181	181	0	18	0.099448	17	88	185
500kb	N62_alt_outgroup	all	181	181	0	0	0.000000	105	521	908
500kb	N62_alt_outgroup	transversions_only	181	181	0	28	0.154696	12	83	171
500kb	N62_primary	all	181	181	0	0	0.000000	115	531	911
500kb	N62_primary	transversions_only	181	181	0	24	0.132597	12	84	170
1Mb	Columbea_alt_outgroup	all	91	91	0	0	0.000000	261	936	1356
1Mb	Columbea_alt_outgroup	transversions_only	91	91	0	1	0.010989	30	140	235
1Mb	Columbea_primary	all	91	91	0	0	0.000000	279	981	1357
1Mb	Columbea_primary	transversions_only	91	91	0	1	0.010989	34	146	243
1Mb	N61_alt_outgroup	all	91	91	0	0	0.000000	249	1044	1674
1Mb	N61_alt_outgroup	transversions_only	91	91	0	1	0.010989	28	167	317
1Mb	N61_primary	all	91	91	0	0	0.000000	261	1100	1710
1Mb	N61_primary	transversions_only	91	91	0	1	0.010989	27	177	341
1Mb	N62_alt_outgroup	all	91	91	0	0	0.000000	268	992	1444
1Mb	N62_alt_outgroup	transversions_only	91	91	0	1	0.010989	43	161	273
1Mb	N62_primary	all	91	91	0	0	0.000000	271	1029	1482
1Mb	N62_primary	transversions_only	91	91	0	1	0.010989	38	165	268

## N61 regional localization

The complete interval summaries, including pooled counts and defined-window ranges, are in `window_region_summary.tsv`. For N61 primary/all-sites, pooled D remains negative in all three intervals: approximately -0.258, -0.330, and -0.251 at 100 kb; -0.258, -0.331, and -0.249 at 250 kb; -0.253, -0.331, and -0.249 at 500 kb; and -0.253, -0.247, and -0.249 at 1 Mb for PNAS_1, PNAS_2, and PNAS_3 respectively. PNAS_2 has fewer windows and lower information, especially at fine resolution. This supports a distributed interval-level shift rather than a single isolated bin, but no individual boundary is inferred from visualization alone.

## Resolution recommendation

Use 250 kb as the main visualization compromise: 100-kb all-site tracks retain useful detail but have roughly 27--31% low-information windows across primary designs, while 100-kb transversion-only tracks are almost entirely low-information. At 250 kb, all-site tracks have about 96% defined windows and far fewer low-information bins, while retaining more localization than 500 kb or 1 Mb. Keep 100 kb as an exploratory view and rely on 500 kb or larger for sparse transversion-only interpretation.

## Preservation

All rows in `pooled_result_preservation.tsv` must have status PASS. These checks compare pooled ABBA/BABA counts against the original summary table.

## Interpretation

Finer windows increase spatial detail but produce more low-information bins. Local extreme D values in low-information windows should not be interpreted as strong evidence. The N61 concentration assessment should use the tracks together with informative-site counts and the existing region-level pooled estimates; no boundary is inferred from visualization alone.
