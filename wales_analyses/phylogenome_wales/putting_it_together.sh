# BLAST them
# blastp -db CEGMA_reference_data/248.prots.fa -query SNAP-predicted/SNAP.export.aa.filtered.cutoff.TAIR.KOG.fa -max_hsps 1 -max_target_seqs 1 -outfmt "6 sacc qacc length evalue pident gaps" -num_threads 8|grep At > SNAP-predicted/SNAP.export.aa.filtered.cutoff.TAIR.KOG.CEGMA.hits.out
# blastx -db CEGMA_reference_data/248.prots.fa -query SNAP-predicted/SNAP.export.nt.filtered.cutoff.TAIR.KOG.fa -max_hsps 1 -max_target_seqs 1 -outfmt "6 sacc qacc length evalue pident gaps" -num_threads 8|grep At > SNAP-predicted/SNAP.export.nt.filtered.cutoff.TAIR.KOG.CEGMA.hits.out
# look at output, 
# first only KOG-identified genes (this is ALL of them)...
# those with matching KOG are good candidates...
# grep '|KOG' SNAP-predicted/SNAP.export.aa.filtered.cutoff.TAIR.KOG.CEGMA.hits.out 
# [no]:	At2g47020___KOG2726	518d7e70-0930-4986-a550-933a1bb61adf_Basecall_Alignment_templat-snap.1|AT4G27490|KOG1068	37	0.80	35.14	0
# [no]:	At2g43400___KOG2415	c7537a4e-a0ff-4a84-8b2a-871af8b68e3d_Basecall_Alignment_templat-snap.2|ATCG00120|KOG1353	65	4.7	23.08	0
# [YES]:	At5g24400___KOG3147	35d2492f-3dd8-4fbf-beaf-200dda979327_Basecall_Alignment_templat-snap.13|AT5G24400|KOG3147	69	0.005	24.64	0
# [no]:	At1g61010___KOG1137	52a4c646-b81d-456f-8dcf-f9b296fa07d6_Basecall_Alignment_templat-snap.6|AT5G16210|KOG0211	55	1.8	34.55	10
# [no]:	At1g13060___KOG0175	3427a1fa-af7c-41e4-b8ca-2a82c6cd066f_Basecall_Alignment_templat-snap.1|AT3G57190|KOG2726	146	0.008	26.71	27
# [YES]:	At3g24200___KOG3855	3e71121b-55c5-4d35-9626-5e2211c64b5d_Basecall_Alignment_templat-snap.6|AT3G24200|KOG3855	103	0.014	27.18	9
# [no]:	At3g10330___KOG1597	0db01e85-0099-4a7d-ba2d-2292f30a0808_Basecall_Alignment_templat-snap.3|AT2G02480|KOG0989	63	6.4	28.57	3
# [YES]:	At5g05470___KOG2916	04342120-a40f-4938-97fb-48a7fb6e18d5_Basecall_Alignment_templat-snap.1|AT5G05470|KOG2916	155	3e-06	30.97	39
# [YES]:	At1g77440___KOG0180	a05d3fe5-f157-45a8-b91a-f26a242af73a_Basecall_Alignment_templat-snap.1|AT1G21720|KOG0180	131	2e-15	41.22	21
# [YES]:	At3g58180___KOG0567	ae46393d-b55a-4ad4-97e1-4707194ad934_Basecall_Alignment_templat-snap.2|AT3G58180|KOG0567	126	1e-08	36.51	11
# [YES]:	At5g65720___KOG1549	ec1ac720-4b7b-4937-b817-5aa93b248603_Basecall_Alignment_templat-snap.1|AT5G65720|KOG1549	124	0.093	22.58	26
# [YES]:	At4g21110___KOG3404	3540415c-15f4-4e15-98d8-482949206861_Basecall_Alignment_templat-snap.4|AT4G21110|KOG3404	101	0.16	27.72	12
# [no]:	At3g62940___KOG2606	31cba0b4-d32b-4c6f-a05a-aef52bb73fd4_Basecall_Alignment_templat-snap.1|AT5G11110|KOG0853	85	2.4	23.53	8
# [no]:	At5g26675___KOG2519	08445df9-f81d-425e-8db7-432cbc09b632_Basecall_Alignment_templat-snap.1|AT1G69460|KOG1691	34	0.77	38.24	4
# [no]:	At1g44900___KOG0477	8b3b7350-53df-45d1-8f18-457bfcac8fee_Basecall_Alignment_templat-snap.2|ATCG00900|KOG3291	34	4.4	44.12	2
# [YES]:	At1g79210___KOG0181	2ed29331-4a23-4e67-9f98-7960b2d15d52_Basecall_Alignment_templat-snap.1|AT1G79210|KOG0181	134	3e-11	37.31	27
# [YES]:	At4g09320___KOG0888	1cf9abcb-d93c-46cd-b269-8eb83232daed_Basecall_Alignment_templat-snap.1|AT4G09320|KOG0888	151	4e-17	40.40	36
# [YES]:	At1g63780___KOG2781	c2bdc675-784e-48ba-8269-921bd9d5af76_Basecall_Alignment_templat-snap.3|AT1G63780|KOG2781	182	2e-19	37.36	67
# [YES]:	At2g25610___KOG0233	d4d7dbd6-2f3c-4085-9417-8ff1e909a423_Basecall_Alignment_templat-snap.1|AT4G32530|KOG0233	87	0.23	33.33	17
# [no]:	At4g31480___KOG1058	529535b0-d47c-4e09-983a-ea186f684268_Basecall_Alignment_templat-snap.1|AT1G26940|KOG0880	27	1.9	33.33	0

# etc..

# then those non-KOG-ID'd reads _may_ be good / usable if cutoff is low enough e.g.
# [definitely maybe]:		AtMi099___KOG1353	48c98bb5-8cf0-463e-b087-bbc74af1d938_Basecall_Alignment_templat-snap.2|AT2G07698|no_kog_value	120	2e-07	35.00	5
# [maybe]:				At5g20890___KOG0363	ee1d1fed-ca32-4912-b546-f865074ed1e1_Basecall_Alignment_templat-snap.1|AT5G40520|no_kog_value	140	0.004	23.57	23
# [no]:					At4g36480___KOG1358	121041c6-0242-4302-95f4-23afe6ea3ceb_Basecall_Alignment_templat-snap.1|ATCG01110|no_kog_value	58	0.28	31.03	6
# [no]:					At5g49650___KOG2531	f4ad2df3-9234-47c2-b20f-1dac5bcd7697_Basecall_Alignment_templat-snap.3|ATCG01000|no_kog_value	43	3.7	32.56	4
# [no]:					At5g26675___KOG2519	1e1a4c10-4839-4fc2-be6d-8ca60ecb47d7_Basecall_Alignment_templat-snap.2|AT1G10390|no_kog_value	29	2.2	41.38	5
# [no]:					At2g21390___KOG0292	6cb5d119-a8ad-4c67-a4dc-b6e9ff75ba1d_Basecall_Alignment_templat-snap.2|AT1G21630|no_kog_value	29	0.94	41.38	0
# [no]:					At1g64550___KOG0062	82216d1a-e320-4d9c-9d76-965aa1e237dc_Basecall_Alignment_templat-snap.3|AT1G51405|no_kog_value	30	1.1	40.00	0
# [no way]:				At4g26870___KOG0556	a6ce2b5e-4127-406c-9ca2-802c14d060e1_Basecall_Alignment_templat-snap.1|AT4G38010|no_kog_value	36	1.8	33.33	0
# [no]:					At3g54670___KOG0018	fd647e51-f95b-4d28-a9a2-f7aca7dd898c_Basecall_Alignment_templat-snap.3|AT5G19350|no_kog_value	40	2.1	25.00	0
# [maybe maybe]:			At1g56075___KOG0469	2e1a7f1e-5d1f-4d6e-a203-eab81030d948_Basecall_Alignment_templat-snap.1|AT2G23790|no_kog_value	105	0.24	24.76	12

# etc..


# so make trees now with:
# [YES]:	At1g77440___KOG0180	a05d3fe5-f157-45a8-b91a-f26a242af73a_Basecall_Alignment_templat-snap.1|AT1G21720|KOG0180	131	2e-15	41.22	21
# [YES]:	At1g79210___KOG0181	2ed29331-4a23-4e67-9f98-7960b2d15d52_Basecall_Alignment_templat-snap.1|AT1G79210|KOG0181	134	3e-11	37.31	27
# [YES]:	At2g25610___KOG0233	d4d7dbd6-2f3c-4085-9417-8ff1e909a423_Basecall_Alignment_templat-snap.1|AT4G32530|KOG0233	87	0.23	33.33	17
# [YES]:	At3g58180___KOG0567	ae46393d-b55a-4ad4-97e1-4707194ad934_Basecall_Alignment_templat-snap.2|AT3G58180|KOG0567	126	1e-08	36.51	11
# [YES]:	At4g09320___KOG0888	1cf9abcb-d93c-46cd-b269-8eb83232daed_Basecall_Alignment_templat-snap.1|AT4G09320|KOG0888	151	4e-17	40.40	36
# [YES]:	At5g65720___KOG1549	ec1ac720-4b7b-4937-b817-5aa93b248603_Basecall_Alignment_templat-snap.1|AT5G65720|KOG1549	124	0.093	22.58	26
# [YES]:	At1g63780___KOG2781	c2bdc675-784e-48ba-8269-921bd9d5af76_Basecall_Alignment_templat-snap.3|AT1G63780|KOG2781	182	2e-19	37.36	67
# [YES]:	At5g05470___KOG2916	04342120-a40f-4938-97fb-48a7fb6e18d5_Basecall_Alignment_templat-snap.1|AT5G05470|KOG2916	155	3e-06	30.97	39
# [YES]:	At5g24400___KOG3147	35d2492f-3dd8-4fbf-beaf-200dda979327_Basecall_Alignment_templat-snap.13|AT5G24400|KOG3147	69	0.005	24.64	0
# [YES]:	At4g21110___KOG3404	3540415c-15f4-4e15-98d8-482949206861_Basecall_Alignment_templat-snap.4|AT4G21110|KOG3404	101	0.16	27.72	12
# [YES]:	At3g24200___KOG3855	3e71121b-55c5-4d35-9626-5e2211c64b5d_Basecall_Alignment_templat-snap.6|AT3G24200|KOG3855	103	0.014	27.18	9


# run from ~/Downloads/phylogenome_wales/SNAP_KOG_great/
# additional gneomes: 
#	../C.reinhardtii.cegma.aa.fa
#	../O.sativa.cegma.aa.fa
#	../P.trichocarpa.cegma.aa.fa
#
# first run:
#	rename_CEGMA_taxa.py
# to rename the KOG0000.aln to KOG0000.aln.renamed with constant taxon names
#
# then run 
#	interleave_prediction_with_CEGMA_reference_data.py --fasta [../C.reinhardtii.cegma.aa.fa|../O.sativa.cegma.aa.fa|../P.trichocarpa.cegma.aa.fa] --output --output CEGMA_reference_data/248alignments/KOG0188.aln.renamed.added --name [arg] --cegma CEGMA_reference_data/248alignments/KOG0188.aln.renamed.added
# for each of the addiional genomes to interleave them
#
# then run e.g
#	python interleave_prediction_with_CEGMA_reference_data.py --fasta SNAP_KOG_maybe/SNAP.export.aa.filtered.cutoff.TAIR.KOG.fa.selected --output CEGMA_reference_data/248alignments/KOG0188.aln.renamed.added.TEST --name A_thal_AAsnap__ --cegma CEGMA_reference_data/248alignments/KOG0188.aln.renamed.added
# for each KOG in the list to add prediction to alignments.

raxml='/Applications/Phylogenetics/RAxML/RAxML-7.2.8-ALPHA/raxmlHPC'

#loops, yay
for f in *KOG.fa; 
	do 
# name the files
ALIGNMENT=$f.aligned.fa
TRIMMED=$f.aligned.trimal.0.9.phy

 muscle -in $f -out $ALIGNMENT -quiet 
 # convert to phylip for RAxML 
 #python /Users/joeparker/Downloads/phylogenome_wales/fasta_to_phylip.py --fasta $f.aligned.fa --output $f.aligned.phy; 
# convert to phylip and trim at the same time
trimal -in $ALIGNMENT -out $TRIMMED -phylip_m10 -gt 0.8
 # RAxML, dayhoff matrix, pthreads=8
# /Applications/Phylogenetics/RAxML/RAxML-7.2.8-ALPHA/raxmlHPC -s $f.aligned.phy -n $f -m PROTCATDAYHOFF; 
 /Applications/Phylogenetics/RAxML/RAxML-7.2.8-ALPHA/raxmlHPC-PTHREADS -T 8 -s $TRIMMED -m PROTCATDAYHOFF -n $f  
 mv RAxML_bestTree.$f $f.tre; 
 rm RAxML*; 
 done 
 
 # then cat the RAxmL outputs together
 # add NEXUS formatting
 # stick in treeannotator for a quick MJR consensue
 
 #### Update 2016 12 05
 # Have checked that the SNAP transcripts (export.tx) are in frame and they are.. 

# Including all the AA genes with length ≥ 100AAs:
# egrep -e '[0-9a-zA-Z\_]{1,}\t[0-9]{3,}(\t[0-9\.e\-]{1,}){3}' ~/Downloads/phylogenome_wales/SNAP-predicted/SNAP.export.aa.filtered.cutoff.TAIR.KOG.CEGMA.hits.out > ~/Downloads/phylogenome_wales/SNAP.export.aa.filtered.cutoff.TAIR.KOG.CEGMA.hits.100aa.out
# gives 125 hits
# of which 10 have KOGs already
# of which 9 have matching KOG tags (the same ones listed above)