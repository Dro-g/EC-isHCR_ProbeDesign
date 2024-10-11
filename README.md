# Pipeline of split-initiator probe design

## 0. environment
We have verified that our R files work in the following environments.
R version 4.3.3 (2024-02-29)

## I. Workflow of the pipeline
1. Find candidate regions for generating probes for the target mRNA (First_ver2.R).
2. Identify similarities in all RNAs except the target (Blastn).
3. Add the initiator sequence (Third_ver1.R).
 
## II. Glossary
 
 <p align="center">
<img src="ImagesREADME/fig1.png" width="66%">  
</p>
<p align="center">
 Figure 1
</p>

 
- The ***probe region*** comprises a pair of probe binding sites and a ***spacer***. Considering that the probe binding site is 25 nucleotide (nt) and the spacer is 2 nt, the probe region is 52 nt.
- The probe for the 5′ mRNA region is Probe1 (***P1***), and the probe for the 3′ mRNA region is Probe2 (***P2***).
- Each probe contains part of the ***initiator*** sequence, and the neighboring probes complete the initiator sequence.
 
## III. Criteria for the mRNA binding sites of probes
-  The GC content should be with 40–60% (45–55% is recommended)<sup>1</sup>.
- In the probe, the “AA” bases connect the RNA binding site and the initiator as a linker. Because binding of the linker to target mRNA can alter the staining efficiency, the probe region at which the linker binds to mRNA (probe region including “T” in the spacer) is excluded as a default.
- Designing with CDS of mRNA is recommended<sup>2</sup>.
 
## IV. Overview of the First_ver2.R program
 
<p align="center">
<img src="ImagesREADME/fig2.png" width="66%">  
</p>
<p align="center">
 Figure 2
</p>

1. As shown in Figure 2, the 1st-52nd sequence of the target mRNA (green letters in Figure 2) is checked to meet the requirements.
2. If the region meets the criteria, it is selected as a candidate Probe Region, and the sequence is shifted by 52 bases towards the 3′ end and examined. If it does not meet the requirements, the sequence is shifted by 1 base and examined.
3. Step 2 is repeated until a satisfactory number of probe regions is obtained.
 
## V. Procedure for the First_ver2.R program
1. The DNA sequence of the target mRNA is downloaded and saved as a text file (.txt).
2. “First_ver2.R” is opened, and the file name is written with an absolute path in “LoadFileName”.
3. The number of necessary probe regions is written in “CandidateNum”. The minimum number of probe sets in the previous study is five and increasing the number of probe sets improves the signal-to-noise ratio<sup>3</sup>. Because the program searches from the 5′ end, the probe regions are biased to the 5′. Therefore, a larger number of probes should be set.
4. The entire program is run.
5. The candidate sequences of the probe regions are output in fasta format. Information about the candidate probe regions is also output in text format. Details of the information file are below.

***Probe_Region_Sense*** is the mRNA sequence of the probe region. Sense indicates the same sequence as that of the target mRNA. 
***PRS GC Content*** is the GC content of the Probe Region:
***P1_GC Content*** is the GC content of P1,
***P2_GC Content*** is the GC content of P2,
***StartBp_num*** indicates the location of the first base of the probe region, and
***EndBp_num*** indicates the location of the end base.

If an insufficient number of probe regions are obtained, the following steps are taken:
- Relieve GC content condition. Set `gc_min` to **40** and `gc_max` to **60 in lines 11–12** of First_ver2.R.
- Relieve linker (spacer) condition. Remove
` & result[[7]] != "T" & result[[8]] != "T"` 
**in line 66** of First_ver2.R.
- Search the probe region in the UTR.
 
## VI. About Blastn
Only an adjacent pair of probes can form the initiator sequence and induce the amplification of the hairpin DNA<sup>2</sup>; therefore, we performed Blastn on the candidate probe regions rather than on individual probes. Because the probe regions are sense sequences, we searched for +/+ fields in the Blastn results. The fasta file of the probe region output in procedure 5, section V, can be used for the Blastn.
 
## VII. Procedure of the Third_ver1.R program
This program generates P1 and P2 by reverse-complementation and splitting of the probe region, and the sections are conjugated parts of the initiator sequence.
 
1. A table of the necessary probe regions is prepared. The rows containing unnecessary probe regions are removed from the information file (.txt) output in procedure 5, section V.
2. Third_ver1.R is opened, and the R files for S45 or A161 hairpin DNA are available.
3. The table file name (.txt) is written with an absolute path in “RegionFile”.
4. The entire program is run.
5. The Probe table is output as a CSV file.
 
## VIII. References
1. https://nepagene.jp/wp-content/uploads/ISHpalette_probe-design_v1-j.pdf
2. https://sites.google.com/view/in-situ-shhcr/faq
3. Choi HMT, Schwarzkopf M, Fornace ME, Acharya A, Artavanis G, Stegmaier J, Cunha A, Pierce NA. Third-generation _in situ_ hybridization chain reaction: multiplexed, quantitative, sensitive, versatile, robust. Development 2018; 145:dev165753.
 
 
 
