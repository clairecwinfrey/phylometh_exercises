#!/bin/bash

#RAxML tutorial
#March 14, 2020


################################################################################
#                            GETTING STARTED
################################################################################
#BINGAMMA option tells RAXML that we are using binary data amd we want to use the
#GAMMA model of rate heterogeneity
#Simple ML search on binary data:
#-n is arbitrary file name appendix
#T1 makes it so that all output files will be names RAxML_FileTypeName.T1 (each run /
#needs unique name)
raxmlHPC -m BINGAMMA -p 12345 -s binary.phy -n T1 \

#Can also use CAT, a lower memory option, for datasets with at least 50-100 taxa
raxmlHPC-SSE3 -m BINCAT -p 12345 -s binary.phy -n T2 \

#RAxML uses randomized stepwise parsimony trees each time, so to get same result each time
#set fixed random number seed.
raxmlHPC -m BINGAMMA -p 12345 -s binary.phy -n T3 \

#or pass in a starting tree:
raxmlHPC -m BINGAMMA -t startingTree.txt -s binary.phy -n T4 \

#or look for best tree (e.g. 20 ML searches on 20 random trees below)
raxmlHPC -m BINGAMMA -p 12345 -s binary.phy -# 20 -n T5 \

############DNA########################
raxmlHPC -m GTRGAMMA -p 12345 -s dna.phy -# 20 -n T6 \
##odd that the best tree so clearly messes up relationships b/w these taxa
##final GAMMA-based score of best tree was -376.281437

###########PROTEINS######################
#1. using WAG model (has pre-defined base frequency):
raxmlHPC -m PROTGAMMAWAG -p 12345 -s protein.phy -# 20 -n T7 \
#again, phylo relationships between these taxa are very wrong...
#final GAMMA-based score of best tree -434.406886

#2. using empirical base frequencies drawn from alignment (i.e. basically counting
#the frequency of the occurrence of the amino acids), e.g. JTTF:
raxmlHPC -m PROTGAMMAJTTF -p 12345 -s protein.phy -n T8 \
#Final GAMMA-based Score of best tree -407.554591

#3. can also use CAT approximation of rate heterogenity again:
raxmlHPC -m PROTCATJTTF -p 12345 -s protein.phy -n T9 \
#Final GAMMA-based Score of best tree -407.554603

###TIPS for choosing protein model. Generate reasonable reference tree, compute likelihood scores,
#under all available models. Then use model that results in best score under GAMMA.
#To do this automatically in RAxML:
raxmlHPC -p 12345 -m PROTGAMMAAUTO -s protein.phy -n AUTO \
#Best-scoring ML tree written to: /Users/clairewinfrey/Desktop/Spring_2020/Phylometh_2020/RAxML/RAxML_bestTree.AUTO
#Final GAMMA-based Score of best tree -401.563151

#GTR models with protein (best to only do w/ large datasets b/c 189 rates (b/w all amino acids) in matrix):
raxmlHPC -p 12345 -m PROTGAMMAGTR -s protein.phy -n GTR \
#Final GAMMA-based Score of best tree -362.656788

#Conduct ML search on multi-state morpho. datasets
raxmlHPC -p 12345 -m MULTIGAMMA -s  multiState.phy -n T10 \
#Final GAMMA-based Score of best tree -226.253232

#-K option lets you choose different models available for multistate characters
#e.g. ORDERED, MK, OR GTR (but default is GTR)
raxmlHPC -p 12345 -m MULTIGAMMA -s  multiState.phy -K MK -n T11 \
#Final GAMMA-based Score of best tree -323.415298
raxmlHPC -p 12345 -m MULTIGAMMA -s  multiState.phy -K ORDERED -n T12 \
#Final GAMMA-based Score of best tree -345.774763

################################################################################
#                            BOOTSTRAPPING
################################################################################
#Generate 20 ML trees on distinct starting trees and print best-scoring ML tree.
raxmlHPC -m GTRGAMMA -p 12345 -# 20 -s dna.phy -n T13 \
#Final GAMMA-based Score of best tree -376.281437

#Get support values for the tree by conducting a bootstrap search.
#-b sets random number seed and -# is number of bootstrap replicates.
raxmlHPC -m GTRGAMMA -p 12345 -b 12345 -# 100 -s dna.phy -n T14 \

#RAxML_bootstrap.T14 was name of last printed file, and we can use it to draw
#bipartitions on the best ML tree. Specifically, RAxML_bipartitions.T15 are support
# values assigned to node and RAxML_bipartitionsBranchLabels.T15 are support values
#assigned to branches of the tree)
###BECAUSE THIS TREE IS UNROOTED, THE CORRECT REPRESENTATION IS THE ONE W/ support
#VALUES ASSIGNED TO BRANCHES AND NOT NODES OF THE TREE!
raxmlHPC -m GTRCAT -p 12345 -f b -t RAxML_bestTree.T13 -z RAxML_bootstrap.T14 -n T15 \

####using bootstrapping methods to build consensus- three ways:
#1. strict consensus
raxmlHPC -m GTRCAT -J STRICT -z RAxML_bootstrap.T14 -n T16 \

#2. majority rule
raxmlHPC -m GTRCAT -J MR -z RAxML_bootstrap.T14 -n T17 \

#3. extended majority rule:
raxmlHPC -m GTRCAT -J MRE -z RAxML_bootstrap.T14 -n T18 \

################################################################################
#                        RAPID BOOTSTRAPPING
################################################################################
#Rapid bootstrapping is one order of magnitude faster than standard algorithm
#use -x instead of -b for random number and can choose different models
raxmlHPC -m GTRGAMMA -p 12345 -x 12345 -# 100 -s dna.phy -n T19 \

#complete analysis in one step: 100 rapid Bootstrap analyses, 20 ML searches,
#and return best ML tree (with support values):
raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s dna.phy -n T20 \
#Final ML Optimization Likelihood: -376.379064

################################################################################
#                         PARTITIONED ANALYSES
################################################################################
#To conduct partitioned analyses, pass info about partitions to RAxML using
#plain text file that is passed via the -q parameter. In the example below, the
#file simpleDNApartition partitions the alignment into two regions (DNA, p1=1-30),
#DNA, p2= 31-60. Here, different parameters, base frequencies, and evol rates in
#GTR matrix are estimated independently for every partition.

raxmlHPC -m GTRGAMMA -p 12345 -q simpleDNApartition.txt -s dna.phy -n T21 \
##Final GAMMA-based Score of best tree -375.308100

#-M makes it so that estimates a separate set of branch lengths for every partition:
raxmlHPC -M -m GTRGAMMA -p 12345 -q simpleDNApartition.txt -s dna.phy -n T22 \
#Final GAMMA-based Score of best tree -359.885537

#more elaborate partitioning by 1st, 2nd, and 3rd codon position:
#Partition file:
#DNA, p1=1-60\3,2-60\3
#DNA, p2=3-60\3
raxmlHPC -m GTRGAMMA -p 12345 -q dna12_3.partition.txt -s dna.phy -n T23 \
#Final GAMMA-based Score of best tree -357.233681

#For datasets that contain DNA and proteins:
#partition file: DNA, p1 = 1- 50; WAG, p2 = 51-110
raxmlHPC -m GTRGAMMA -p 12345 -q dna_protein_partitions.txt -s dna_protein.phy -n T24 \

#This is equal to T24 above
raxmlHPC -m PROTGAMMAWAG -p 12345 -q dna_protein_partitions.txt -s dna_protein.phy -n T25 \

#More on partition files:BIN is binary (e.g. BIN, p1=1-50); MULTI is multi-partitioned
#-K  allows you to specify your sub. model for multi-state regions

################################################################################
#                        SECONDARY STRUCTURE MODELS
################################################################################

#To specify secondary structure models for RNA alignment: 1) read in plain RNA alignment
#2) use -S to pass in an additional text file that specifies which RNA alignment estimates
#need to be grouped together.

#Example of plain text file using standard bracket notation. If, for example, DNA
#test alignment has 60 sites, secondary structure file needs to contain a string of 60
#characters:

#..................((.......))............(.........)........
#periods mean a normal RNA site; brackets (parenthesis?) mean stems; <>, {}, and [],
#can be used to indicate psuedo knots.

#use -A option to specify different number-state models (e.g. S6A, S6C, etc.)

raxmlHPC -m GTRGAMMA -p 12345 -S secondaryStructure.txt -s dna.phy -n T26 \
#Final GAMMA-based Score of best tree -365.275413

################################################################################
#                       USING THE PTHREADS VERSION
################################################################################
#USING THE PTHREADS VERSION
#NEVER RUN PTHREADS WITH MORE THREADS THAN YOU HAVE CORES!

#To do so: use correct executable *(e.g. (raxmlHPC-PTHREADS or raxmlHPC-PTHREADS-SSE3)
#then use -T to specify number of threads to use.

#E.g. my laptop has four cores:
raxmlHPC-PTHREADS -T 4 -p 12345 -m PROTGAMMAWAG -s protein.phy -n T27 \

################################################################################
#                             CONSTRAINT TREES
################################################################################
#To compute liklihood of a tree, the tree must be resolved.

#You can use constraint trees to test hypotheses of monophyly by comparing trees
#generated by different constraints using liklihood-based significance tests.

#If a constraint tree is specified that does not place all taxa, these taxa can
#potentially be inserted anywhere in the tree. You have to specify, e.g.
#"((A1,A2,A3,A4,A5),(B1,B2),(X1,X2,...,X17);" if you DON'T want say X1, X2, and X17
#to go within clades A and/or B.

#RAxML has: 1) binary backbone constraint trees and  2) multifurcating constraint trees

raxmlHPC-SSE3 -P12345 -r backboneConstraint.txt -m GTRGAMMA -s dna.phy -n T28 -p 5 \


#For multifurcating constraint trees, monophyletic structure is maintained, but
#taxa within mono. clade may be moved around.

raxmlHPC-SSE3 -p 12345 -g multiConstraint.txt -m GTRGAMMA -s dna.phy -n T29 \
#Final GAMMA-based Score of best tree -381.853476

########################## CONSTRAINING SISTER TAXA ############################

#example of how to constrain potential sister taxa to monophyly by putting multi-furcating
#backboen constraint into a text file (here, monophylyConstraint.txt)

raxmlHPC-SSE3 -p 12345 -g monophylyConstraint.txt -m GTRGAMMA -s dna.phy -n T29monophyly \

#can force bad monophyly too, but ML score will be bad! Can then use liklihood tests to
#determine whether two hypotheses of monophyly yield significantly different LnL scores.

################################################################################
#                               OUTGROUPS
################################################################################
#outgroups allow you to root the tree at the branch which leads to that outgroup.

#specify a single outgroup (i.e mouse):
raxmlHPC-SSE3 -p 12345 -o Mouse -m GTRGAMMA -s dna.phy -n T30 \

#With mouse and rats as outgroups:
raxmlHPC-SSE3 -p 12345 -o Mouse,Rat -m GTRGAMMA -s dna.phy -n T31 \
