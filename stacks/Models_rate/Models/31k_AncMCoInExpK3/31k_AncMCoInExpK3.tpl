//Parameters for the coalescence simulation program : simcoal.exe
3 samples to simulate :
//Population effective sizes (number of genes)
NPOP0$
NPOP1$
NPOP2$
//Samples sizes and samples age 
8
12
19
//Growth rates: negative growth implies population expansion
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 0 0
0 0 0
0 0 0
//Migration matrix 1
0 0 MIG120$
0 0 0
MIG012$ 0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
5 historical event
TRESIZE$ 0 0 0 RESIZE_0$ 0 0
TRESIZE$ 1 1 0 RESIZE_1$ 0 0
TRESIZE$ 2 2 0 RESIZE_2$ 0 0
TDIV1$ 1 2 1 RESIZE_3$ 0 1
172 2 0 1 RESIZE_4$ 0 1
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 MUTRATE$