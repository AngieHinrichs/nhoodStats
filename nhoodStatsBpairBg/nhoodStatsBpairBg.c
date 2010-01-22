/* nhoodStatsBpairBg - Compute avg neighborhood corr for genes across window 
 * sizes and p-val vs rand permuts. */
/* Copyright (c) 2010 Danielle G. Lemay and Angie S. Hinrichs. All rights reserved. */
#include "common.h"
#include "hash.h"
#include "linefile.h"
#include "obscure.h"
#include "options.h"
#include "portable.h"
#include "sqlNum.h"

static char const rcsid[] = "$Id: Exp $";

// Globals for command line options:
int permCount = 5;

void usage()
/* Explain usage and exit. */
{
errAbort(
"nhoodStatsBpairBg - Compute average neighborhood correlation (ANC) for\n"
"                    clusters of genes across window sizes, and p-values of\n"
"                    those ANCs vs ANCs calculated using shuffled pairwise\n"
"                    correlation values.\n"
"usage:\n"
"   nhoodStatsBpairBg corr geneStarts chr# minBlkSize maxBlkSize resultDir\n"
"options:\n"
"   -perms=N    Number of random permutation runs for calculating p-vals\n"
"               (default %d)\n",
permCount);
}

// Note: if this becomes >= 64k-1 (because we use max+1 as error/extra), 
// can't use unsigned short for corr indexes anymore (use int):
#ifndef MAX_GENES_IN_GENOME
#define MAX_GENES_IN_GENOME ((64 * 1024) - 2)
#endif

#ifndef MAX_GENES_IN_CLUSTER
#define MAX_GENES_IN_CLUSTER 1024
#endif

static struct optionSpec options[] = {
    {"perms", OPTION_INT},
    {NULL, 0},
};

// TODO: make geneSets a linked list + ANC.  Then clustersANCs doesn't statically alloc.

/****************************** Data structures ******************************/
/* Two structs are defined here in order to compute and store clusters and 
 * their ANCs:
 *
 * struct geneSet: This represents a collection of genes as start coordinates 
 *                 and indices into the pairwise correlation values 2-D array, 
 *                 and the array that maps IDs to gene names.
 *                 GeneSets are used in two different ways:
 *                 1. the set of all genes on a chrom / all genes in genome
 *                 2. the set of genes in a cluster
 *
 * struct clustersANCs: This represents a collection of gene clusters and 
 *                      their ANCs, for a particular window size.
 *
 * For each window size, clusters and their ANCs are computed once and stored
 * in a struct clustersANCs.
 *
 * In a second pass to estimate the p-value of ANCs, we iterate through the
 * window sizes and clusters, recomputing ANCs using shuffled IDs (hence 
 * shuffled pairwise correlation values) and storing that randomized 
 * distribution.
 *
 * Since each cluster's ANC must be compared to the entire distribution of
 * randomized ANCs, and there are many of them (number of clusters in all
 * windows * number of permutations), a fast comparison method is important
 * for performance.  Therefore the array of randomized ANCs is sorted in
 * ascending order using qsort, and binary search is used to identify each
 * ANC's place in that order.
 */

struct geneSet
/* A set of genes with indexes and start coords. */
    {
    int n;                 // # of genes in set
    int maxN;              // allocated size of arrays
    unsigned short *cixs;  // each gene's index into geneIDs and corrValues
    unsigned int *starts;  // each gene's start coord.
    };

struct clustersANCs
/* Cluster geneSets and neighborhood ANC distribution for a single window 
 * size. */
    {
    int winSize;                // window size for clustering
    int n;                      // number of clusters found
    int maxN;                   // allocated size of arrays
    float *ancs;                // each cluster's ANC
    struct geneSet **clusters;  // each cluster's gene starts and indexes
    };

struct geneSet *geneSetNew(int maxGenes)
/* Alloc & return empty geneSet that can hold up to maxGenes genes. */
{
struct geneSet *gs;
AllocVar(gs);
gs->n = 0;
gs->maxN = maxGenes;
AllocArray(gs->cixs, maxGenes);
AllocArray(gs->starts, maxGenes);
return gs;
}

void geneSetAdd(struct geneSet *gs, unsigned short cix, unsigned int start)
/* Add a gene's index and start coord. */
{
if (gs->n >= gs->maxN)
    errAbort("geneSetAdd: add reallocation to geneSet code or increase max from %d", gs->maxN);
gs->cixs[gs->n] = cix;
gs->starts[gs->n] = start;
gs->n++;
}


int readCorrFile(struct lineFile *corrLf, char ***retGeneIDs,
		 float ***retCorrValues, struct hash **retGeneIDToCix)
/* Read in gene IDs and 2-D table of floating point gene pair correlation values.
 * Make hash of gene IDs to correlation table indexes (cix). */
{
// First line of correllation file contains gene ID words:
int lineSize;
char *line;
if (! lineFileNext(corrLf, &line, &lineSize))
    lineFileAbort(corrLf, "Failed to read first line");
char *geneIDStr = cloneStringZ(line, lineSize);
stripChar(geneIDStr, '"');
char **geneIDs;
AllocArray(geneIDs, MAX_GENES_IN_GENOME+1);
int geneCount = chopByWhite(geneIDStr, geneIDs, MAX_GENES_IN_GENOME+1);
if (geneCount > MAX_GENES_IN_GENOME)
    errAbort("Recompile with larger MAX_GENES_IN_GENOME (%d)", MAX_GENES_IN_GENOME);
verbose(1, "Got %d gene IDs\n", geneCount);

// Subsequent lines contain geneID followed by geneCount floating point values. 
// Store in 2-D array (geneCount x geneCount):
struct hash *geneIDToCix = hashNew((int)(log(MAX_GENES_IN_GENOME) / log(2) + 1.5));
char *words[MAX_GENES_IN_GENOME+1];
int wordCount;
float **corrValues;
AllocArray(corrValues, geneCount);
int geneCix = 0;
while ((wordCount = lineFileChop(corrLf, words)) > 0)
    {
    if (wordCount != geneCount+1)
	lineFileAbort(corrLf, "Expected geneCount+1=%d words, got %d", geneCount+1, wordCount);
    stripChar(words[0], '"');
    hashAddInt(geneIDToCix, geneIDs[geneCix], geneCix);
    AllocArray(corrValues[geneCix], geneCount);
    if (! sameString(words[0], geneIDs[geneCix]))
	lineFileAbort(corrLf, "Row geneIDs are not in same order as header gene IDs "
		"(expected %s, got %s)", geneIDs[geneCix], words[0]);
    int j;
    for (j = 0;  j < geneCount;  j++)
	corrValues[geneCix][j] = sqlFloat(words[j+1]);
    geneCix++;
    }
if (geneCix != geneCount)
    lineFileAbort(corrLf,
	    "Expected geneCount=%d lines of float values after gene ID line, but got %d lines",
	     geneCount, geneCix);
verbose(1, "Read %dx%d values\n", geneCount, geneCount);
*retGeneIDs = geneIDs;
*retCorrValues = corrValues;
*retGeneIDToCix = geneIDToCix;
return geneCount;
}

void readGeneFile(struct lineFile *geneLf, char *chrNum, struct hash *geneIDToCix, int maxIDs,
		  struct geneSet **retGenesOnChrom, struct geneSet **retGenesInGenome)
/* Read file of [gene ID, chr# (a la Ensembl), start coordinate]; build two sets of genes,
 * one for genes on chrNum and one for all genes in the genome. */
{
struct geneSet *genesOnChrom = geneSetNew(maxIDs), *genesInGenome = geneSetNew(maxIDs);
int wordCount;
char *words[4];
while ((wordCount = lineFileChop(geneLf, words)) > 0)
    {
    lineFileExpectWords(geneLf, 3, wordCount);
    stripChar(words[0], '"');
    unsigned short cix = (unsigned short)hashIntValDefault(geneIDToCix, words[0],
							   MAX_GENES_IN_GENOME+1);
    if (cix == MAX_GENES_IN_GENOME+1)
	lineFileAbort(geneLf, "genes ID %s not found in corr file", words[0]);
    unsigned int start = sqlUnsigned(words[2]);
    geneSetAdd(genesInGenome, cix, start);
    if (sameString(words[1], chrNum))
	geneSetAdd(genesOnChrom, cix, start);
    }
*retGenesOnChrom = genesOnChrom;
*retGenesInGenome = genesInGenome;
}

float calcANC(float **corrValues, unsigned short *cixs, int n)
/* Average the n gene indices' pairwise corrValues to get the 
 * average neighborhood correlation (ANC).  */
{
int i, j;
double total = 0.0;
for (i = 0;  i < n-1;  i++)
    for (j = i+1;  j < n;  j++)
	total += corrValues[cixs[i]][cixs[j]];
return (float)(2.0 * total / (float)(n * (n-1)));
}


void clusterOneWinSize(int winSize, float **corrValues, struct geneSet *gs,
		       struct clustersANCs *dist)
/* Compute clusters of genes in gs and their average neighborhood correlation (ANC) values
 * and store in pre-allocated dist. */
{
dist->n = 0;
dist->maxN = (gs->n * 2);
AllocArray(dist->ancs, dist->maxN);
AllocArray(dist->clusters, dist->maxN);

int startI = 0, endI = 0;
while (startI < gs->n)
    {
    unsigned int startPos = gs->starts[startI];
    unsigned int endPos = gs->starts[endI];
    while ((endPos < (startPos + winSize)) && (endI < gs->n-1))
	endPos = gs->starts[++endI];
    int clusterSize = endI - startI;
    if (clusterSize >= 2)
	{
	struct geneSet *newCluster = geneSetNew(clusterSize);
	int i;
	for (i = 0;  i < clusterSize;  i++)
	    geneSetAdd(newCluster, gs->cixs[startI + i], gs->starts[startI + i]);
	if (dist->n >= dist->maxN)
	    errAbort("dist has more than %d clusters -- add reallocation code here.",
		     dist->maxN);
	dist->clusters[dist->n] = newCluster;
	dist->ancs[dist->n] = calcANC(corrValues, newCluster->cixs, clusterSize);
	dist->n++;
	}
    startI++;
    }
}

struct clustersANCs *getDistros(float **corrValues, struct geneSet *genesOnChrom, 
			       int minBlkSize, int maxBlkSize)
/* For each window size, compute clusters of genes and their average neighborhood 
 * correlations (ANCs). */
{
int winCount = ((maxBlkSize - minBlkSize) / minBlkSize) + 1;
struct clustersANCs *dist;
AllocArray(dist, winCount);
int win, winSize;
for (win = 0, winSize = minBlkSize;  winSize <= maxBlkSize;  win++, winSize += minBlkSize)
    clusterOneWinSize(winSize, corrValues, genesOnChrom, &(dist[win]));
if (win != winCount)
    errAbort("arithmetic error: expected %d windows, got %d", winCount, win);
return dist;
}

int randUpTo(int n)
/* Return a random int from [0..n-1] (from man 3 rand / Numerical Recipes in C). */
{
return (int) ((float)n * (rand() / (RAND_MAX + 1.0)));
}

int *makeShuffleMap(int n)
/* Return an array of size n, containing shuffled [0..n-1]. */
{
int *shuffleMap;
AllocArray(shuffleMap, n);
unsigned char *taken;
AllocArray(taken, n);
int i;
for (i = 0;  i < n;  i++)
    {
    int randIx;
    while (taken[(randIx = randUpTo(n))]);
    shuffleMap[i] = randIx;
    taken[i] = TRUE;
    }
return shuffleMap;
}

int cmpFloatP(const void *a, const void *b)
/* Comparison function for qsort and binary search. */
{
float fa = *(float *)a, fb = *(float *)b;
if (fa < fb)
    return -1;
if (fa > fb)
    return 1;
return 0;
}

float **getRandomANCs(float **corrValues, int geneCount, struct clustersANCs *realDistros,
		      int winCount, int permCount)
/* For permCount iterations, shuffle the IDs of genesInGenome and use realDistros'
 * clusters to compute ANCs using the shuffled IDs for all clusters in each window size.
 * Return a 2-D array of per-window-size null-model ANCs, [win][dist->n], for
 * per-window-size p-value computation. */
{
float **randANCs;
AllocArray(randANCs, winCount);
unsigned short mappedCixs[MAX_GENES_IN_CLUSTER];
int p, win;
for (p = 0;  p < permCount;  p++)
    {
    int *shuffleMap = makeShuffleMap(geneCount);
    for (win = 0;  win < winCount;  win++)
	{
	float *randANCsInWin = randANCs[win];
	struct clustersANCs *winDist = &(realDistros[win]);
	if (p == 0)
	    {
	    AllocArray(randANCsInWin, (permCount * winDist->n));
	    randANCs[win] = randANCsInWin;
	    }
	int permOffset = (p * winDist->n);
	int c;
	for (c = 0;  c < winDist->n;  c++)
	    {
	    struct geneSet *cl = winDist->clusters[c];
	    if (cl->n > MAX_GENES_IN_CLUSTER)
		errAbort("Recompile with MAX_GENES_IN_CLUSTER > %d", MAX_GENES_IN_CLUSTER);
	    int i;
	    for (i = 0;  i < cl->n; i++)
		{
		unsigned short realCix = cl->cixs[i];
		mappedCixs[i] = shuffleMap[realCix];
		}
	    randANCsInWin[permOffset + c] = calcANC(corrValues, mappedCixs, cl->n);
	    }
	randANCs[win] = randANCsInWin;
	}
    }
for (win = 0;  win < winCount;  win++)
    {
    struct clustersANCs *winDist = &(realDistros[win]);
    qsort(randANCs[win], (permCount * winDist->n), sizeof(*randANCs[win]), cmpFloatP);
    }
return randANCs;
}


int findGreaters(float score, float *cmpScores, int cmpScoresCount)
/* Return the number of elements in sorted array cmpScores that are greater than score.
 * This is basically a binary search except we don't require an exact match, just the
 * index of the smallest member of cmpScores that is greater than score. */
{
int leftI = 0, rightI = cmpScoresCount-1;
while (rightI > leftI+1)
    {
    int midI = (rightI + leftI) / 2;
    if (cmpScores[midI] <= score)
	leftI = midI;
    else
	rightI = midI;
    }
// Just in case score is completely outside of the range of cmpScores:
if (cmpScores[leftI] > score)
    rightI = leftI;
if (cmpScores[rightI] <= score)
    rightI++;
return (cmpScoresCount - rightI);
}

void makePValues(FILE *pvalOutF, struct clustersANCs *realDistros, float **randANCs, int permCount,
		 char **geneIDs, int minBlkSize, int maxBlkSize)
/* Print out cluster ANCs, genes and those ANC's p-values (e.g. the proportion
 * of cluster ANCs derived from shuffled correlation values that are greater) */
{
fprintf(pvalOutF, " WindowSize \t BlockNumber \t ANC \t P-value \t NeighborhoodGenes\n");
int win, winSize;
for (win = 0, winSize = minBlkSize;  winSize <= maxBlkSize;  win++, winSize += minBlkSize)
    {
    struct clustersANCs *curDist = &(realDistros[win]);
    int totalANCs = curDist->n * permCount;
    int c, i;
    for (c = 0;  c < curDist->n;  c++)
	{
	struct geneSet *gs = curDist->clusters[c];
	float anc = curDist->ancs[c];
	float pvalue = (float)findGreaters(anc, randANCs[win], totalANCs) / (float)totalANCs;
	int windex = gs->starts[0] / minBlkSize;
	fprintf(pvalOutF, "%d \t %d \t %f \t %f \t", winSize, windex, anc, pvalue);
	for (i = 0;  i < gs->n;  i++)
	    fprintf(pvalOutF, " %s", geneIDs[gs->cixs[i]]);
	fprintf(pvalOutF, "\n");
	}
    }
}


void nhoodStatsBpairBg(char *corrFile, char *geneStartFile, char *chrNum,
		       int minBlkSize, int maxBlkSize, char *resultDir, int permCount)
/* nhoodStatsBpairBg - Compute avg neighborhood corr for genes across window sizes and 
 * p-val vs rand permuts. */
{
// Check all file/dir existence before launching into the real work:
struct lineFile *corrLf = lineFileOpen(corrFile, TRUE);
struct lineFile *geneLf = lineFileOpen(geneStartFile, TRUE);
makeDirsOnPath(resultDir);
char pvalOutFile[2048];
safef(pvalOutFile, sizeof(pvalOutFile), "%s/nhood_pvalues_bpair_c_CHR%s_MIN%d_MAX%d.txt",
      resultDir, chrNum, minBlkSize, maxBlkSize);
FILE *pvalOutF = mustOpen(pvalOutFile, "w");

// Get parallel arrays of gene IDs and pairwise correlation values, and a hash
// of gene IDs to indexes into those arrays, from corrFile.
char **geneIDs;
float **corrValues;
struct hash *geneIDToCix = NULL;
int geneCount = readCorrFile(corrLf, &geneIDs, &corrValues, &geneIDToCix);
lineFileClose(&corrLf);
(void)system("date");

// Get gene starts and indexes into geneIDs and corrValues from geneStartFile.
// FIXME: genesInGenome is currently unused because I'm shuffling the entire set of IDs.
struct geneSet *genesOnChrom, *genesInGenome;
readGeneFile(geneLf, chrNum, geneIDToCix, geneCount, &genesOnChrom, &genesInGenome);
lineFileClose(&geneLf);
if (genesOnChrom->n == 0)
    errAbort("%s contains no genes for %s", geneStartFile, chrNum);

// Compute gene clusters and real ANCs.
struct clustersANCs *realDistros = getDistros(corrValues, genesOnChrom, minBlkSize, maxBlkSize);

// Make the specified number of permutations of pairwise correlation values,
// and use those to calculate ANCs on all clusters, to get per-window-size null model 
// distributions.
int winCount = ((maxBlkSize - minBlkSize) / minBlkSize) + 1;
float **randANCs = getRandomANCs(corrValues, geneCount, realDistros, winCount, permCount);

// Use the per-win null models to estimate p-values for all clusters' ANCs, and write 
// out clusters, ANCs and p-values.
makePValues(pvalOutF, realDistros, randANCs, permCount, geneIDs, minBlkSize, maxBlkSize);
carefulClose(&pvalOutF);
(void)system("date");
}

int main(int argc, char *argv[])
/* Process command line. */
{
optionInit(&argc, argv, options);
if (argc != 7)
    usage();
permCount = optionInt("perms", permCount);
nhoodStatsBpairBg(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), argv[6], permCount);
return 0;
}
