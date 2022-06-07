/* Header treevolve.h     */
/* (c) Nick Grassly 1997  */
/* Dept. Zoology, Oxford. */

#ifndef NETEVOLVE_H
#define NETEVOLVE_H_

#include "Plugin.h"
#define MAX_NUMBER_REGIMES 10
#define MAX_NUMBER_SUBPOPS 100

#define MIN_NE 1.0e-15
#define BIG_NUMBER 1.0e+30
#define MAX_RATE_CATS 32
#define BRACKET 0.1
#define PROGRAM_NAME "SerialNetEvolve"
#define VERSION_NO 1.0
#define oops(s) { perror((s)); exit(EXIT_FAILURE); }
#define MALLOC(s,t) if(((s) = malloc(t)) == NULL) { oops("error: malloc() "); }
#define TOLERANCE 1.0e-10

/* Period parameters */  
//#define N 624
//#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */
#define VERSION_NO 1.0
/* Tempering parameters */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

#define FOR(i,n) for(i=0; i<n; i++)
#define FPN(file) fputc('\n', file)

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

#define BASE        2    /* base of floating point arithmetic */
#define DIGITS     53    /* no. of digits to the base BASE in the fraction */
#define MAXITER    30    /* max. no. of iterations to converge */

#define pos(i,j,n)      ((i)*(n)+(j))


typedef struct complexx complexx;
struct complexx {
   double re, im;
};
#define csize(a) (fabs(a.re)+fabs(a.im))

typedef struct Node Node;
struct Node {
	Node *daughters[2];
	double time;
	short *ancestral;
	char *sequence;
	short type;/*0=ca 1=re 2=tip 3=looked at*/
	int cutBefore;
	int deme;
	Node *next;
	Node *previous;
	Node *L_father;	/* Added by PB for MaxQueue/slimTree */
	Node *R_father;	/* Added by PB for MaxQueue/slimTree */
	Node *MaxNext;	/* Added by PB for MaxQueue/slimTree */
	char ID[10];	/* Added by PB */
	int Period;		/* Added by PB */
	int sampled;/* Added by PB 9/17/04*/
	Node *jumpright;/* Added by PB for slimTree*/
	Node *jumpleft;/* Added by PB for slimTree*/
};

typedef struct Regime Regime;
struct Regime {
	double t, N, e, r, m;
	int d;
};
typedef struct NodeInterval NodeInterval;
struct NodeInterval{
	double Start;
	double End;
	int Qty;
	int Time0;
	int Rec;
	int LookedAt;
	int Isolated;
};

enum {
        NoRates,
        CodonRates,
        GammaRates,
        DiscreteGammaRates
};
enum { A, C, G, T };
enum { F84, HKY, REV, numModels };


class NetEvolvePlugin : public Plugin {
	private:
	Node ***NodesPerPeriod; /* Added by PB 12/07/05 for conservative sampling */


int numBases, sampleSize;
double mutRate;
int ki[MAX_NUMBER_SUBPOPS], K, noRE, noCA, noMI;
double globTime, factr;
double genTimeInvVar;

int pCounter;



int SSizeperP=0; /* Sample Size per period*/
int PeriodsWanted=0; /* Number of Sampling Periods */
int StartAt=0; /* First Period to start sampling */
int ClassicTV=0;
int NoClock=0; /* Clock */
int outputFile=0;
int currPeriod2;
double iNodeProb=0;
int range=2;
int repNo;
char OFile[100];

Node *first, *avail;
char *ModelNames[3]={
	"F84",
	"HKY",
	"REV"
};

double intervalDis=0.0;
double prevTime=0.0;

int noPeriods;
Regime history[MAX_NUMBER_REGIMES];
int numIter, haploid;
int outputCoTimes;
char coTimeFile[50];
double gammaShape;
int numCats, rateHetero;
double catRate[MAX_RATE_CATS];
int currDeme;
int tipNo;

double freqRate[MAX_RATE_CATS];

char *nucleotides="ACGT";
double matrix[MAX_RATE_CATS][16];
double vector[4];
double *siteRates;
short *categories;

/* Quick fix Globals*/
int PeriodsChosen; /* Counts total number of time periods already chosen*/
 
int Pick;
int TotalRec;

Node **OUTStorage; /* To know where the BKP Nodes are, to free them later */
int slimNodes;
int slimIDs;

int MaxQty;
int numSampled;

typedef struct NodeInterval NodeInterval;
NodeInterval *Interval;
int IntervalPeriods;
int LastPeriod;

Node *Head_MaxQueue;
double probRE, probMI, Pn;
double lambda, N, migrationRate, recRate;


int model;

double freq[4], addFreq[4];
double freqR, freqY, freqAG, freqCT;
double freqA, freqC, freqG, freqT;
double tstv, kappa;

double beta, beta_A_R, beta_A_Y;
double tab1A, tab2A, tab3A;
double tab1C, tab2C, tab3C;
double tab1G, tab2G, tab3G;
double tab1T, tab2T, tab3T;
double mu, mu_kappa_1;

int i, j, l,result, currPeriod;
double maxRecRate;
char ppath[90],  pbase[100], temp[100];
FILE *coTimes;

double Rmat[6];
double Qij[16], Cijk[256], Root[4];
int mr;
int z_rndu=137;
unsigned w_rndu=13757;
unsigned long mt[624]; /* the array for the state vector  */
int mti=625; /* mti==N+1 means mt[N] is not initialized */

/*prototypes*/
	public:
	void input(std::string);
	void run();
	void output(std::string);
void Stack(Node *going);
long CalcGi(int deMe);
void Recombine(double t, int deme);
void Coalesce(double t, int deme);
void Migration(int deme, int numDemes);
/* -------------- prototypes ------------- */
void PrintTitle();
void PrintUsage();
void PrintParams();
int ImplementRegime(Regime* pp);
void ReadParams(std::string);
void ReadPeriod(FILE*);
int CheckPeriods(Regime *pp);
void ReadUntil(FILE *fv, char stopChar, char *what);
int ConstRoutine(Regime *pp);
double GenerateTimeEnd(Regime *pp);


int ConstSubRoutine(Regime *pp);
double GenTimeEndSub(Regime *pp);
void DistributeGenes(Regime *pp);
void SetModel(int model);
void SetCategories();
void RateHetMem();

/* protoypes */
void SeqEvolve();
void SeqEvolve(char*, char*);
void SeqEvolveSSTV();
void Mutate(Node *node,  int *nodeNo);
char SetBase(double *P);
void RandomSequence(char *seq);
void MutateSequence(char *seq, double len);
//void WriteSequence(Node *node, char *P);
void WriteSequence(int sampleNo, char *P,Node *node);
int EpiRoutine(Regime *pp);
double epifunc(double t);
double GenerateTimeEpi();
double approxFunc();
int EpiSubRoutine(Regime *pp);
double epiSubfunc(double t);
double GenTimeEpiSub(Regime *pp);
double approxFuncSub();
void InstantCoalesce(double t, Regime *pp);
void PrintNetworkOrTree(FILE *fv, Node *node, double maxTime);
void PrintNode(FILE *treeFile, Node *node, double maxTime);
void StorePeriods1b(Node *node);
void Enqueue(Node *maxNode, Node *Child);
void PicknAssignID(Node *node,  int *nodeNo, int BKP, int recParentBKP);
void EnlargeIntervalDimension(int Period);
void ReID(Node *node, int Number);
	void FixInternalNodeClockTime(Node *node, double time);
	void FixClockPeriods();
void heapify_Nodelist ( Node **list , int newnode );
void heapsort_Nodelist ( Node **list, int last );
Node *FirstNodePop();
void PickPeriods2();
Node *NodePop(Node *prev);
void SetMatrix(double *matrix, double len);
void SetVector(double *vector, short base, double len);
int EigenREV (double Root[], double Cijk[]);
void SetSeed (unsigned long seed);
double rndu (void);

double rndgamma (double s);
double zeroinA(double ax, double bx);//, double (*f)(double x));
double zeroinB(double ax, double bx);//, double (*f)(double x));
int DiscreteGamma (double freqK[], double rK[], double alfa, double beta, int K, int median);
double rndgamma1 (double s);
double rndgamma2 (double s);
double LnGamma (double alpha);
double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);
double PointNormal (double prob);
double PointChi2 (double prob, double v);
int abyx (double a, double x[], int n);
int xtoy (double x[], double y[], int n);
int matinv( double x[], int n, int m, double space[]);
int eigen(int job, double A[], int n, double rr[], double ri[],
          double vr[], double vi[], double w[]);
void balance(double mat[], int n, int *low, int *hi, double scale[]);
void unbalance(int n, double vr[], double vi[], int low, int hi,
               double scale[]);
int realeig(int job, double mat[], int n,int low, int hi, double valr[],
            double vali[], double vr[], double vi[]);
void elemhess(int job, double mat[], int n, int low, int hi, 
            double vr[], double vi[], int work[]);

complexx myc(double re,double im);
complexx conj (complexx a);
complexx cplus (complexx a, complexx b);
complexx cminus (complexx a, complexx b);
complexx cby (complexx a, complexx b);
complexx cdiv (complexx a,complexx b);
complexx cexp (complexx a);
complexx cfactor (complexx x, double a);
int cxtoy (complexx x[], complexx y[], int n);
int cmatby (complexx a[], complexx b[], complexx c[], int n,int m,int k);
int cmatout (FILE * fout, complexx x[], int n, int m);
int cmatinv( complexx x[], int n, int m, double space[]);
void SaveNodeinPeriodsArray(Node *dec1, double t);
void SeqEvolveSSTV (char * pbase, char *ppath);
void FixSlimAncestors(Node *node);
void StorePeriods1();
Node* BuildSlimTree(Node*);
void PickPeriods();
	};

#endif /* _TREEVOLVE_H_ */
