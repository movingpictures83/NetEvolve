/*---------------------------------------*/
/* TREEVOLVE v1.3 - main source file     */
/*                                       */
/*  The coalescent process with:         */
/*      recombination                    */
/*      population subdivision           */
/*      exponential growth               */
/*                                       */
/* 1997 (c) Nick Grassly                 */
/*          Dept. of Zoology             */
/*          Oxford. OX1 3PS              */
/*      nicholas.grassly@zoo.ox.ac.uk    */
/*      http://evolve.zoo.ox.ac.uk/      */
/*---------------------------------------*/
/* 2006 NetEvolve Modifications by		 */
/*		Patricia Buendia				 */
/*		http://biorg.cs.fiu.edu/ 		 */
/*---------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include "PluginProxy.h"
#include "NetEvolvePlugin.h"
#ifdef __MWERKS__
#include <console.h>
#endif	

//#define PROGRAM_NAME "SerialNetEvolve"
#define VERSION_NO 1.0


/* -------------- variables --------------- */
//int numBases, sampleSize;
//double mutRate;
//int SSizeperP=0; /* Sample Size per period*/
//int PeriodsWanted=0; /* Number of Sampling Periods */
//int StartAt=0; /* First Period to start sampling */
//int ClassicTV=0;
//int NoClock=0; /* Clock */
//int pCounter;
//int outputFile=0;
//int currPeriod2;
//double iNodeProb=0;
//int range=2;
//int repNo;
//char OFile[100];
//Node ***NodesPerPeriod;

//Node *first, *avail;
//char *ModelNames[numModels]={
//	"F84",
//	"HKY",
//	"REV"
//};
//int ki[MAX_NUMBER_SUBPOPS], K, noRE, noCA, noMI;
//double globTime, factr;
//double genTimeInvVar;

//double intervalDis=0.0;
//double prevTime=0.0;

//static int noPeriods;
//static Regime history[MAX_NUMBER_REGIMES];
//static int numIter, haploid;
//static int outputCoTimes;
//static char coTimeFile[50];

/* -------------- prototypes ------------- */
//static void PrintTitle();
//static void PrintUsage();
//static void PrintParams();
//int ImplementRegime(Regime* pp);
//void ReadParams(int argc, char **argv);
//void ReadPeriod();
//int CheckPeriods(Regime *pp);
//void ReadUntil(FILE *fv, char stopChar, char *what);


#define oops(s) { perror((s)); exit(EXIT_FAILURE); }
#define MALLOC(s,t) if(((s) = malloc(t)) == NULL) { oops("error: malloc() "); }

/*-----------------------------------------*/
void NetEvolvePlugin::input(std::string inputfile)
{
/*int i, j, l,result, currPeriod;
double maxRecRate;
char ppath[90],  pbase[100], temp[100];
FILE *coTimes;
*/
/*#ifdef __MWERKS__
	argc=ccommand(&argv);
#endif	*/

	if (setvbuf( stderr, NULL, _IOLBF, 512 ) ||
		setvbuf( stdout, NULL, _IOLBF, 512 ) ) {
		fprintf(stderr, "Failed to buffer stdout and stderr\n");
		exit(0);
	}
	//long now;
	//time(&now);
	//SetSeed(time(NULL));
	SetSeed(1234);
	avail=NULL;
	PrintTitle();
	ReadParams(inputfile);
	maxRecRate=0;
	for(i=0;i<noPeriods;i++){
		if(CheckPeriods(&history[i])==0) exit(0);
		if(history[i].r>0 && maxRecRate==0)
			maxRecRate=history[i].r;
		else if(maxRecRate<history[i].r)
			maxRecRate=history[i].r;
	}
	PrintParams();
}

void NetEvolvePlugin::run() {
	if(StartAt==0 && PeriodsWanted>0 && SSizeperP>0){
		if(maxRecRate==0)
			range=4;
		else{
			range=8+(int)(log10(maxRecRate/mutRate));
			if(range<2)
				range=2;
		}
		sampleSize=PeriodsWanted*SSizeperP;
		MALLOC(NodesPerPeriod, sizeof(Node **) * (sampleSize*range));/* nodes in slimTree */
		for (i=0; i < sampleSize*range; i++) 
			MALLOC(NodesPerPeriod[i], sizeof(Node *) * (sampleSize*range));/* (recombination adds nodes:How many? should be prop. to rec rate) */

	}
	else if(PeriodsWanted==0 ||  SSizeperP==0)
		ClassicTV=1;
}

void NetEvolvePlugin::output(std::string outfile) {
	if(outputCoTimes){
		if( (coTimes=fopen(coTimeFile, "w"))==NULL){
			fprintf(stderr, "Failed to open coalescent times file\n\n");
			exit(0);
		}
		fprintf(stderr, "Coalescent times output to file: %s\n\n\n", coTimeFile);
	}
	RateHetMem();
	SetModel(model);
	for(repNo=1;repNo<=numIter;repNo++){
		//if (maxRecRate>0)
		//	printf ("Creating Recombinant Network #%d... Please wait\n",repNo);
		//else
		//	printf ("Creating Tree #%d... Please wait\n",repNo); Error: Prints to stdout
		if(!ClassicTV){
			for (i=0; i < sampleSize*range; i++) {
				for(j=0;j<sampleSize*range ;j++)
					NodesPerPeriod[i][j]=NULL;
			}
		}

		K=sampleSize;
		first=FirstNodePop();/*memory allocation for first node*/
		first->type=2;
		first->time=0.0;
		first->sampled=1;
		first->Period=1;
		first->deme=0;
		for(i=1;i<sampleSize;i++){
			first=NodePop(first);/*memory allocation for rest of sample at t=0*/
			first->type=2;
			first->time=0.0;
			first->sampled=1;
			first->Period=1;
			first->deme=0;
		}
		globTime=0.0;
		noRE=noCA=noMI=0;
		currPeriod2=0;
		currPeriod=0;
		pCounter=0;
		result=1;
		while(currPeriod<noPeriods && result==1){
			result=ImplementRegime(&history[currPeriod]);
			fprintf(stderr, "Coalescent calculations finished for period %d\n", (currPeriod+1));
			if(result==0){
				fprintf(stderr, "All coalescences have now occurred at time %f\n\n", globTime);
				if(outputCoTimes)
					fprintf(coTimes, "%f\n", globTime);
			}else
				fprintf(stderr, "No. of extant lineages = %d\n", K);
			currPeriod++;
		}
		fprintf(stderr, "No. of coalescent events    = %d\n", noCA);
		fprintf(stderr, "No. of recombination events = %d\n", noRE);
		fprintf(stderr, "No. of migrations           = %d\n", noMI);
		if(outputFile){
			j=strlen(OFile);
			ppath[0]='\0';
			l=0;
			for(i=0;i<j;i++){
				temp[l++]=OFile[i];
				if(OFile[i]=='\\'){
					temp[l]='\0';
					strcat(ppath,temp);
					l=0;
				} else if(OFile[i]=='.'){
					temp[l-1]='\0';
					strcpy(pbase,ppath);
					strcat(pbase,temp);
					break;
				}

			}
			sprintf(temp,"%s%d.txt",pbase,repNo);
			remove(temp);
		}
		else{
			pbase[0]='\0';
			ppath[0]='\0';

		}
		if(ClassicTV)
			SeqEvolve(pbase,ppath);
		else
			SeqEvolveSSTV(pbase,ppath);
		fprintf(stderr, "Sequence evolution finished for replicate %d\n\n", repNo);
	}
	if(outputCoTimes)
		fclose(coTimes);
	if (!ClassicTV) {
		for (i = 0; i < PeriodsWanted*2; i++) 
			free(NodesPerPeriod[i]);
		free(NodesPerPeriod);
	}

	fprintf(stderr, "it's all over! -\n	sequences written to file\n");
	exit(1);
}

int NetEvolvePlugin::ImplementRegime(Regime* pp)
{
	if(pp->e==0.0){
		if(pp->d==1)
			return(ConstRoutine(pp));
		else
			return(ConstSubRoutine(pp));
	}else{
		if(pp->d==1)
			return(EpiRoutine(pp));
		else
			return(EpiSubRoutine(pp));
	}
}

void NetEvolvePlugin::PrintTitle()
{
	fprintf(stderr, "%s %.2f - Serial Coalescent Simulation and Sequence Evolution\n",PROGRAM_NAME ,VERSION_NO);
	fprintf(stderr, "--------------------------------------------------------------------------\n\n");
}

void NetEvolvePlugin::PrintUsage()
{
	fprintf(stderr, "Usage: %s -options <PARAMETER.FILE >OUPUT_FILE\n", PROGRAM_NAME);
	fprintf(stderr, "\tFor a description of options available please\n");
	fprintf(stderr, "\trefer to manual included in package\n\n");
	exit(0);
}

void NetEvolvePlugin::ReadParams(std::string infile)
{
	int  j;
	char ch, str[255], *str2;
	//char *P;
	
	numBases=500;
	sampleSize=20;
	mutRate=0.0000001;

	numCats=1;
	rateHetero=NoRates;
	catRate[0]=1.0;
	gammaShape=1.0;	
	Rmat[0]=Rmat[1]=Rmat[2]=Rmat[3]=Rmat[4]=1.0;

	freq[0]=freq[1]=freq[2]=freq[3]=0.25;
	tstv=2.0;
	model=F84;
	numIter=1;
	haploid=0;
	outputCoTimes=0;
	genTimeInvVar=1.0;
	FILE* fp =fopen(infile.c_str(), "r");
	
	if(feof(fp)){
		fprintf(stderr, "Unable to read parameters from fp\n");
		exit(0);
	}
	//str2="BEGIN TVBLOCK";
	fgets(str, 255, fp);
	str2=strstr(str,"TVBLOCK");

        // GN added next 3 statements
        //j = strlen(str);
       // str[j - 1] = '\0';
       // str[j - 2] = '\0';
	//if(strcmp(str, str2)!=0){//Problems in Linux with newline chars
	if(str2==NULL){
		fprintf(stderr, "Input does not contain a paramters block: %s\n",str);
		exit(0);
	}
	ch=fgetc(fp);
	while(isspace(ch))
		ch=fgetc(fp);
	while(ch=='['){
		ReadUntil(fp, ']', "closing bracket");
		ch=fgetc(fp);
		while(isspace(ch))
			ch=fgetc(fp);
	}
	noPeriods=0;
	
	while(!feof(fp)){
		if(ch=='*'){
			ReadPeriod(fp);
			noPeriods++;
		}else{
			ch=toupper(ch);
			switch (ch) {
				case 'A':
					if (rateHetero==CodonRates) {
						fprintf(stderr, "You can only have codon rates or gamma rates not both\n");
						exit(0);
					}
					if (rateHetero==NoRates)
						rateHetero=GammaRates;
					if (fscanf(fp, "%lf", &gammaShape)!=1 || gammaShape<=0.0) {
						fprintf(stderr, "Bad Gamma Shape\n");
						exit(0);
					}
				break;
				case 'B':
					if (fscanf(fp, "%lf", &genTimeInvVar)!=1) {
						fprintf(stderr, "Bad compound generation time parameter\n\n");
						PrintUsage();
					}
				break;
				case 'C':
					if (rateHetero==GammaRates) {
						fprintf(stderr, "You can only have codon rates or gamma rates not both\n");
						exit(0);
					}
					numCats=3;
					rateHetero=CodonRates;
					if (fscanf(fp, "%lf,%lf,%lf", &catRate[0], &catRate[1], &catRate[2])!=3) {
						fprintf(stderr, "Bad codon-specific rates\n\n");
						PrintUsage();
					}
				break;
				case 'D':
					if (fscanf(fp, "%d", &StartAt)!=1) {
						fprintf(stderr, "Bad [Starting Period] wanted number\n\n");
						PrintUsage();
					}
				break;
				case 'F':
					if (fscanf(fp, "%lf,%lf,%lf,%lf", &freq[0], &freq[1], 
													&freq[2], &freq[3])!=4) {
						fprintf(stderr, "Bad Base Frequencies\n\n");
						PrintUsage();
					}
				break;
				case 'G':
					if (rateHetero==CodonRates) {
						fprintf(stderr, "You can only have codon rates or gamma rates not both\n");
						exit(0);
					}
					
					rateHetero=DiscreteGammaRates;
					if ((fscanf(fp, "%d", &numCats))!=1 || numCats<2 || numCats>MAX_RATE_CATS) {
						fprintf(stderr, "Bad number of Gamma Categories\n");
						exit(0);
					}
				break;
				case 'H':
					haploid=1;
					while(!isspace(ch))
						ch=getc(fp);
				break;
				case 'I':
					if (fscanf(fp, "%lf", &iNodeProb)!=1|| iNodeProb<0) {
						fprintf(stderr, "Bad [Internal Nodes Sampling Probability] number\n\n");
						PrintUsage();
					}
				break;
				case 'K':
					NoClock=1;
					while(!isspace(ch))
						ch=getc(fp);
				break;
				
				case 'L':
					if (fscanf(fp, "%d", &numBases)!=1 || numBases<1) {
						fprintf(stderr, "Bad sequence length\n\n");
						PrintUsage();
					}
				break;
				case 'N':
					if (fscanf(fp, "%d", &numIter)!=1 || numIter <1) {
						fprintf(stderr, "Bad number of replicates\n\n");
						PrintUsage();
					}
				break;
				case 'M':
					if (fscanf(fp, "%lf", &intervalDis)!=1 || intervalDis <0) {
						fprintf(stderr, "Bad symmetrical interval distance\n\n");
						PrintUsage();
					}
				break;
				case 'O':
					outputCoTimes=1;
					ch=fgetc(fp);
					if(isspace(ch))
						strcpy(coTimeFile, "coalescent.times");
					else{
						j=0;
						do{
							coTimeFile[j]=ch;
							j++;
							ch=fgetc(fp);
						}while(!isspace(ch));
						coTimeFile[j]='\0';
					}
				break;
				case 'P':
					if (fscanf(fp, "%d", &PeriodsWanted)!=1 || PeriodsWanted <2) {
						fprintf(stderr, "Bad [Periods wanted] number\n\n");
						PrintUsage();
					}
				break;
				case 'R':
					if (fscanf(fp, "%lf,%lf,%lf,%lf,%lf,%lf", &Rmat[0], &Rmat[1], 
													&Rmat[2], &Rmat[3], &Rmat[4], &Rmat[5])!=6) {
						fprintf(stderr, "Bad general rate matrix\n\n");
						PrintUsage();
					}
					if (Rmat[5]!=1.0) {
						for (j=0; j<5; j++) 
							Rmat[j]/=Rmat[5];
						Rmat[5]=1.0;
					}
				break;
				case 'S':
					if (fscanf(fp, "%d", &sampleSize)!=1 || sampleSize<1) {
						fprintf(stderr, "Bad sample size\n\n");
						PrintUsage();
					}
				break;
				case 'U':
					if (fscanf(fp, "%lf", &mutRate)!=1) {
						fprintf(stderr, "Bad mutation rate\n\n");
						PrintUsage();
					}
				break;
				case 'T':
					if (fscanf(fp, "%lf", &tstv)!=1) {
						fprintf(stderr, "Bad Ts/Tv ratio\n\n");
						PrintUsage();
					}
				break;
				case 'V':
					model=-1;
					fgets(str, 4, fp);
					for (j=F84; j<numModels; j++) {
						if (strcmp(str, ModelNames[j])==0)
							model=j;
					}
					if (model==-1) {
						fprintf(stderr, "Unknown Model: %s\n\n", str);
						PrintUsage();
						exit(0);
					}
				break;
				case 'Z':
					if (fscanf(fp, "%d", &SSizeperP)!=1|| SSizeperP <2) {
						fprintf(stderr, "Bad [Sampling Size Per Period] number\n\n");
						PrintUsage();
					}
				break;
				case 'X':
					outputFile=1;
					ch=fgetc(fp);
					if(isspace(ch))
						outputFile=0;
					else{
						j=0;
						if(ch=='"'){
							ch=fgetc(fp);
							do{
								OFile[j]=ch;
								j++;
								ch=fgetc(fp);
							}while(ch!='"');
						}
						else
							fprintf(stderr, "Output file needs to be enclosed in quotation marks\n\n");

					
						OFile[j]='\0';
					}

				break;

				default :
					fprintf(stderr, "Incorrect parameter: %c\n\n", ch);
					PrintUsage();
				break;
			}
		}
		ch=fgetc(fp);
		while(isspace(ch) && !feof(fp))
			ch=fgetc(fp);
		while(ch=='['){
			ReadUntil(fp, ']', "closing bracket");
			ch=fgetc(fp);
			while(isspace(ch))
				ch=fgetc(fp);
		}
	}
		
}

void NetEvolvePlugin::ReadPeriod(FILE* fp)
{
int periodNo, i;
char ch, str[255];
	
	ch=fgetc(fp);
	if(ch!='P'){
		fprintf(stderr, "Asterisks denote period data\n");
		exit(0);
	}
	while(!isspace(ch))
		ch=fgetc(fp);
	if(fscanf(fp, "%d", &periodNo)!=1){
		fprintf(stderr, "Periods must have numerical label\n");
		exit(0);
	}
	periodNo--;
	history[periodNo].t=-1.0;	/* set defaults */
	history[periodNo].N=1000000;
	history[periodNo].e=0.0;
	history[periodNo].d=1;
	history[periodNo].m=0.0;
	history[periodNo].r=0.0;
	ch=fgetc(fp);
	while(isspace(ch))
		ch=fgetc(fp);
	while(ch=='['){
		ReadUntil(fp, ']', "closing bracket");
		ch=fgetc(fp);
		while(isspace(ch))
			ch=fgetc(fp);
	}
	while(ch!='*'){
		ch=toupper(ch);
		switch (ch) {
			case 'T':
				if (fscanf(fp, "%lf", &history[periodNo].t)!=1) {
					fprintf(stderr, "Bad length of period %d\n\n", periodNo+1);
					PrintUsage();
				}
			break;
			case 'N':
				ch=fgetc(fp);
				i=0;
				while(!isspace(ch)){
					str[i]=ch;
					i++;
					ch=fgetc(fp);
				}
				str[i]='\n';
				if (sscanf(str, "%lf", &history[periodNo].N)!=1) {
					str[0]=toupper(str[0]);
					if(str[0]=='P'){
						history[periodNo].N=((history[periodNo-1].N)*exp(-(history[periodNo-1].e)*(history[periodNo-1].t)));
						if(history[periodNo].N<MIN_NE){
							history[periodNo].N=MIN_NE;
							fprintf(stderr, "Period %d Ne too small. Therefore setting to %e\n", periodNo, MIN_NE);
						}
					}else{
						fprintf(stderr, "Bad population size for period %d\n\n", periodNo+1);
						PrintUsage();
					}
				}
			break;
			case 'E':
				if (fscanf(fp, "%lf", &history[periodNo].e)!=1) {
					fprintf(stderr, "Bad exponential growth for period %d\n\n", periodNo+1);
					PrintUsage();
				}
			break;
			case 'D':
				if (fscanf(fp, "%d", &history[periodNo].d)!=1) {
					fprintf(stderr, "Bad # sub populations for period %d\n\n", periodNo+1);
					PrintUsage();
				}
			break;
			case 'M':
				if (fscanf(fp, "%lf", &history[periodNo].m)!=1) {
					fprintf(stderr, "Bad migration rate for period %d\n\n", periodNo+1);
					PrintUsage();
				}
			break;
			case 'R':
				if (fscanf(fp, "%lf", &history[periodNo].r)!=1) {
					fprintf(stderr, "Bad recombination rate for period %d\n\n", periodNo+1);
					PrintUsage();
				}
			break;
			default :
				fprintf(stderr, "Incorrect period parameter: %c\n\n", ch);
				PrintUsage();
			break;
		}
		
		ch=fgetc(fp);
		while(isspace(ch) && !feof(fp))
			ch=fgetc(fp);
		while(ch=='['){
			ReadUntil(fp, ']', "closing bracket");
			ch=fgetc(fp);
			while(isspace(ch))
				ch=fgetc(fp);
		}
	}
	while(!feof(fp) && !isspace(ch))
		ch=fgetc(fp);
	
}


void NetEvolvePlugin::ReadUntil(FILE *fv, char stopChar, char *what)
{
	char ch;
	
	ch=fgetc(fv);
	while (!feof(fv) && ch!=stopChar) 
		ch=fgetc(fv);

	if (feof(fv) || ch!=stopChar) {
		fprintf(stderr, "%s missing", what);
		exit(0);
	}
}

int NetEvolvePlugin::CheckPeriods(Regime *pp)
{

	if(pp->d==1 && pp->m!=0.0){
		fprintf(stderr, "The migration rate has no meaning if the number of subpopulations is one\n");
		return 0;
	}
	return 1;	
}

void NetEvolvePlugin::PrintParams()
{
int i;
double t, N;

	fprintf(stderr, "sequence length      = %d\n", numBases);
	fprintf(stderr, "sample size          = %d\n", sampleSize);
	fprintf(stderr, "mutation rate u      = %e\n", mutRate);
	fprintf(stderr, "number of replicates = %d\n", numIter);
	fprintf(stderr, "substitution model   = %s\n", ModelNames[model]);
	if (rateHetero==CodonRates) {
		fprintf(stderr, "Codon position rate heterogeneity:\n");
		fprintf(stderr, "    rates = 1:%f 2:%f 3:%f\n", catRate[0], catRate[1], catRate[2]);
	} else if (rateHetero==GammaRates) {
		fprintf(stderr, "Continuous gamma rate heterogeneity:\n");
		fprintf(stderr, "    shape = %f\n", gammaShape);
	} else if (rateHetero==DiscreteGammaRates) {
		fprintf(stderr, "Discrete gamma rate heterogeneity:\n");
		fprintf(stderr, "    shape = %f, %d categories\n", gammaShape, numCats);
	} else
		fprintf(stderr, "Rate homogeneity of sites.\n");
	fprintf(stderr, "Model=%s\n", ModelNames[model]);
	if (model==F84)
		fprintf(stderr, "  transition/transversion ratio = %G (K=%G)\n", tstv, kappa);
	else if (model==HKY)
		fprintf(stderr, "  transition/transversion ratio = %G (kappa=%G)\n", tstv, kappa);
	else if (model==REV) {
		fprintf(stderr, "  rate matrix = gamma1:%7.4f alpha1:%7.4f  beta1:%7.4f\n", Rmat[0], Rmat[1], Rmat[2]);
		fprintf(stderr, "                                beta2:%7.4f alpha2:%7.4f\n", Rmat[3], Rmat[4]);
		fprintf(stderr, "                                              gamma2:%7.4f\n", Rmat[5]);
	}
	if(haploid){
		fprintf(stderr, "haploid model implemented\n");
		factr=2.0;
	}else{
		fprintf(stderr, "diploid model implemented\n");
		factr=4.0;
	}
	fprintf(stderr, "Generation time / variance in offspring number = %f\n", genTimeInvVar);
	if(genTimeInvVar==1.0)
		fprintf(stderr, "\t- corresponds to Wright-Fisher model of reproduction\n");
	fprintf(stderr, "\nPopulation Dynamic Periods:\n");
	fprintf(stderr, "---------------------------\n");
	t=0.0;
	for(i=0;i<noPeriods;i++){
		fprintf(stderr, "Period %d\n", i+1);
		if(history[i].t>0.0)
			fprintf(stderr, "Length: %f\n",  history[i].t);
		else
			fprintf(stderr, "Period running until final coalescence\n",  history[i].t);
		fprintf(stderr, "Time at start: %f\n", t);
		if(history[i].d>1){
			fprintf(stderr, "Population subdivided into %d demes\n", history[i].d);
			fprintf(stderr, "Deme size at t = %f is: %f\n", t,history[i].N);
			fprintf(stderr, "Migration rate: %e\n", history[i].m);
			if(history[i].e!=0.0){
				fprintf(stderr, "Exponential growth (decline backwards) at rate: %f\n", history[i].e);
				if(history[i].t>0.0){
					N=((history[i].N)*exp(-(history[i].e)*(history[i].t)));
					fprintf(stderr, "Expected deme size at end of period: %f\n", N);
				}
			}
		}else{
			fprintf(stderr, "Population panmictic with size %f\n", history[i].N);
			if(history[i].e!=0.0){
				fprintf(stderr, "Exponential growth (decline backwards) at rate: %f\n", history[i].e);
				if(history[i].t>0.0){
					N=((history[i].N)*exp(-(history[i].e)*(history[i].t)));
					fprintf(stderr, "Expected population size at end of period: %f\n", N);
				}
			}
		}
		fprintf(stderr, "Recombination rate: %e\n\n", history[i].r);
		fprintf(stderr, "---------------------------\n");
		t+=history[i].t;
	}
}


long NetEvolvePlugin::CalcGi(int deMe)
{
int i, j, k, posn;
long count;
short *ptr1, *ptr2;
Node *nptr;

	count=0;
	nptr=first;
	for(i=0;i<K;i++){
		if(nptr->deme==deMe){/* i.e. if in same deme add Gi */
			ptr1=nptr->ancestral;
			posn=0;
			while(posn<numBases && *ptr1==0){
				ptr1++;
				posn++;
			}
			
			if(posn<(numBases-1)){
				for(j=(posn+1);j<numBases;j++){
					ptr2=ptr1;
					for(k=j;k<numBases;k++){
						ptr2++;
						if(*ptr2==1){
							count++;
							break;
						}
					}
					ptr1++;
				}
			}
		}
		nptr=nptr->next;/*move to next gamete*/
	}
	return count;
}
/********************************************************************************************************/
/*** PB: An array of pointers to tree nodes needs to be allocated with size Periods*2 and SsizePP*2   ***/
/*** Each row identifies a sampling time period. Nodes are added backwards, so last periods first     ***/
/********************************************************************************************************/
void NetEvolvePlugin::SaveNodeinPeriodsArray(Node *dec1, double t)
{
int i, j, newP;
int end_iNodes;
double rnd;

	if(dec1->time==0){
		dec1->daughters[0]=NULL;
		if(pCounter>=SSizeperP ){
			pCounter=0;
			currPeriod2++;
		}
		newP=currPeriod2;
		dec1->Period=currPeriod2;
		pCounter++;

	}
	else{
		if(dec1->type==1)
			newP=dec1->daughters[0]->Period+1;
		else{
			if(dec1->daughters[0]->Period>dec1->daughters[1]->Period)
				newP=dec1->daughters[0]->Period+1;
			else
				newP=dec1->daughters[1]->Period+1;
		}
		dec1->Period=newP;

	}
	i=0;end_iNodes=0;
	while (NodesPerPeriod[newP][i]!=NULL) {/* How many internal nodes*/
		if(NodesPerPeriod[newP][i++]->daughters[0]!=NULL)
			end_iNodes++;
	}
	i=end_iNodes;
	if(newP<=currPeriod2 && dec1->time!=0){/* To have all internal nodes to the left in the array*/
		j=end_iNodes+SSizeperP;
		for(j;j>end_iNodes;j--)
			NodesPerPeriod[newP][j]=NodesPerPeriod[newP][j-1];
		
	}
	else{
		while (NodesPerPeriod[newP][i]!=NULL) 
			i++;
	}

	//if(dec1->time==0 && (NoClock ||i==0||(i==end_iNodes && iNodeProb==0))){/* For Leaves only */
	if(dec1->time==0 && (NoClock ||i==0||i==end_iNodes )){/* For Leaves only */
		rnd=rndu();
		if(!NoClock && newP>0){/* for first leaf */
				if(iNodeProb>0){
					if(NodesPerPeriod[newP+1][0]!=NULL)/* new time should be below that of next period first internal node time*/
						t=NodesPerPeriod[newP+1][0]->time-1;
				}
				if(intervalDis>0)
					dec1->time=prevTime+intervalDis;
				else
					dec1->time=prevTime+rnd*(t-prevTime);
		}
		else
			dec1->time=rnd*t;
		if(dec1->time>t)
			printf("Problem with time asiggnmetn!\n");
		prevTime=dec1->time;
	}
	else if(!NoClock){/* If Clock */ 
		if(	i>0){
			if(dec1->time==0)
				dec1->time=NodesPerPeriod[newP][i-1]->time; /*  leaves have same time like other leaves*/
			if(dec1->time>t)
				printf("Problem with time asiggnmetn!\n");
		}
	}
			
	NodesPerPeriod[newP][i]=dec1;
}

/*-------------------------------------------------------------------------------*/
void NetEvolvePlugin::Recombine(double t, int deme)
{
int i, picked, sum1, sum2;
double rnd;
Node *rec, *anc1, *anc2;

	do{
		do{
			rnd=rndu();
			picked=(int) ( rnd*(ki[deme]) );
		}while(picked==ki[deme]);
		
		rec=first;
		i=0;
		while(rec->deme!=deme)/*pick recombinant from correct deme*/
			rec=rec->next;
		while(i<picked){
			rec=rec->next;
			while(rec->deme!=deme)
				rec=rec->next;
			i++;
		}
		
		do{
			rnd=rndu();
			picked=(int) (rnd*numBases);
		}while(picked==numBases || picked==0);/*cuts sequence at one of 
												m-1 possible points*/
		sum1=sum2=0;
		for(i=0;i<picked;i++)			/*checks ancestral tuples*/
			sum1+=rec->ancestral[i];
		for(i=picked;i<numBases;i++)
			sum2+=rec->ancestral[i];
			
	}while(sum1==0 || sum2==0);
	if(!ClassicTV)
		SaveNodeinPeriodsArray(rec, t);


	
	anc1=first=NodePop(first);			/*memory allocation*/
	anc2=first=NodePop(first);			/*memory allocation*/
	rec->previous->next=rec->next;		/*maintain loop    */
	rec->next->previous=rec->previous;	/*maintain loop    */
	
	anc1->daughters[0]=rec;				/*point two ancestral gametes to recombinant*/
	anc2->daughters[0]=rec;
	anc1->daughters[1]=anc2;			/*point RE ancestors to each other using spare pointer*/
	anc2->daughters[1]=anc1;
	anc1->time=t;						/*record time of event*/
	anc2->time=t;
	anc1->sampled=0;						/*record time of event*/
	anc2->sampled=0;
	anc1->Period=0;						/*record period*/
	anc2->Period=0;
	anc1->type=1;						/*recombinant node*/
	anc2->type=1;
	anc1->cutBefore=picked;				/*record cut posn*/
	anc2->cutBefore=-1;
	for(i=0;i<picked;i++){				/*sets ancestral tuples*/
		anc1->ancestral[i]=rec->ancestral[i];
		anc2->ancestral[i]=0;
	}
	for(i=picked;i<numBases;i++){
		anc2->ancestral[i]=rec->ancestral[i];
		anc1->ancestral[i]=0;
	}
	anc1->deme=deme;					/*record which deme the ancestors are in*/
	anc2->deme=deme;
	if(!ClassicTV){						/* Serial Sampling */
		SaveNodeinPeriodsArray(anc1, t);
		SaveNodeinPeriodsArray(anc2, t);
	}
	else if(NoClock){					/* No Clock */
		if(rec->time==0){
				rnd=rndu();
				rec->time=rnd*t;
		}

	}

}



/*---------------------------------------------------------------------------------------*/
void NetEvolvePlugin::Coalesce(double t, int deme)
{
int i, picked1, picked2;
double rnd;
Node *dec1, *dec2;
short *p, *q, *r;


	do{						/*choose CA candidates from correct deme*/
		rnd=rndu();
		picked1=(int) ( rnd*(ki[deme]) );
	}while(picked1==ki[deme]);
	do{
		rnd=rndu();
		picked2=(int) ( rnd*(ki[deme]) );
	}while(picked2==ki[deme] || picked2==picked1);
	
	dec1=dec2=first;
	i=0;
	while(dec1->deme!=deme)/*pick CA candidates from loop*/
		dec1=dec1->next;
	while(i<picked1){
		dec1=dec1->next;
		while(dec1->deme!=deme)
			dec1=dec1->next;
		i++;
	}
	i=0;
	while(dec2->deme!=deme)/*pick CA candidates from loop*/
		dec2=dec2->next;
	while(i<picked2){
		dec2=dec2->next;
		while(dec2->deme!=deme)
			dec2=dec2->next;
		i++;
	}
	if(!ClassicTV){ /*PB: Modify time of leaves and save in periods array*/
		if(dec1->time==0)
			SaveNodeinPeriodsArray(dec1, t);
		if(dec2->time==0)
			SaveNodeinPeriodsArray(dec2, t);
	}
	else if(NoClock){/* No Clock but classic Treevolve*/
		if(dec1->time==0){
				rnd=rndu();
				dec1->time=rnd*t;
		}
		if(dec2->time==0){
				rnd=rndu();
				dec2->time=rnd*t;
		}


	}


	
	first=NodePop(first);/*memory allocation NB before removing dec1 & dec2*/
	first->type=0;/* type CA */
	first->daughters[0]=dec1;
	first->daughters[1]=dec2;
	first->time=t;
	first->sampled=0;					
	first->Period=0;						/*record period*/
	first->deme=deme;
	if(!ClassicTV)
		SaveNodeinPeriodsArray(first, t);

	
	dec1->previous->next=dec1->next;/* loop maintenance */
	dec1->next->previous=dec1->previous;
	dec2->previous->next=dec2->next;
	dec2->next->previous=dec2->previous;
	
	p=first->ancestral;
	q=dec1->ancestral;
	r=dec2->ancestral;
	for(i=0;i<numBases;i++){
		*p= (*q) | (*r);/* ORs the ancestral states*/
		p++;
		q++;
		r++;
	}
}
/*------------------------------*/
void NetEvolvePlugin::Migration(int deme, int numDemes)
{
double rnd;
int i, migrant, recipDeme;
Node *P;

	do{				/* pick migrant */
		rnd=rndu();
		migrant=((int) (rnd*(ki[deme])) );
	}while(migrant==ki[deme]);
	do{				/* pick recipient deme */
		rnd=rndu();
		recipDeme=( (int) (rnd*numDemes) );
	}while(recipDeme==deme);
	
	P=first;
	i=0;
	while(i<migrant || P->deme!=deme){	/*pick migrant from loop*/
		if(P->deme==deme)
			i++;
		P=P->next;
	}
	P->deme=recipDeme;
	ki[deme]--;
	ki[recipDeme]++;
}

int NetEvolvePlugin::ConstRoutine(Regime *pp)
{
double t, rnd;
int i;
Node *nptr;

	nptr=first;
	for(i=0;i<K;i++){	/* make all genes in same deme */
		nptr->deme=0;
		nptr=nptr->next;
	}
	ki[0]=K;
	t=0.0;
	do{
		if(!NoClock && intervalDis>0.0 && pCounter>=SSizeperP-1 && (t+globTime)<prevTime+intervalDis)
			t=prevTime+intervalDis;
		t+=GenerateTimeEnd(pp);
		if( t > (pp->t) && pp->t > 0.0){
			globTime+=(pp->t);
			return 1;
		}
	//fprintf(stderr, "probRE=%lf K=%d\n",probRE, K);
		rnd=rndu();
		if(rnd<probRE){	/*RE event*/
			Recombine(t+globTime, 0);
			K++;
			ki[0]++;
			noRE++;
		}
		else{			/*CA event*/
			Coalesce(t+globTime, 0);
			K--;
			ki[0]--;
			noCA++;
		}
	}while(K>1);
	/* Detective Search for Recomb rate Bug */
	//	pp->r=0.0;

	globTime+=t;
	return 0;
}


double NetEvolvePlugin::GenerateTimeEnd(Regime *pp) /* for constant pop size and no subdivision */
{
double zz, rnd, t, GR;
long gi;

	do {
		rnd=rndu();
		printf("%lf\n", rnd);
	}
	while(rnd==0.0);
	
	if(pp->r!=0.0){
		gi=CalcGi(0);
		GR=((double) gi) * factr * (pp->N) * genTimeInvVar * (pp->r);
	}else
		GR=0.0;gi=0;
	zz=(double) ( GR + (K*(K-1)) ); /* times are exponentially distributed                */
	probRE=GR/zz;					/* with parameter genTimeInvVar*factr*N*r + (K*(K-1)) */
	t=((-log(rnd))/zz) * factr * genTimeInvVar * (pp->N);
	return t;
}



int NetEvolvePlugin::ConstSubRoutine(Regime *pp)
{
double t, rnd;

	DistributeGenes(pp);
	t=0.0;
	do{
		t+=GenTimeEndSub(pp);
		if( t > (pp->t) && pp->t > 0.0){
			globTime+=(pp->t);
			return 1;
		}
		rnd=rndu();
		if(rnd<probRE){ 				/* RE event */
			Recombine(t+globTime, currDeme);
			K++;
			ki[currDeme]++;
			noRE++;
		}else{			
			if(rnd<(probRE+probMI)){	/* MI event */
				Migration(currDeme, pp->d);
				noMI++;
			}else{						/* CA event */
				Coalesce(t+globTime, currDeme);
				K--;
				ki[currDeme]--;
				noCA++;
			}
		}
	}while(K>1);
	globTime+=t;
	return 0;
}


/*void NetEvolvePlugin::DistributeGenes(Regime *pp)
{
int i, j;
Node *nptr;

	for(i=0;i<(pp->d);i++)
		ki[i]=0;
	nptr=first;
	for(i=0;i<K;i++){
		do
			j=((int) (rndu()*(pp->d)));
		while(j==(pp->d));
		ki[j]++;
		nptr->deme=j;
		nptr=nptr->next;
	}
}*/

static double NetEvolvePlugin::GenTimeEndSub(Regime *pp)
{
int shortest;
long gi;
double GR, MK, zz, t1, t2, rnd;

	t1=BIG_NUMBER;
	for(currDeme=0;currDeme<(pp->d);currDeme++){
		if(ki[currDeme]!=0){
			do
				rnd=rndu();
			while(rnd==0.0);
			if(pp->r!=0.0){
				gi=CalcGi(currDeme);
				GR=( ((double) gi) * factr * genTimeInvVar * (pp->N) * (pp->r) );
			}else
				GR=0.0;gi=0;
			MK=factr * genTimeInvVar * (pp->N) * (pp->m) * (ki[currDeme]);
			zz=GR+MK+( (ki[currDeme])*(ki[currDeme]-1) );
			probRE=GR/zz;
			probMI=MK/zz;
			t2=(-log(rnd))/zz;
			if(t2<t1){
				t1=t2;
				shortest=currDeme;
			}
		}
	}
	currDeme=shortest;
	if(pp->r!=0.0){
		gi=CalcGi(currDeme);
		GR=( ((double) gi) * factr * genTimeInvVar * (pp->N) * (pp->r) );
	}else
		GR=0.0;gi=0;
	MK=factr * genTimeInvVar * (pp->N) * (pp->m) * (ki[currDeme]);
	zz=GR+MK+( (ki[currDeme])*(ki[currDeme]-1) );
	probRE=GR/zz;
	probMI=MK/zz;
	return (t1 * factr * genTimeInvVar * (pp->N) );
}
void NetEvolvePlugin::ReID(Node *node, int Number)
{
	char T[5], NewID[8];
	if (node->Period<10)
		sprintf(T, "00%d",node->Period);
	else if(node->Period<100)
		sprintf(T, "0%d",node->Period);
	else
		sprintf(T, "%d",node->Period);


	sprintf(NewID,"%s.%d",T, Number);
	strcpy(node->ID,NewID);

}

void NetEvolvePlugin::FixInternalNodeClockTime(Node *node, double time)
{/* For internal nodes when clock is used*/
	double childTime;
	int j, period;

	if(node->daughters[0]!=NULL){
		if(node->daughters[0]->time>time ){
			j=0;
			period=LastPeriod-node->daughters[0]->Period;
			while( NodesPerPeriod[period][j]->daughters[0]!=NULL )
				j++;
			childTime=NodesPerPeriod[period][j]->time; /* time of leaves from daughter period*/
			childTime=childTime+rndu()*(time-childTime);
			FixInternalNodeClockTime(node->daughters[0],childTime);
		}
	}
	if(node->daughters[1]!=NULL && node->type!=3 && node->type!=1){
		if(node->daughters[1]->time>time){
			j=0;
			period=LastPeriod-node->daughters[1]->Period;
			while( NodesPerPeriod[period][j]->daughters[0]!=NULL )
				j++;
			childTime=NodesPerPeriod[period][j]->time; /* time of leaves from daughter period*/
			childTime=childTime+rndu()*(time-childTime);
			FixInternalNodeClockTime(node->daughters[1],childTime);
		}
	}
	node->time=time;



}

void NetEvolvePlugin::FixClockPeriods()
{
int i,j;
int period;
	for(i=0;i<LastPeriod;i++){
			j=0;
			while( NodesPerPeriod[i][j]!=NULL){
				period=(int)(((NodesPerPeriod[LastPeriod][0]->time+200)-NodesPerPeriod[i][j]->time)/100);
				if(NodesPerPeriod[i][j]->type!=2){
					if(period>=NodesPerPeriod[i][j]->daughters[0]->Period)
						period=NodesPerPeriod[i][j]->daughters[0]->Period-1;
					if(NodesPerPeriod[i][j]->type==0){
						if(period>=NodesPerPeriod[i][j]->daughters[1]->Period)
							period=NodesPerPeriod[i][j]->daughters[1]->Period-1;
					}
				}
				NodesPerPeriod[i][j]->Period=period;
				j++;
			}
	}
	NodesPerPeriod[LastPeriod][0]->Period=0;
}

void NetEvolvePlugin::PickPeriods2()
{
	int  	 Period, i,j, iNodes,Isolated;
	double iNodeProb2,prev_rndPick,rndPick,rnd ;

		/* Go Period per period */
		NodesPerPeriod[LastPeriod][0]->Period=0;
		for(i=0;i<LastPeriod;i++){
			Period=LastPeriod-i;
			j=0;iNodes=0;
			while( NodesPerPeriod[i][j]!=NULL && j<sampleSize){
				if(NodesPerPeriod[i][j]->daughters[0]!=NULL)
					iNodes++;
				j++;
			}	
			j=0;
			while(NodesPerPeriod[i][j]!=NULL && j<iNodes+SSizeperP){
				if(iNodeProb>1){
					NodesPerPeriod[i][j]->sampled=1;
					numSampled++;
				}
				else{
					NodesPerPeriod[i][j]->sampled=0;
					if((iNodeProb==0 || iNodes==0) && NodesPerPeriod[i][j]->daughters[0]==NULL){
							NodesPerPeriod[i][j]->sampled=1;/* Sample from leaves only*/
							numSampled++;
					}
				}
				NodesPerPeriod[i][j]->Period=Period;
				j++;		
			}
			if(iNodeProb<=1){
				if(iNodes==0||iNodeProb==0)
					iNodeProb2=0;
				else
					iNodeProb2=(iNodes*iNodeProb)/j;
					

				if(iNodeProb2>0){/* random sampling or strategic sampling*/
					j=0;
					Isolated=0;
					prev_rndPick=0;
					if(i<PeriodsWanted){/* Sample */
						while(Isolated<SSizeperP){
								rnd=rndu();
								j=0;
								rndPick=0;
								while( j<iNodes+SSizeperP){
									if(j<iNodes)
										rndPick+=iNodeProb2/iNodes;
									else
										rndPick+=(1-iNodeProb2)/SSizeperP;
									if(prev_rndPick<rnd && rnd<rndPick ){/* Sample and break*/
										if(NodesPerPeriod[i][j]->sampled==1)
											break;
										if(j<iNodes && !NoClock) {/* internal node: rare event*/
											FixInternalNodeClockTime(NodesPerPeriod[i][j],NodesPerPeriod[i][iNodes]->time);
											if (NodesPerPeriod[i][j]->type==1)
												NodesPerPeriod[i][j]->daughters[1]->time=NodesPerPeriod[i][j]->time;
										}
										NodesPerPeriod[i][j]->sampled=1;
										Isolated++;
										numSampled++;
										break;
									}
									prev_rndPick=rndPick;
									j++;
								}
						}
					}
				}
			}
		}
		if(!NoClock)
			FixClockPeriods();

}

void NetEvolvePlugin::PicknAssignID(Node *node,  int *nodeNo, int BKP, int recParentBKP)
{
	double r_s,threshold,effectiveQty;
	int Period;
	
	if(StartAt>0){/* Heavy Sampling from large tree */
		Period=node->Period;
		if(Period==LastPeriod)
			effectiveQty=Interval[Period].Qty-Interval[Period].LookedAt;
		else
			effectiveQty=Interval[Period].Qty-Interval[Period].Time0-Interval[Period].LookedAt;
		if(Period==LastPeriod || node->time!=0.0)	
			Interval[Period].LookedAt++;
	}


	(*nodeNo)++;
	ReID(node,*nodeNo);

	if(StartAt>0){
		if(Period!=LastPeriod && NoClock==1 && node->time==0.0)
				node->sampled=0;
		else if(Interval[Period].Isolated>-1 && Interval[Period].Isolated <SSizeperP){
				r_s=rndu();
					//if(recParentBKP>-1 && (recParentBKP<numBases/8||recParentBKP>(numBases*5)/8))
					//	threshold=0.75;/* For seqs whose children will recombine with breakpoint at ends */
					//else{
						//if(BKP==0)
							threshold=(double)(SSizeperP-Interval[Period].Isolated)/effectiveQty;/* Sample x per time period */
						//else {
						//	if(BKP>numBases/8 && BKP<(numBases*5)/8)
						//		threshold=0.75; /* Recomb child with BKP in the middle */
						//	else
						//		threshold=1; /* BKP at ends */
						//}
					//}
					if ( r_s <threshold  ) { /* sample this sequence */
						node->sampled=1;
						Interval[Period].Isolated++;
						numSampled++;
					}
					else
						node->sampled=0;
		}
		else
			node->sampled=0;
	}
	


}




void NetEvolvePlugin::EnlargeIntervalDimension(int Period)
{
	int i;
	NodeInterval *tmp;
	/* Check if Interval array large enough */
	if (Period+1>=IntervalPeriods){
	  if ((tmp = realloc(Interval, sizeof(NodeInterval) * (IntervalPeriods + 10))) == NULL) {
		fprintf(stderr, "ERROR: realloc failed");
		exit(0);
	  }
	  else{ /* Initialize new array positions*/
		Interval = tmp;
		for(i=IntervalPeriods;i<IntervalPeriods+10;i++){
			Interval[i].Qty=0;
			Interval[i].Time0=0;
			Interval[i].LookedAt=0;
		}
		IntervalPeriods=IntervalPeriods+10;

	  }
	}

}



void NetEvolvePlugin::StorePeriods1b(Node *node)
{
	int Period;
		if(node->L_father!=NULL){
			if(node->L_father->type==1){
				if (node->L_father->Period>node->R_father->Period)
					Period=node->L_father->Period+1;
				else
					Period=node->R_father->Period+1;
			}
			else
					Period=node->L_father->Period+1;
		} else
			Period=0;
		if(LastPeriod<Period) LastPeriod=Period;

		Interval[Period].Qty++;
		if(node->time==0.0)
			Interval[Period].Time0++;

		
		node->Period=Period;
		if(node->L_father!=NULL){
			if (node->L_father->type==1){ /* Recombinant child */
				Interval[Period].Rec=1; /* This interval contains Recombinants */
				TotalRec++;
			}
		}
		
		EnlargeIntervalDimension(Period);

		if(node->type==1){ /* Is recomb father?*/
			if(node->daughters[0]->Period!=-1){/* First time here */
				if (node->cutBefore!=-1){
					node->daughters[0]->L_father=node;/*Only to check fathers Period */
					node->daughters[0]->R_father=node->daughters[1];
				}
				else{
					node->daughters[0]->L_father=node->daughters[1];
					node->daughters[0]->R_father=node;
				}			
				node->daughters[0]->Period=-1;/* stop */
			}
			else
				StorePeriods1b(node->daughters[0]);/* go on */

		}else if(node->type==0){ /* Not recomb */
			node->daughters[0]->L_father=node;/*Only to check fathers Period */
			node->daughters[1]->L_father=node;
			StorePeriods1b(node->daughters[0]);
			StorePeriods1b(node->daughters[1]);
		}
		
	


}

void NetEvolvePlugin::Enqueue(Node *maxNode, Node *Child)
{
	Node *SortSearch, *tmp;


	SortSearch=maxNode;
	while(SortSearch->MaxNext !=NULL){
		if (SortSearch->MaxNext->time > Child->time  )
				SortSearch=SortSearch->MaxNext;
		else
			break;
	}
	tmp=SortSearch->MaxNext; /* Could be NULL */
	SortSearch->MaxNext=Child;
	Child->MaxNext=tmp;

	if(SortSearch==maxNode)
		Head_MaxQueue=Child;


}

void NetEvolvePlugin::StorePeriods1()
{

	Node *maxNode;


	while (Head_MaxQueue != NULL){

		maxNode=Head_MaxQueue;
		Head_MaxQueue=Head_MaxQueue->MaxNext;
		if(maxNode->L_father->Period==LastPeriod)
		{/*If maxNode.parent and maxNode in same interval*/
			Interval[LastPeriod].End=maxNode->time+1;
			LastPeriod++;
			Interval[LastPeriod].Start=maxNode->time;
			Interval[LastPeriod].Rec=0;
		}

		Interval[LastPeriod].Qty++;


		if(MaxQty<Interval[LastPeriod].Qty)
			MaxQty=Interval[LastPeriod].Qty;

		if (maxNode->L_father->type==1){ /* Recombinant child */
			Interval[LastPeriod].Rec=1; /* This interval contains Recombinants */
			TotalRec++;
		}

		maxNode->Period=LastPeriod;

		/*Enqueue maxNode.child1, maxNode.child2(if not recombinant)*/
		if (maxNode->type==1 ){/*recombinant father */
			if ( maxNode->daughters[0]->Period!=-1){/*First time here */
				if (maxNode->cutBefore!=-1){
					maxNode->daughters[0]->L_father=maxNode;/*Only to check fathers Period */
					maxNode->daughters[0]->R_father=maxNode->daughters[1];
				}
				else{
					maxNode->daughters[0]->L_father=maxNode->daughters[1];
					maxNode->daughters[0]->R_father=maxNode;
				}
				maxNode->daughters[0]->Period=-1;
				Enqueue(maxNode,maxNode->daughters[0]); 
			}
				
		} 
		else if(maxNode->type!=2)
			maxNode->daughters[0]->R_father=NULL;
	

		
		if (maxNode->type==0){ /*Common Ancestor */
			maxNode->daughters[0]->L_father=maxNode;/*Only to check fathers Period */
			maxNode->daughters[1]->L_father=maxNode;
			Enqueue(maxNode,maxNode->daughters[0]); 
			Enqueue(maxNode,maxNode->daughters[1]);
		}

		EnlargeIntervalDimension(LastPeriod);

	}/* loop while (Head_MaxQueue != NULL)*/
}

void NetEvolvePlugin::FixSlimAncestors(Node *node)
{

	if(node->type==3||node->type==1){/* Recomb  Parent*/
		if(node->jumpleft!=NULL){
			if(node->cutBefore!=-1){/* Right parent*/
				node->jumpleft->L_father=node;
				FixSlimAncestors(node->jumpleft);
			}
			else
				node->jumpleft->R_father=node;
		}
	}
	else if(node->type==0){/* CA */
		if(node->jumpleft!=NULL){
			node->jumpleft->L_father=node;
			FixSlimAncestors(node->jumpleft);
		}
		if(node->jumpright!=NULL){
			node->jumpright->L_father=node;
			FixSlimAncestors(node->jumpright);
		}
	}
}

Node * NetEvolvePlugin::BuildSlimTree(Node *node)
{
	int i;

	node->jumpleft=NULL;
	node->jumpright=NULL;
	if(node->type==3||node->type==1){/* left=3 right=1 recomb parent -> traverse if other path not visited yet*/
		if(node->daughters[1]->jumpleft==node->daughters[1])/* Self-loop to check if first encounter*/
			node->jumpleft=BuildSlimTree(node->daughters[0]);
		else
			node->jumpleft=node->daughters[1]->jumpleft;
	}
	else if (node->type==0) {/*not tip*/
		node->jumpleft=BuildSlimTree(node->daughters[0]);
		node->jumpright=BuildSlimTree(node->daughters[1]);
		if (node->jumpleft==NULL && node->jumpright!=NULL){/* For tree printing */
				node->jumpleft=node->jumpright;
				node->jumpright=NULL;
		}
		if(node->jumpright==NULL && node->jumpleft!=NULL){
			if(node->jumpleft->sampled==0 && (node->jumpleft->type!=3 && node->jumpleft->type!=1) ){
				/* Remove left non-sampled node from path */
				if(!ClassicTV){
					for(i=0;i<slimNodes;i++){
						if(	OUTStorage[i]==node->jumpleft)
							break;
					}
					if(i<slimNodes-1){
						for(i;i<slimNodes-1;i++)
							OUTStorage[i]=OUTStorage[i+1];

					}
					slimNodes--;
				}
				node->jumpright=node->jumpleft->jumpright;
				node->jumpleft=node->jumpleft->jumpleft;

			}
		}

		
	}

	if(!node->sampled){
		if(node->type==2)/* tip can't check for jumpleft !=NULL */
			return(NULL);
		else if(node->type==1 || node->type==3){/* recomb node can't check jumpright!=NULL */
			if( node->jumpleft==NULL)
				return(NULL);
		}
		else if(node->jumpleft==NULL || node->jumpright==NULL ){
			if(node->jumpleft!=NULL)
				return(node->jumpleft);
			else if(node->jumpright!=NULL)
				return(node->jumpright);
			else
				return(NULL);
		}
	}

	if(!ClassicTV){/* This was used for large tree sampling */
		OUTStorage[slimNodes++]=node;
		ReID(node,slimIDs++);
	}
	else if (node->Period==0)
		ReID(node,slimIDs++);

	return(node);

}

void NetEvolvePlugin::PickPeriods()
{
	int j, SearchOn;
	div_t modzero;


	PeriodsChosen=0;
	SearchOn=0;
	Pick=((LastPeriod-StartAt)/PeriodsWanted);/* Sample every Pick Periods for a total of PeriodsWanted sampling periods */
	if(Pick==0)Pick=1;
	for (j=0; j<IntervalPeriods; j++) {
		if(j>StartAt && PeriodsChosen<PeriodsWanted){
			modzero = div(j,Pick);
			if(modzero.rem==0 ){
				if (Interval[j].Qty-Interval[j].Time0<SSizeperP){
					SearchOn++;
					Interval[j].Isolated =-1;
				}
				else{
					if(SearchOn>0)
						SearchOn--;
					Interval[j].Isolated=0;
					PeriodsChosen++;
				}

			}
			else if(SearchOn>0 && ((Interval[j].Qty-Interval[j].Time0)>=SSizeperP)){
				if(SearchOn>0)
					SearchOn--;
				Interval[j].Isolated=0;
				PeriodsChosen++;
			}
			else
				Interval[j].Isolated=-1;
		}
		else
			Interval[j].Isolated=-1;
	}

}

void NetEvolvePlugin::PrintNode(FILE *treeFile, Node *node, double maxTime)
{
	double len2Parent;



	if (node->jumpleft==NULL || node->cutBefore>-1 ){/* leaf or left rec parent*/
		/*rec left parent becomes leaf when encountered */
		if(node->sampled==0)
			fputc('~', treeFile);				
		fprintf(treeFile, "%s", node->ID); 
	}
	else  {
		if( node->sampled==1){/* internal sampled Node*/
			/* treat as leaf with 0-length branch*/
			fputc('(', treeFile);	
			fprintf(treeFile, "%s:0,", node->ID);
		}
		if(node->jumpright!=NULL ||node->type==1 || node->type==3 ){

			fputc('(', treeFile);				
			PrintNode(treeFile,  node->jumpleft, maxTime);
			fputc(',', treeFile);		/* Recombination? */	
			if(node->type==1 || node->type==3) {
				if(node->daughters[1]->sampled==0)
					fputc('~', treeFile);
				/* If here: right parent: print left parent */
				fprintf(treeFile, "%s#%d:0", node->daughters[1]->ID,node->daughters[1]->cutBefore);
			}
			else
				PrintNode(treeFile,  node->jumpright, maxTime);
			fputc(')', treeFile);	
			if( node->sampled==1){/* internal sampled Node*/
				/* treat as leaf with 0-length branch*/
				fprintf(treeFile, ":0");
				fputc(')', treeFile);	

			}
		}
		else{
			PrintNode(treeFile,  node->jumpleft, maxTime);
			//if(node->jumpleft->jumpleft==NULL)/* If internal sampled parent of leaf*/
			fputc(')', treeFile);	
		}

	}
	
	/* Calculate tree lengths */
	len2Parent=(node->L_father->time-node->time)/maxTime;
	if (len2Parent<0)
		fprintf (stderr,"negative branch length:%s\n",node->ID);
	if(len2Parent==0)
		fprintf(treeFile, ":0");
	else
		fprintf(treeFile, ":%lf", len2Parent);
	


}

void NetEvolvePlugin::PrintNetworkOrTree(FILE *fv, Node *node, double maxTime)
{
	if(node->jumpleft==NULL||node->jumpright==NULL){
		if(node->jumpleft!=NULL)
			PrintNode(fv, node->jumpleft, maxTime);
		else if(node->jumpright!=NULL)
			PrintNode(fv, node->jumpright, maxTime);

	}
	else{
		fputc('(', fv);	
		PrintNode(fv, node->jumpleft, maxTime);
		fputc(',', fv);			
		PrintNode(fv, node->jumpright, maxTime);
		fprintf(fv, ");\n");
	}
}


void NetEvolvePlugin::SeqEvolveSSTV (char * pbase, char *ppath)
{
	int i, j, l;
	char *P;
	char  Anc[20],name[100],space[2]; 
	FILE *ance_f, *out_phy, *treeFile ;
	int n; /* Actual Number of Sampled Seqs */
	Node *PNode;

	remove("TestOutput.txt");
	space[0]=' ';
	space[1]='\0';
	numSampled=0;

	tipNo=0;
	RandomSequence(first->sequence);
	SetCategories();
/* MY Code starts Here: */
	std::string base1 = std::string(PluginManager::prefix())+"/"+std::string(pbase);
	std::string base2 = std::string(PluginManager::prefix())+"/out";
	if(outputFile)
		sprintf(name,"%s%d.txt",base1.c_str(),repNo);
	else
		sprintf(name,"%s%d.txt",base2.c_str(), repNo);
	if(repNo==1){
		fprintf(stdout, "%s %.2f - Serial Coalescent Simulation and Sequence Evolution\n",PROGRAM_NAME ,VERSION_NO);
		fprintf(stdout, "--------------------------------------------------------------------------\n\n");
	}
	if(outputFile)
		fprintf(stdout, "Data set replicate #%d is written to %s, tree to %s%d.tre \n",repNo,name,pbase,repNo);
	else
		fprintf(stdout, "Data set replicate #%d is written to %s, tree to tv%d.tre in program's folder.\n",repNo,name,repNo);
	if (freopen (name,"w",stdout)==NULL){
		fprintf(stderr, "Error opening file: '%s'\n",name);
		exit(0);
	}

	std::string base3 = std::string(PluginManager::prefix())+"/"+std::string(ppath);
	std::string base4 = std::string(PluginManager::prefix())+"/A";
	if(outputFile)
		sprintf(name,"%sA%d.txt",base3.c_str(),repNo);
	else
		sprintf(name,"%s%d.txt",base4.c_str(), repNo);
	if ( (ance_f=fopen(name, "w"))==NULL ) {
		fprintf(stderr, "Error opening file: '%s'\n",name);
		exit(0);
	}


	if(outputFile)
		sprintf(name,"%s%d.phy",base1.c_str(),repNo);
	else
		sprintf(name,"%s%d.phy",base2.c_str(), repNo);
	if ( (out_phy=fopen(name, "w"))==NULL ) {
		fprintf(stderr, "Error opening file: '%s'\n",name);
		exit(0);
	}
	std::string base5 = std::string(PluginManager::prefix())+"tv";
	if(outputFile)
		sprintf(name,"%s%d.tre",base1.c_str(),repNo);
	else
		sprintf(name,"%s%d.tre",base5.c_str(), repNo);
	if ( (treeFile=fopen(name, "w"))==NULL ) {
		fprintf(stderr, "Error opening tree file: '%s'\n",name);
		exit(0);
	}


	MaxQty=0;
	TotalRec=0;
	LastPeriod=0;
	first->L_father=NULL;
	first->R_father=NULL;
	/*** StorePeriods code ***/
	if(StartAt>0){/* Heavy Sampling from large tree */
		IntervalPeriods=10;

		MALLOC(Interval, sizeof(NodeInterval) * sampleSize);
		for(i=0;i<sampleSize;i++){
			Interval[i].Qty=0;
			Interval[i].Time0=0;
			Interval[i].LookedAt=0;
			Interval[i].Isolated=0;
		}
		if(NoClock==0){
			Interval[LastPeriod].Start=first->time;
			Interval[LastPeriod].Qty=1;
			
			first->Period=LastPeriod;
			first->daughters[0]->L_father=first;
			first->daughters[1]->L_father=first;
			if (first->daughters[0]->time >first->daughters[1]->time){
				Head_MaxQueue=first->daughters[0];
				Head_MaxQueue->MaxNext=first->daughters[1];
				first->daughters[1]->MaxNext=NULL;
			}else{
				Head_MaxQueue=first->daughters[1];
				Head_MaxQueue->MaxNext=first->daughters[0];
				first->daughters[0]->MaxNext=NULL;
			
			}
			StorePeriods1();
			Interval[LastPeriod].End=0;
		}else/* Store Periods without clock assumption */
			StorePeriods1b(first);
		PickPeriods();
	}
	else{
		LastPeriod=first->Period;
		PickPeriods2();
	}

	/*** End StorePeriods code ***/



	sprintf(Anc,"000.0");
	strcpy(first->ID,Anc);



	/* Mutate & pick isolated sequences */
	l=1;
	n=0;
	Mutate(first,&n);

	Stack(first);

	MALLOC(OUTStorage, sizeof(Node *) * (sampleSize * range));/* nodes in slimTree - should be prop.to rec rate */
	slimNodes=0;slimIDs=0;
	PNode=BuildSlimTree(first);
	if(PNode->jumpleft==NULL && PNode->jumpright==NULL)
		printf("problem with jumps!\n");
	FixSlimAncestors(PNode);
	PrintNetworkOrTree(treeFile,first,first->time);
	/* Use OutNode to print: first heapsort it */
	slimNodes--;
	for(i=0;i<slimNodes;i++)
		heapify_Nodelist(OUTStorage,i);
	heapsort_Nodelist(OUTStorage,slimNodes-1);

	fprintf(out_phy,"%d %d I O\nO 1 \n",numSampled+1,numBases);
	fprintf(out_phy, "%s     ", first->ID);/* Root or First internal node  */
	P=first->sequence;
	for (j=0; j<numBases; j++) {
		fputc(nucleotides[*P], out_phy);
		P++;
	}
	fputc('\n', out_phy);

	/* Now Print */
	fprintf(ance_f, "%s;%lf \n", first->ID,first->time,Pick);
	for (i=0;i<slimNodes;i++){

		if(OUTStorage[i]->L_father->type==3 || OUTStorage[i]->L_father->type==1)/* recomb has two ancestors*/
			fprintf(ance_f, "%s;%s|%d|%s;%lf;%d \n", OUTStorage[i]->ID,OUTStorage[i]->L_father->ID, OUTStorage[i]->L_father->cutBefore,OUTStorage[i]->R_father->ID, OUTStorage[i]->time, OUTStorage[i]->sampled);
		else
			fprintf(ance_f, "%s;%s;%lf;%d \n", OUTStorage[i]->ID,OUTStorage[i]->L_father->ID,OUTStorage[i]->time,OUTStorage[i]->sampled);

		/* Only the sampled */
		if(OUTStorage[i]->sampled==1){
			strcpy(name,OUTStorage[i]->ID);
			for(j=strlen(name)+1;j<11;j++)
				strcat(name,space);
			fprintf(out_phy,"%s",name);
			fprintf(stdout, ">%s\n",OUTStorage[i]->ID);
			P = OUTStorage[i]->sequence;
			for (l=0; l<numBases; l++) {
				fputc(nucleotides[*P], stdout);
				fputc(nucleotides[*P], out_phy);
				P++;
			
			}
			fputc('\n', stdout);
			fputc('\n', out_phy);
		}
	}
	

	fclose(ance_f);
	fclose(out_phy);
	fclose(treeFile);
	if(!ClassicTV)
		fclose (stdout);	
	free(OUTStorage);
	if(StartAt>0)
		free(Interval);
}

void NetEvolvePlugin::SeqEvolve (char * pbase, char *ppath)
{
	int n;
	FILE * treeFile;
	Node *PNode;
	char name[100];

	if(outputFile){
		sprintf(name,"%s%d.txt",pbase,repNo);
	
		fprintf(stdout,"Process: 90%%\n");
		if (freopen (name,"w",stdout)==NULL){
			fprintf(stderr, "Error opening file: '%s'\n",name);
			exit(0);
		}
	}
	tipNo=0;
	RandomSequence(first->sequence);
	SetCategories();
	fprintf(stdout, " %d %d\n", sampleSize, numBases);
	n=0;
	Mutate(first,&n);
	Stack(first);

	if(outputFile)
		sprintf(name,"%s%d.tre",pbase,repNo);
	else
		sprintf(name,"tv%d.tre",repNo);
	if ( (treeFile=fopen(name, "w"))==NULL ) {
		fprintf(stderr, "Error opening tree file: '%s'\n",name);
		exit(0);
	}
	slimIDs=tipNo;
	PNode=BuildSlimTree(first);
	if(PNode->jumpleft==NULL && PNode->jumpright==NULL)
		printf("problem with jumps!\n");
	FixSlimAncestors(PNode);
	PrintNetworkOrTree(treeFile,first,first->time);
	fclose(treeFile);

}

void NetEvolvePlugin::Mutate(Node *node, int *nodeNo)
{
int i, Dir, BKP,cutBefore;//, expNo, posn;
//float xm;
char  *p1, *p2, *child;//*point,
//double rnd;


BKP=0;

	node->jumpleft=node; /* to avoid traversing tree twice in BuildSlimTree */
	node->jumpright=node;

	if(node->type==1){						/*i.e. recombination event */
		if(node->daughters[1]->type==1)		/* first encounter */
			node->type=3;
		else{								/* second encounter */
			p1=node->sequence;
			p2=node->daughters[1]->sequence;/* remember this is spare for re and points to other parent */
			child=node->daughters[0]->sequence;
			if(node->cutBefore!=-1){ 
				Dir=1;
				BKP=node->cutBefore;
			}else{
				Dir=0;
				BKP=node->daughters[1]->cutBefore;
			}

			for(i=0;i<numBases;i++){
				if(node->ancestral[i]==1)
					*child=*p1;
				else
					*child=*p2;
				child++;
				p1++;
				p2++;
				
			}
			*child='\0';
			MutateSequence(node->daughters[0]->sequence, (mutRate*((double) node->time - (double) node->daughters[0]->time)));

			//sprintf(Ance1,"%s-R%d%R-%s", node->daughters[1]->ID, BKP, node->ID);
			if(!ClassicTV)
				PicknAssignID(node->daughters[0], nodeNo, BKP,-1);
			if(node->daughters[0]->type==2){
				if(ClassicTV){
					tipNo++;
					WriteSequence(tipNo, node->daughters[0]->sequence,node->daughters[0]);
				//WriteSequence(node->daughters[0], node->daughters[0]->sequence);
				}
			}else
				Mutate(node->daughters[0], nodeNo);

		
			Stack(node->daughters[1]);		/* stack up memory */
			if(node->daughters[0]->type!=3)
				Stack(node->daughters[0]);
		}
	}else{									/* a coalescent event */

		/* we're moving to the daughter Period */


		memcpy(node->daughters[0]->sequence, node->sequence, (numBases+1));
		MutateSequence(node->daughters[0]->sequence, (mutRate*((double) node->time - (double) node->daughters[0]->time)));
		cutBefore=-1;
		if(node->daughters[0]->type==1||node->daughters[0]->type==3 ){
			if(node->daughters[0]->cutBefore>-1 )
				cutBefore=node->daughters[0]->cutBefore;
			else if (node->daughters[0]->daughters[1]->cutBefore>-1)
				cutBefore=node->daughters[0]->daughters[1]->cutBefore;
		}
		if(!ClassicTV)
			PicknAssignID(node->daughters[0],  nodeNo, 0, cutBefore);
		
		if(node->daughters[0]->type==2 ){/* it's a leaf */
			if(ClassicTV){
				tipNo++;
				WriteSequence(tipNo, node->daughters[0]->sequence,node->daughters[0]);
			}
		}else
			Mutate(node->daughters[0],  nodeNo);
		if(node->daughters[0]->type!=3)
			Stack(node->daughters[0]);		/* stack up memory */
		
		memcpy(node->daughters[1]->sequence, node->sequence, (numBases+1));
		MutateSequence(node->daughters[1]->sequence, (mutRate*((double) node->time - (double) node->daughters[1]->time)));
		cutBefore=-1;
		if(node->daughters[1]->type==1 || node->daughters[1]->type==3){
			if(node->daughters[1]->cutBefore>-1 )
				cutBefore=node->daughters[1]->cutBefore;
			else if (node->daughters[1]->daughters[1]->cutBefore>-1)
				cutBefore=node->daughters[1]->daughters[1]->cutBefore;
		}
		if(!ClassicTV)
			PicknAssignID(node->daughters[1], nodeNo, 0,cutBefore);
		if(node->daughters[1]->type==2 ){
			if(ClassicTV){
				tipNo++;
				WriteSequence(tipNo, node->daughters[1]->sequence,node->daughters[1]);
			}

		}else 
			Mutate(node->daughters[1],  nodeNo);

			
		if(node->daughters[1]->type!=3)
			Stack(node->daughters[1]);		/* stack up memory */
	}

}

void NetEvolvePlugin::RateHetMem()
{
	if (rateHetero==GammaRates){
		if( (siteRates=(double*)calloc(numBases, sizeof(double)))==NULL){
			fprintf(stderr, "Out of memory allocating site specific rates (RateHetMem)\n");
			exit(0);
		}
	}else if (rateHetero==DiscreteGammaRates){
		if( (categories=(short*)calloc(numBases, sizeof(short)))==NULL){
			fprintf(stderr, "Out of memory allocating discrete gamma categories (RateHetMem)\n");
			exit(0);
		}
	}
}


void NetEvolvePlugin::SetCategories()
{
	int  i;//, cat,j;
	double sumRates;
	
	if (rateHetero==CodonRates) {
		sumRates=catRate[0]+catRate[1]+catRate[2];
		if (sumRates!=3.0) {
			catRate[0]*=3.0/sumRates;
			catRate[1]*=3.0/sumRates;
			catRate[2]*=3.0/sumRates;
		}
	} else if (rateHetero==GammaRates) {
		for (i=0; i<numBases; i++)
			siteRates[i]=rndgamma(gammaShape) / gammaShape;

	} else if (rateHetero==DiscreteGammaRates) {
		DiscreteGamma(freqRate, catRate, gammaShape, gammaShape, numCats, 0);
		for (i=0; i<numBases; i++)
			categories[i]=(int)(rndu()*numCats);
	}
}


char NetEvolvePlugin::SetBase(double *P)
{
	char j;
	double r;
	
	r=rndu();
	for (j='\0'; r>(*P) && j<'\3'; j++) P++;
	return j;
}


void NetEvolvePlugin::RandomSequence(char *seq)
{
	int i;
	char *P;
	
	P=seq;
	for (i=0; i<numBases; i++) {
		*P=SetBase(addFreq);
		P++;
	}
	*P='\0';
}

void NetEvolvePlugin::MutateSequence(char *seq, double len)
{
	int i, cat;
	double *Q;
	short *R;
	char *P;
	
	P=seq;
	
	switch (rateHetero) {
		case GammaRates:
			Q=siteRates;
			
			for (i=0; i<numBases; i++) {
				SetVector(vector, *P, (*Q)*len);
				*P=SetBase(vector);
				P++;
				Q++;
			}
		break;
		case DiscreteGammaRates:
			for (i=0; i<numCats; i++)
				SetMatrix(matrix[i], catRate[i]*len);
			
			R=categories;
			for (i=0; i<numBases; i++) {
				*P=SetBase(matrix[*R]+(*P<<2));
				P++;
				R++;
			}
		break;
		case CodonRates:
			for (i=0; i<numCats; i++)
				SetMatrix(matrix[i], catRate[i]*len);
			
			for (i=0; i<numBases; i++) {
				cat=i%3;
				*P=SetBase(matrix[cat]+(*P<<2));
				P++;
			}
		break;
		case NoRates:
			SetMatrix(matrix[0], len);
			
			for (i=0; i<numBases; i++) {
				*P=SetBase(matrix[0]+(*P<<2));
				P++;
			}
		break;
	}
}

//void WriteSequence(int sampleNo, char *P)
//void WriteSequence(Node *node, char *P)
//{
//int j;
//	fprintf(stdout, "sample%s\n", node->ID);
//	for (j=0; j<numBases; j++) {
//		fputc(nucleotides[*P], stdout);
//		P++;
//	}
//	fputc('\n', stdout);
//}


void NetEvolvePlugin::WriteSequence(int sampleNo, char *P, Node *node)
{
int j;
	ReID(node,sampleNo);

	fprintf(stdout, "%-8s ", node->ID);
	for (j=0; j<numBases; j++) {
		fputc(nucleotides[*P], stdout);
		P++;
	}
	fputc('\n', stdout);

}
int NetEvolvePlugin::EpiRoutine(Regime *pp)
{
Node *nptr;
int i;
double t, rnd, tiMe;
	fprintf(stderr, "In EpiRoutine\n");

	nptr=first;
	for(i=0;i<K;i++){	/* make all genes in same deme */
		nptr->deme=0;
		nptr=nptr->next;
	}
	ki[0]=K;
	
	recRate=pp->r;
	lambda=pp->e;
	N=pp->N;
	t=0.0;
	do{
		tiMe=GenerateTimeEpi();
		if(t+tiMe  > (pp->t) && pp->t > 0.0 ){
			globTime+=(pp->t);
			fprintf(stderr, "Population size at t = %1f is %f\n", globTime, (N*exp(-lambda*((pp->t)-t))));
			return 1;
		}
		if(!NoClock && intervalDis>0.0 && pCounter>=SSizeperP-1 && (t+globTime+tiMe)<prevTime+intervalDis)
			t=prevTime+intervalDis;
		N=(N*exp(-lambda*tiMe));
		t+=tiMe;
		rnd=rndu();
		if(rnd<probRE){	/*RE event*/
			Recombine(t+globTime, 0);
			K++;
			ki[0]++;
			noRE++;
		}
		else{			/*CA event*/
			Coalesce(t+globTime, 0);
			K--;
			ki[0]--;
			noCA++;
		}
	}while(K>1);
	globTime+=t;
	fprintf(stderr, "Population size at t = %1f is %f\n", globTime, N);
	return 0;
}

double NetEvolvePlugin::GenerateTimeEpi()
{
double t, brackets[2], P;

	do
		Pn=rndu();
	while(Pn==0.0);
	t=approxFunc();
	brackets[0]=0.0;
	brackets[1]=t;
	P=epifunc(brackets[1]);
	while(P>0.0){
		brackets[1]+=brackets[1]*BRACKET;
		P=epifunc(brackets[1]);
	}
	t=zeroinA(brackets[0], brackets[1]);//, epifunc);
	return t;
}

double NetEvolvePlugin::approxFunc()	
{
double t;
long gi;
	if(recRate!=0.0)
		gi=CalcGi(0);
	else
		gi=0;
	t= -(log(Pn)) / (((double) gi * recRate ) + ( ((double)(K*(K-1))) / (factr*genTimeInvVar*N) ) );
	return t;
	
}

double NetEvolvePlugin::epifunc(double t)	
{
double a, b, temp, diff;
long gi;
	if(recRate!=0.0){
		gi=CalcGi(0);
		a=(double) recRate*gi*t;
	}else
		a=0.0;gi=0;
	temp=( ( (double) (K*(K-1)) ) / (lambda*factr*genTimeInvVar*N) );
	b= (exp(lambda*t) - 1.0) * temp;
	diff=(exp(-a-b))-Pn;
	if(a!=0.0)
		probRE=a/(a+b);
	else
		probRE=0.0;
	return diff;
}
int NetEvolvePlugin::EpiSubRoutine(Regime *pp)
{
double t, tiMe, rnd, Ntemp;
//int i;

	DistributeGenes(pp);
	lambda=pp->e;
	N=pp->N;
	migrationRate=pp->m;
	recRate=pp->r;
	t=0.0;
	do{
		tiMe=GenTimeEpiSub(pp);
		if(t+tiMe  > (pp->t) && pp->t > 0.0 ){	/* End of this Regime */
			Ntemp=((pp->N)*exp(-lambda*((pp->t)-t)));
			if(Ntemp<MIN_NE){
				InstantCoalesce(t, pp);
				return 0;
			}
			globTime+=(pp->t);
			return 1;
		}
		Ntemp=(N*exp(-lambda*tiMe));
		if(Ntemp<MIN_NE){
			InstantCoalesce(t, pp);
			return 0;
		}
		N=Ntemp;
		t+=tiMe;
		rnd=rndu();
		if(rnd<probRE){ 				/* RE event */
			Recombine(t+globTime, currDeme);
			K++;
			ki[currDeme]++;
			noRE++;
		}else{			
			if(rnd<(probRE+probMI)){	/* MI event */
				Migration(currDeme, pp->d);
				noMI++;
			}else{						/* CA event */
				Coalesce(t+globTime, currDeme);
				K--;
				ki[currDeme]--;
				noCA++;
			}
		}
	}while(K>1);
	globTime+=t;
	fprintf(stderr, "Population size at t = %1f is %f\n", globTime, N);
	return 0;
			
}

static void NetEvolvePlugin::DistributeGenes(Regime *pp)
{
int i, j;
Node *nptr;

	for(i=0;i<(pp->d);i++)
		ki[i]=0;
	nptr=first;
	for(i=0;i<K;i++){
		do
			j=((int) (rndu()*(pp->d)));
		while(j==(pp->d));
		ki[j]++;
		nptr->deme=j;
		nptr=nptr->next;
	}
}
void NetEvolvePlugin::InstantCoalesce(double t, Regime *pp)
{
double tiMe;
Node *nptr;
int i;

	fprintf(stderr, "Each subpopulation has reached 1.0 before final coalescence,\n");
	fprintf(stderr, "with the number of genes at %d.\n", K);
	fprintf(stderr, "All genes are therefore instantly coalescing\n");
	tiMe= -(log(MIN_NE/N)) / lambda;
	N=MIN_NE;
	t+=tiMe;
	nptr=first;
	for(i=0;i<K;i++){
		nptr->deme=0;
		nptr=nptr->next;
	}
	for(i=0;i<pp->d;i++)
		ki[i]=0;
	ki[0]=K;
	while(K>1){
		Coalesce(t+globTime, 0);
		K--;
		ki[0]--;
		noCA++;
	}
}


double NetEvolvePlugin::GenTimeEpiSub(Regime *pp)
{
double t1, t2, brackets[2], P;
int shortest;
double a, b, m, temp;
long gi;

	t1=HUGE_VAL;
	for(currDeme=0;currDeme<pp->d;currDeme++){
		if(ki[currDeme]!=0){
			do
				Pn=rndu();
			while(Pn==0.0);
			if(migrationRate!=0.0)
				t2=approxFuncSub();/* exact if ki=1 (but inf. if m=0.0) */
			else
				t2=BIG_NUMBER;
			
			if(ki[currDeme]!=1){
				brackets[0]=0.0;
				brackets[1]=t2;
				P=epiSubfunc(brackets[1]);
				while(P>0.0){
					brackets[1]+=brackets[1]*BRACKET;
					P=epiSubfunc(brackets[1]);
				}
				t2=zeroinB(brackets[0], brackets[1]);//, epiSubfunc);
			}
			
			if(t2<t1){
				t1=t2;
				shortest=currDeme;
			}
		}
	}	
	currDeme=shortest;
	if(recRate!=0.0){
		gi=CalcGi(currDeme);
		a=(double) recRate*gi*t1;
	}else
		a=0.0;gi=0;
	
	m=(double) migrationRate*ki[currDeme]*t1;
	if(ki[currDeme]!=1){
		temp=( ( (double) ((ki[currDeme])*(ki[currDeme]-1)) ) / (lambda*factr*genTimeInvVar*N) );
		b= temp*(exp(lambda*t1) - 1.0);
	}else
		b=0.0;
	probRE=a/(a+m+b);
	probMI=m/(a+m+b);
	return t1;
}

double NetEvolvePlugin::approxFuncSub()
{
double t, zz;
long gi;
	
	if(recRate!=0.0)
		gi=CalcGi(currDeme);
	else
		gi=0;
	zz=( ((double) (ki[currDeme]*(ki[currDeme]-1)))/(factr*genTimeInvVar*N));
	t= -(log(Pn)) / (((double) gi * recRate) + (migrationRate*ki[currDeme]) + zz);
	return t;

}


double NetEvolvePlugin::epiSubfunc(double t)
{
double a, b, m, temp, diff;
long gi;
	if(recRate!=0.0){
		gi=CalcGi(currDeme);
		a=(double) recRate*gi*t;
	}else
		a=0.0;gi=0;
	m=(double) migrationRate*ki[currDeme]*t;
	temp=( ( (double) ((ki[currDeme])*(ki[currDeme]-1)) ) / (lambda*factr*genTimeInvVar*N) );
	b= temp*(exp(lambda*t) - 1.0);
	diff=(exp(-a-m-b))-Pn;
	return diff;
}
/* Convert any random array to heap array*/
void NetEvolvePlugin::heapify_Nodelist ( Node **list , int newnode )
{
	int done=0, dad;
	Node *temp;
	double nNodeID, dadID;

	dad = ( newnode - 1 ) / 2;

	while( dad >= 0 && !done)
	{
		nNodeID=atof(list[newnode]->ID);
		dadID=atof(list[dad]->ID);
		if(nNodeID> dadID)
		 {
			   temp = list[newnode];
			   list[newnode]=list[dad];
			   list[dad]=temp;
			   newnode = dad;
			   dad = dad/2;
		 }
		 else
			done = 1;

	}
}
/* Heap sort*/
void NetEvolvePlugin::heapsort_Nodelist ( Node **list, int last )
{
	int i;
	Node *temp;
	while(last > 0 )
	{
		temp = list[0];
		list[0]=list[last];
		list[last]= temp;
		for(i=0;i<last;i++)
			heapify_Nodelist(list,i);
		last--;
	}
}

Node * NetEvolvePlugin::FirstNodePop()
{
Node *child;
int i;

	if(avail==NULL){
		if ( (child=(Node*)malloc(sizeof(Node)))==NULL ) {
			fprintf(stderr, "Out of memory allocating Node (NodePop)\n");
			exit(0);
		}
		if ( (child->ancestral=(short*)malloc((sizeof(short))*numBases))==NULL ) {
			fprintf(stderr, "Out of memory allocating ancestral gamete (IndiPop)\n");
			exit(0);
		}
		if ( (child->sequence=(char*)malloc((sizeof(char))*(numBases+1)))==NULL ) {
			fprintf(stderr, "Out of memory allocating sequence (IndiPop)\n");
			exit(0);
		}
	}else{
		child=avail;
		avail=avail->next;
	}
	for(i=0;i<numBases;i++)
		child->ancestral[i]=1;
	child->next=child;
	child->previous=child;
	child->cutBefore=-1;
	return child;
}
/*---------------------------------------------------------------------------------------*/
Node * NetEvolvePlugin::NodePop(Node *prev)
{
Node *child;
int i;
	if(avail==NULL){
		if ( (child=(Node*)malloc(sizeof(Node)))==NULL ) {
			fprintf(stderr, "Out of memory allocating Node (NodePop)\n");
			exit(0);
		}
		if ( (child->ancestral=(short*)malloc((sizeof(short))*numBases))==NULL ) {
			fprintf(stderr, "Out of memory allocating ancestral gamete (IndiPop)\n");
			exit(0);
		}
		if ( (child->sequence=(char*)malloc((sizeof(char))*(numBases+1)))==NULL ) {
			fprintf(stderr, "Out of memory allocating sequence (IndiPop)\n");
			exit(0);
		}
	}else{
		child=avail;
		avail=avail->next;
	}

	for(i=0;i<numBases;i++)
		child->ancestral[i]=1;
		
	child->next=prev;
	prev->previous->next=child;
	child->previous=prev->previous;
	prev->previous=child;
	child->cutBefore=-1;
	return child;
}
/*----------------------------------------------------------------------------*/
void NetEvolvePlugin::Stack(Node *going)
{
	going->next=avail;
	avail=going;
}
/*************************************/
void NetEvolvePlugin::SetModel(int theModel)
{	
	int i,j,k;
	double xi, xv, F84temp1, F84temp2,  sumFreq;
	double freqAG, freqCT, freqA2, freqC2, freqG2, freqT2;
	
	model=theModel;

	sumFreq=freq[A]+freq[C]+freq[G]+freq[T];
	if (sumFreq!=1.0) {
		freq[A]*=1.0/sumFreq;
		freq[C]*=1.0/sumFreq;
		freq[G]*=1.0/sumFreq;
		freq[T]*=1.0/sumFreq;
	}
	freqA=freq[A];
	freqC=freq[C];
	freqG=freq[G];
	freqT=freq[T];
	addFreq[A]=freqA;
	addFreq[C]=addFreq[A]+freqC;
	addFreq[G]=addFreq[C]+freqG;
	addFreq[T]=addFreq[G]+freqT;
	
	if (model==F84 || model==HKY) {
		freqR=freqA+freqG;
		freqY=freqC+freqT;
		freqAG=freqA*freqG;
		freqCT=freqC*freqT;
		
		tab1A=freqA*((1/freqR)-1);
		tab2A=(freqR-freqA)/freqR;
		tab3A=freqA/freqR;
		tab1C=freqC*((1/freqY)-1);
		tab2C=(freqY-freqC)/freqY;
		tab3C=freqC/freqY;
		tab1G=freqG*((1/freqR)-1);
		tab2G=(freqR-freqG)/freqR;
		tab3G=freqG/freqR;
		tab1T=freqT*((1/freqY)-1);
		tab2T=(freqY-freqT)/freqY;
		tab3T=freqT/freqY;
	}
		
	switch (model) {
		case F84:
			freqAG=freq[A]*freq[G];
			freqCT=freq[C]*freq[T];
			freqA2=freq[A]*freq[A];
			freqC2=freq[C]*freq[C];
			freqG2=freq[G]*freq[G];
			freqT2=freq[T]*freq[T];
			F84temp1=freqA2+freqC2+freqG2+freqT2;
			F84temp2=((freqA2/freqR)+(freqC2/freqY)+(freqG2/freqR)+(freqT2/freqY));

			xi=freqR*freqY*(freqR*freqY*tstv-freqAG-freqCT);	
			xv=(freqCT*freqR)+(freqAG*freqY);
			kappa=xi/xv;
				
			mu=-1.0/((1-F84temp1)+(kappa*(1-F84temp2)));
			mu_kappa_1=mu*(kappa+1);
		break;
		case HKY:
			kappa=(tstv*freqR*freqY)/(freqAG+freqCT);
			beta=-1.0/(2*(freqR*freqY + kappa*(freqAG+freqCT)));
			
			beta_A_R=beta*(1.0+freqR*(kappa-1));
			beta_A_Y=beta*(1.0+freqY*(kappa-1));
		break;
		case REV:
			k=0;
			for (i=0; i<3; i++) 
				for (j=i+1; j<4; j++)
		      		if (i*4+j!=11)
						Qij[i*4+j]=Qij[j*4+i]=Rmat[k++];
			Qij[3*4+2]=Qij[2*4+3]=1.0;
			
			for (i=0; i<4; i++) 
				for (j=0; j<4; j++) 
					Qij[i*4+j] *= freq[j];
					
			mr=0;		
			for (i=0; i<4; i++) {
				Qij[i*4+i]=0; 
				Qij[i*4+i]=-(Qij[i*4]+Qij[i*4+1]+Qij[i*4+2]+Qij[i*4+3]); 
				mr-=freq[i]*Qij[i*4+i];
			}
			
			EigenREV(Root, Cijk);
			
		/* calculate mean ts/tv ratio */
			mr=2*(freq[3]*Qij[3*4+1]+freq[0]*Qij[0*4+2]);
			tstv=mr/(1-mr);
		break;
	}
}


#define PIJ_SAME_A freqA+(tab1A*aa)+(tab2A*bbR)
#define PIJ_TS_A freqA+(tab1A*aa)-(tab3A*bbR)
#define PIJ_TV_A freqA*(1-aa)

#define PIJ_SAME_C freqC+(tab1C*aa)+(tab2C*bbY)
#define PIJ_TS_C freqC+(tab1C*aa)-(tab3C*bbY)
#define PIJ_TV_C freqC*(1-aa)
	
#define PIJ_SAME_G freqG+(tab1G*aa)+(tab2G*bbR)
#define PIJ_TS_G freqG+(tab1G*aa)-(tab3G*bbR)
#define PIJ_TV_G freqG*(1-aa)
	
#define PIJ_SAME_T freqT+(tab1T*aa)+(tab2T*bbY)
#define PIJ_TS_T freqT+(tab1T*aa)-(tab3T*bbY)
#define PIJ_TV_T freqT*(1-aa)	
	

void NetEvolvePlugin::SetMatrix(double *matrix, double len)
{	
	double aa, bbR, bbY;
	int i,j,k;
	double expt[4];
	double *P;
	
	if (model==REV) {
/* P(t)ij = SUM Cijk * exp{Root*t}
*/
		P=matrix;
		if (len<1e-6) { 
			for (i=0; i<4; i++) 
				for (j=0; j<4; j++) {
					if (i==j)
						*P=1.0;
					else 	
						*P=0.0;
					P++;
				}
			return; 
		}
		
		for (k=1; k<4; k++) 
			expt[k]=exp(len*Root[k]);
		for (i=0; i<4; i++) 
			for (j=0; j<4; j++) {
				(*P)=Cijk[i*4*4+j*4+0];
				for (k=1; k<4; k++)
					(*P)+=Cijk[i*4*4+j*4+k]*expt[k];
				P++;
			}
	} else {
		switch (model) {
			case F84:
				aa=exp(mu*len);
				bbR=bbY=exp(mu_kappa_1*len);
			break;
			case HKY:
				aa=exp(beta*len);
				bbR=exp(beta_A_R*len);
				bbY=exp(beta_A_Y*len);
			break;
		}
		matrix[0]=PIJ_SAME_A;
		matrix[1]=PIJ_TV_C;
		matrix[2]=PIJ_TS_G;
		matrix[3]=PIJ_TV_T;

		matrix[4]=PIJ_TV_A;
		matrix[5]=PIJ_SAME_C;
		matrix[6]=PIJ_TV_G;
		matrix[7]=PIJ_TS_T;
		
		matrix[8]=PIJ_TS_A;
		matrix[9]=matrix[1];  /* PIJ_TV_C */
		matrix[10]=PIJ_SAME_G;
		matrix[11]=matrix[3]; /* PIJ_TV_T */
		
		matrix[12]=matrix[4]; /* PIJ_TV_A */
		matrix[13]=PIJ_TS_C;
		matrix[14]=matrix[6]; /* PIJ_TV_G */
		matrix[15]=PIJ_SAME_T;
	}	

	
/* the rows are cumulative to help with picking one using
   a random number */
	matrix[1]+=matrix[0];
	matrix[2]+=matrix[1];
	matrix[3]+=matrix[2]; /* This should always be 1.0... */

	matrix[5]+=matrix[4];
	matrix[6]+=matrix[5];
	matrix[7]+=matrix[6]; /* ...but it is easier to spot bugs... */
	
	matrix[9]+=matrix[8];
	matrix[10]+=matrix[9];
	matrix[11]+=matrix[10]; /* ...though less efficient... */
	
	matrix[13]+=matrix[12];
	matrix[14]+=matrix[13];
	matrix[15]+=matrix[14]; /* ...but probably not much. */
}


void NetEvolvePlugin::SetVector(double *vector, short base, double len)
{	
	double aa, bbR, bbY;
	int i,j,k;
	double expt[4];
	double *P;

	if (model==REV) {
		P=vector;
		if (len<1e-6) { 
			for (i=0; i<4; i++) {
				if (i==base)
					*P=1.0;
				else 	
					*P=0.0;
				P++;
			}
			return; 
		}
		for (k=1; k<4; k++) 
			expt[k]=exp(len*Root[k]);

		for (j=0; j<4; j++) {
			(*P)=Cijk[base*4*4+j*4+0];
			for (k=1; k<4; k++)
				(*P)+=Cijk[base*4*4+j*4+k]*expt[k];
			P++;
		}
		
		vector[1]+=vector[0];
		vector[2]+=vector[1];
		vector[3]+=vector[2];
	} else {
		switch (model) {
			case F84:
				aa=exp(mu*len);
				bbR=bbY=exp(mu_kappa_1*len);
			break;
			case HKY:
				aa=exp(beta*len);
				bbR=exp(beta_A_R*len);
				bbY=exp(beta_A_Y*len);
			break;
		}
		
		switch (base) {
			case 0:
				vector[0]=PIJ_SAME_A;
				vector[1]=PIJ_TV_C+vector[0];
				vector[2]=PIJ_TS_G+vector[1];
				vector[3]=PIJ_TV_T+vector[2];
			break;
			case 1:
				vector[0]=PIJ_TV_A;
				vector[1]=PIJ_SAME_C+vector[0];
				vector[2]=PIJ_TV_G+vector[1];
				vector[3]=PIJ_TS_T+vector[2];
			break;
			case 2:
				vector[0]=PIJ_TS_A;
				vector[1]=PIJ_TV_C+vector[0];
				vector[2]=PIJ_SAME_G+vector[1];
				vector[3]=PIJ_TV_T+vector[2];
			break;
			case 3:
				vector[0]=PIJ_TV_A;
				vector[1]=PIJ_TS_C+vector[0];
				vector[2]=PIJ_TV_G+vector[1];
				vector[3]=PIJ_SAME_T+vector[2];
			break;
		}
	}
}


/* Everything below is shamelessly taken from Yang's Paml package */


int NetEvolvePlugin::EigenREV (double Root[], double Cijk[])
{
/* freq[] is constant
*/
	int i,j,k;
	double U[16], V[16], T1[16], T2[16];

	abyx (1/mr, Qij, 16);

	if ((k=eigen (1, Qij, 4, Root, T1, U, V, T2))!=0) {
		fprintf(stderr, "\ncomplexx roots in EigenREV");
		exit(0);
	}
	xtoy (U, V, 16);
	matinv (V, 4, 4, T1);
	for (i=0; i<4; i++) 
   		for (j=0; j<4; j++) 
   			for (k=0; k<4; k++) 
   				Cijk[i*4*4+j*4+k] = U[i*4+k]*V[k*4+j];
	return (0);
}
   
int NetEvolvePlugin::abyx (double a, double x[], int n)
{ int i; for (i=0; i<n; x[i]*=a,i++) ;  return(0); }
int NetEvolvePlugin::xtoy (double x[], double y[], int n)
{ int i; for (i=0; i<n; y[i]=x[i],i++) ;  return(0); }
int NetEvolvePlugin::matinv( double x[], int n, int m, double space[])
{
/* x[n*m]  ... m>=n
*/
   register int i,j,k;
   int *irow=(int*) space;
   double ee=1.0e-20, t,t1,xmax;
   double det=1.0;

   for (i=0; i<n; i++)  {
      xmax = 0.;
      for (j=i; j<n; j++) {
	 if (xmax < fabs(x[j*m+i]))  {
	    xmax = fabs( x[j*m+i] );
	    irow[i] = j;
	 }
      }
      det *= xmax;
      if (xmax < ee)   {
	 printf("\nDet becomes zero at %3d!\t\n", i+1);
	 return(-1);
      }
      if (irow[i] != i) {
	 for (j=0; j<m; j++) {
	    t = x[i*m+j];
	    x[i*m+j] = x[irow[i] * m + j];
	    x[ irow[i] * m + j] = t;
	 }
      }
      t = 1./x[i*m+i];
      for (j=0; j<n; j++) {
	 if (j == i) continue;
	 t1 = t*x[j*m+i];
	 for (k=0; k<n; k++)  x[j*m+k] -= t1*x[i*m+k];
	 x[j*m+i] = -t1;
      }
      for (j=0; j<m; j++)   x[i*m+j] *= t;
      x[i*m+i] = t;
   }                            /* i  */
   for (i=n-1; i>=0; i--) {
      if (irow[i] == i) continue;
      for (j=0; j<n; j++)  {
	 t = x[j*m+i];
	 x[j*m+i] = x[ j*m + irow[i] ];
	 x[ j*m + irow[i] ] = t;
      }
   }
   return (0);
}

/***********************************************************
*  This eigen() works for eigenvalue/vector analysis
*         for real general square matrix A
*         A will be destroyed
*         rr,ri are vectors containing eigenvalues
*         vr,vi are matrices containing (right) eigenvectors
*
*              A*[vr+vi*i] = [vr+vi*i] * diag{rr+ri*i}
*
*  Algorithm: Handbook for Automatic Computation, vol 2
*             by Wilkinson and Reinsch, 1971
*             most of source codes were taken from a public domain
*             solftware called MATCALC.
*  Credits:   to the authors of MATCALC
*
*  return     -1 not converged
*              0 no complexx eigenvalues/vectors
*              1 complexx eigenvalues/vectors
*  Tianlin Wang at University of Illinois
*  Thu May  6 15:22:31 CDT 1993
***************************************************************/

int NetEvolvePlugin::eigen(int job, double A[], int n, double rr[], double ri[], 
          double vr[], double vi[], double work[])
{    
/*  double work[n*2]: working space
*/
    int low,hi,i,j,k, it, istate=0;
    double tiny=sqrt(pow((double)BASE,(double)(1-DIGITS))), t; 

    balance(A,n,&low,&hi,work);
    elemhess(job,A,n,low,hi,vr,vi, (int*)(work+n));
    if (-1 == realeig(job,A,n,low,hi,rr,ri,vr,vi)) return (-1);
    if (job) unbalance(n,vr,vi,low,hi,work);

/* sort, added by Z. Yang */
   for (i=0; i<n; i++) {
       for (j=i+1,it=i,t=rr[i]; j<n; j++)
           if (t<rr[j]) { t=rr[j]; it=j; }
       rr[it]=rr[i];   rr[i]=t;
       t=ri[it];       ri[it]=ri[i];  ri[i]=t;
       for (k=0; k<n; k++) {
          t=vr[k*n+it];  vr[k*n+it]=vr[k*n+i];  vr[k*n+i]=t;
          t=vi[k*n+it];  vi[k*n+it]=vi[k*n+i];  vi[k*n+i]=t;
       }
       if (fabs(ri[i])>tiny) istate=1;
   }

    return (istate) ;
}

/* complexx funcctions
*/

complexx NetEvolvePlugin::myc (double re,double im)
{
    complexx r;

    r.re = re;
    r.im = im;
    return(r);
}

complexx NetEvolvePlugin::conj (complexx a)
{
    a.im = -a.im;
    return(a);
}

#define csize(a) (fabs(a.re)+fabs(a.im))

complexx NetEvolvePlugin::cplus (complexx a, complexx b)
{
   complexx c;
   c.re = a.re+b.re;  
   c.im = a.im+b.im;   
   return (c);
}

complexx NetEvolvePlugin::cminus (complexx a, complexx b)
{
   complexx c;
   c.re = a.re-b.re;  
   c.im = a.im-b.im;   
   return (c);
}

complexx NetEvolvePlugin::cby (complexx a, complexx b)
{
   complexx c;
   c.re = a.re*b.re-a.im*b.im ;
   c.im = a.re*b.im+a.im*b.re ;
   return (c);
}

complexx NetEvolvePlugin::cdiv (complexx a,complexx b)
{
    double ratio, den;
    complexx c;

    if (fabs(b.re) <= fabs(b.im)) {
        ratio = b.re / b.im;
        den = b.im * (1 + ratio * ratio);
        c.re = (a.re * ratio + a.im) / den;
        c.im = (a.im * ratio - a.re) / den;
    }
    else {
        ratio = b.im / b.re;
        den = b.re * (1 + ratio * ratio);
        c.re = (a.re + a.im * ratio) / den;
        c.im = (a.im - a.re * ratio) / den;
    }
    return(c);
}

complexx NetEvolvePlugin::cexp (complexx a)
{
   complexx c;
   c.re = exp(a.re);
   if (fabs(a.im)==0) c.im = 0; 
   else  { c.im = c.re*sin(a.im); c.re*=cos(a.im); }
   return (c);
}

complexx NetEvolvePlugin::cfactor (complexx x, double a)
{
   complexx c;
   c.re = a*x.re; 
   c.im = a*x.im;
   return (c);
}

int NetEvolvePlugin::cxtoy (complexx x[], complexx y[], int n)
{
   int i;
   FOR (i,n) y[i]=x[i];
   return (0);
}

int NetEvolvePlugin::cmatby (complexx a[], complexx b[], complexx c[], int n,int m,int k)
/* a[n*m], b[m*k], c[n*k]  ......  c = a*b
*/
{
   int i,j,i1;
   complexx t;

   FOR (i,n)  FOR(j,k) {
       for (i1=0,t=myc(0,0); i1<m; i1++)  
           t = cplus (t, cby(a[i*m+i1],b[i1*k+j]));
       c[i*k+j] = t;
   }
   return (0);
}

int NetEvolvePlugin::cmatout (FILE * fout, complexx x[], int n, int m)
{
   int i,j;
   for (i=0,FPN(fout); i<n; i++,FPN(fout)) 
        FOR(j,m) fprintf(fout, "%7.3f%7.3f  ", x[i*m+j].re, x[i*m+j].im);
   return (0);
}

int NetEvolvePlugin::cmatinv( complexx x[], int n, int m, double space[])
{
/* x[n*m]  ... m>=n
*/
   int i,j,k, *irow=(int*) space;
   double xmaxsize, ee=1e-20;
   complexx xmax, t,t1;

   FOR(i,n)  {
       xmaxsize = 0.;
       for (j=i; j<n; j++) {
          if ( xmaxsize < csize (x[j*m+i]))  {
               xmaxsize = csize (x[j*m+i]);
               xmax = x[j*m+i];
               irow[i] = j;
          }
       }
       if (xmaxsize < ee)   {
           printf("\nDet goes to zero at %8d!\t\n", i+1);
           return(-1);
       }
       if (irow[i] != i) {
           FOR(j,m) {
                t = x[i*m+j];
                x[i*m+j] = x[irow[i]*m+j];
                x[ irow[i]*m+j] = t;
           }
       }
       t = cdiv (myc(1,0), x[i*m+i]);
       FOR(j,n) {
           if (j == i) continue;
           t1 = cby (t,x[j*m+i]);
           FOR(k,m)  x[j*m+k] = cminus (x[j*m+k], cby(t1,x[i*m+k]));
           x[j*m+i] = cfactor (t1, -1);
       }
       FOR(j,m)   x[i*m+j] = cby (x[i*m+j], t);
       x[i*m+i] = t;
   }                         
   for (i=n-1; i>=0; i--) {
        if (irow[i] == i) continue;
        FOR(j,n)  {
           t = x[j*m+i];
           x[j*m+i] = x[j*m+irow[i]];
           x[ j*m+irow[i]] = t;
        }
   }
   return (0);
}


void NetEvolvePlugin::balance(double mat[], int n,int *low, int *hi, double scale[])
{
/* Balance a matrix for calculation of eigenvalues and eigenvectors
*/
    double c,f,g,r,s;
    int i,j,k,l,done;
        /* search for rows isolating an eigenvalue and push them down */
    for (k = n - 1; k >= 0; k--) {
        for (j = k; j >= 0; j--) {
            for (i = 0; i <= k; i++) {
                if (i != j && fabs(mat[pos(j,i,n)]) != 0) break;
            }

            if (i > k) {
                scale[k] = j;

                if (j != k) {
                    for (i = 0; i <= k; i++) {
                       c = mat[pos(i,j,n)];
                       mat[pos(i,j,n)] = mat[pos(i,k,n)];
                       mat[pos(i,k,n)] = c;
                    }

                    for (i = 0; i < n; i++) {
                       c = mat[pos(j,i,n)];
                       mat[pos(j,i,n)] = mat[pos(k,i,n)];
                       mat[pos(k,i,n)] = c;
                    }
                }
                break;
            }
        }
        if (j < 0) break;
    }

    /* search for columns isolating an eigenvalue and push them left */

    for (l = 0; l <= k; l++) {
        for (j = l; j <= k; j++) {
            for (i = l; i <= k; i++) {
                if (i != j && fabs(mat[pos(i,j,n)]) != 0) break;
            }
            if (i > k) {
                scale[l] = j;
                if (j != l) {
                    for (i = 0; i <= k; i++) {
                       c = mat[pos(i,j,n)];
                       mat[pos(i,j,n)] = mat[pos(i,l,n)];
                       mat[pos(i,l,n)] = c;
                    }

                    for (i = l; i < n; i++) {
                       c = mat[pos(j,i,n)];
                       mat[pos(j,i,n)] = mat[pos(l,i,n)];
                       mat[pos(l,i,n)] = c;
                    }
                }

                break;
            }
        }

        if (j > k) break;
    }

    *hi = k;
    *low = l;

    /* balance the submatrix in rows l through k */

    for (i = l; i <= k; i++) {
        scale[i] = 1;
    }

    do {
        for (done = 1,i = l; i <= k; i++) {
            for (c = 0,r = 0,j = l; j <= k; j++) {
                if (j != i) {
                    c += fabs(mat[pos(j,i,n)]);
                    r += fabs(mat[pos(i,j,n)]);
                }
            }

            if (c != 0 && r != 0) {
                g = r / BASE;
                f = 1;
                s = c + r;

                while (c < g) {
                    f *= BASE;
                    c *= BASE * BASE;
                }

                g = r * BASE;

                while (c >= g) {
                    f /= BASE;
                    c /= BASE * BASE;
                }

                if ((c + r) / f < 0.95 * s) {
                    done = 0;
                    g = 1 / f;
                    scale[i] *= f;

                    for (j = l; j < n; j++) {
                        mat[pos(i,j,n)] *= g;
                    }

                    for (j = 0; j <= k; j++) {
                        mat[pos(j,i,n)] *= f;
                    }
                }
            }
        }
    } while (!done);
}


/*
 * Transform back eigenvectors of a balanced matrix
 * into the eigenvectors of the original matrix
 */
void NetEvolvePlugin::unbalance(int n,double vr[],double vi[], int low, int hi, double scale[])
{
    int i,j,k;
    double tmp;

    for (i = low; i <= hi; i++) {
        for (j = 0; j < n; j++) {
            vr[pos(i,j,n)] *= scale[i];
            vi[pos(i,j,n)] *= scale[i];
        }
    }

    for (i = low - 1; i >= 0; i--) {
        if ((k = (int)scale[i]) != i) {
            for (j = 0; j < n; j++) {
                tmp = vr[pos(i,j,n)];
                vr[pos(i,j,n)] = vr[pos(k,j,n)];
                vr[pos(k,j,n)] = tmp;

                tmp = vi[pos(i,j,n)];
                vi[pos(i,j,n)] = vi[pos(k,j,n)];
                vi[pos(k,j,n)] = tmp;        
            }
        }
    }

    for (i = hi + 1; i < n; i++) {
        if ((k = (int)scale[i]) != i) {
            for (j = 0; j < n; j++) {
                tmp = vr[pos(i,j,n)];
                vr[pos(i,j,n)] = vr[pos(k,j,n)];
                vr[pos(k,j,n)] = tmp;

                tmp = vi[pos(i,j,n)];
                vi[pos(i,j,n)] = vi[pos(k,j,n)];
                vi[pos(k,j,n)] = tmp;        
            }
        }
    }
}

/*
 * Reduce the submatrix in rows and columns low through hi of real matrix mat to
 * Hessenberg form by elementary similarity transformations
 */
void NetEvolvePlugin::elemhess(int job,double mat[],int n,int low,int hi, double vr[],
              double vi[], int work[])
{
/* work[n] */
    int i,j,m;
    double x,y;

    for (m = low + 1; m < hi; m++) {
        for (x = 0,i = m,j = m; j <= hi; j++) {
            if (fabs(mat[pos(j,m-1,n)]) > fabs(x)) {
                x = mat[pos(j,m-1,n)];
                i = j;
            }
        }

        if ((work[m] = i) != m) {
            for (j = m - 1; j < n; j++) {
               y = mat[pos(i,j,n)];
               mat[pos(i,j,n)] = mat[pos(m,j,n)];
               mat[pos(m,j,n)] = y;
            }

            for (j = 0; j <= hi; j++) {
               y = mat[pos(j,i,n)];
               mat[pos(j,i,n)] = mat[pos(j,m,n)];
               mat[pos(j,m,n)] = y;
            }
        }

        if (x != 0) {
            for (i = m + 1; i <= hi; i++) {
                if ((y = mat[pos(i,m-1,n)]) != 0) {
                    y = mat[pos(i,m-1,n)] = y / x;

                    for (j = m; j < n; j++) {
                        mat[pos(i,j,n)] -= y * mat[pos(m,j,n)];
                    }

                    for (j = 0; j <= hi; j++) {
                        mat[pos(j,m,n)] += y * mat[pos(j,i,n)];
                    }
                }
            }
        }
    }
    if (job) {
       for (i=0; i<n; i++) {
          for (j=0; j<n; j++) {
             vr[pos(i,j,n)] = 0.0; vi[pos(i,j,n)] = 0.0;
          }
          vr[pos(i,i,n)] = 1.0;
       }

       for (m = hi - 1; m > low; m--) {
          for (i = m + 1; i <= hi; i++) {
             vr[pos(i,m,n)] = mat[pos(i,m-1,n)];
          }

         if ((i = work[m]) != m) {
            for (j = m; j <= hi; j++) {
               vr[pos(m,j,n)] = vr[pos(i,j,n)];
               vr[pos(i,j,n)] = 0.0;
            }
            vr[pos(i,m,n)] = 1.0;
         }
      }
   }
}

/*
 * Calculate eigenvalues and eigenvectors of a real upper Hessenberg matrix
 * Return 1 if converges successfully and 0 otherwise
 */
 
int NetEvolvePlugin::realeig(int job,double mat[],int n,int low, int hi, double valr[],
      double vali[], double vr[],double vi[])
{
   complexx v;
   double p=0,q=0,r=0,s=0,t,w,x,y,z=0,ra,sa,norm,eps;
   int niter,en,i,j,k,l,m;
   double precision  = pow((double)BASE,(double)(1-DIGITS));

   eps = precision;
   for (i=0; i<n; i++) {
      valr[i]=0.0;
      vali[i]=0.0;
   }
      /* store isolated roots and calculate norm */
   for (norm = 0,i = 0; i < n; i++) {
      for (j = max(0,i-1); j < n; j++) {
         norm += fabs(mat[pos(i,j,n)]);
      }
      if (i < low || i > hi) valr[i] = mat[pos(i,i,n)];
   }
   t = 0;
   en = hi;

   while (en >= low) {
      niter = 0;
      for (;;) {

       /* look for single small subdiagonal element */

         for (l = en; l > low; l--) {
            s = fabs(mat[pos(l-1,l-1,n)]) + fabs(mat[pos(l,l,n)]);
            if (s == 0) s = norm;
            if (fabs(mat[pos(l,l-1,n)]) <= eps * s) break;
         }

         /* form shift */

         x = mat[pos(en,en,n)];

         if (l == en) {             /* one root found */
            valr[en] = x + t;
            if (job) mat[pos(en,en,n)] = x + t;
            en--;
            break;
         }

         y = mat[pos(en-1,en-1,n)];
         w = mat[pos(en,en-1,n)] * mat[pos(en-1,en,n)];

         if (l == en - 1) {                /* two roots found */
            p = (y - x) / 2;
            q = p * p + w;
            z = sqrt(fabs(q));
            x += t;
            if (job) {
               mat[pos(en,en,n)] = x;
               mat[pos(en-1,en-1,n)] = y + t;
            }
            if (q < 0) {                /* complexx pair */
               valr[en-1] = x+p;
               vali[en-1] = z;
               valr[en] = x+p;
               vali[en] = -z;
            }
            else {                      /* real pair */
               z = (p < 0) ? p - z : p + z;
               valr[en-1] = x + z;
               valr[en] = (z == 0) ? x + z : x - w / z;
               if (job) {
                  x = mat[pos(en,en-1,n)];
                  s = fabs(x) + fabs(z);
                  p = x / s;
                  q = z / s;
                  r = sqrt(p*p+q*q);
                  p /= r;
                  q /= r;
                  for (j = en - 1; j < n; j++) {
                     z = mat[pos(en-1,j,n)];
                     mat[pos(en-1,j,n)] = q * z + p *
                     mat[pos(en,j,n)];
                     mat[pos(en,j,n)] = q * mat[pos(en,j,n)] - p*z;
                  }
                  for (i = 0; i <= en; i++) {
                     z = mat[pos(i,en-1,n)];
                     mat[pos(i,en-1,n)] = q * z + p * mat[pos(i,en,n)];
                     mat[pos(i,en,n)] = q * mat[pos(i,en,n)] - p*z;
                  }
                  for (i = low; i <= hi; i++) {
                     z = vr[pos(i,en-1,n)];
                     vr[pos(i,en-1,n)] = q*z + p*vr[pos(i,en,n)];
                     vr[pos(i,en,n)] = q*vr[pos(i,en,n)] - p*z;
                  }
               }
            }
            en -= 2;
            break;
         }
         if (niter == MAXITER) return(-1);
         if (niter != 0 && niter % 10 == 0) {
            t += x;
            for (i = low; i <= en; i++) mat[pos(i,i,n)] -= x;
            s = fabs(mat[pos(en,en-1,n)]) + fabs(mat[pos(en-1,en-2,n)]);
            x = y = 0.75 * s;
            w = -0.4375 * s * s;
         }
         niter++;
           /* look for two consecutive small subdiagonal elements */
         for (m = en - 2; m >= l; m--) {
            z = mat[pos(m,m,n)];
            r = x - z;
            s = y - z;
            p = (r * s - w) / mat[pos(m+1,m,n)] + mat[pos(m,m+1,n)];
            q = mat[pos(m+1,m+1,n)] - z - r - s;
            r = mat[pos(m+2,m+1,n)];
            s = fabs(p) + fabs(q) + fabs(r);
            p /= s;
            q /= s;
            r /= s;
            if (m == l || fabs(mat[pos(m,m-1,n)]) * (fabs(q)+fabs(r)) <=
                eps * (fabs(mat[pos(m-1,m-1,n)]) + fabs(z) +
                fabs(mat[pos(m+1,m+1,n)])) * fabs(p)) break;
         }
         for (i = m + 2; i <= en; i++) mat[pos(i,i-2,n)] = 0;
         for (i = m + 3; i <= en; i++) mat[pos(i,i-3,n)] = 0;
             /* double QR step involving rows l to en and columns m to en */
         for (k = m; k < en; k++) {
            if (k != m) {
               p = mat[pos(k,k-1,n)];
               q = mat[pos(k+1,k-1,n)];
               r = (k == en - 1) ? 0 : mat[pos(k+2,k-1,n)];
               if ((x = fabs(p) + fabs(q) + fabs(r)) == 0) continue;
               p /= x;
               q /= x;
               r /= x;
            }
            s = sqrt(p*p+q*q+r*r);
            if (p < 0) s = -s;
            if (k != m) {
               mat[pos(k,k-1,n)] = -s * x;
            }
            else if (l != m) {
               mat[pos(k,k-1,n)] = -mat[pos(k,k-1,n)];
            }
            p += s;
            x = p / s;
            y = q / s;
            z = r / s;
            q /= p;
            r /= p;
                /* row modification */
            for (j = k; j <= (!job ? en : n-1); j++){
               p = mat[pos(k,j,n)] + q * mat[pos(k+1,j,n)];
               if (k != en - 1) {
                  p += r * mat[pos(k+2,j,n)];
                  mat[pos(k+2,j,n)] -= p * z;
               }
               mat[pos(k+1,j,n)] -= p * y;
               mat[pos(k,j,n)] -= p * x;
            }
            j = min(en,k+3);
              /* column modification */
            for (i = (!job ? l : 0); i <= j; i++) {
               p = x * mat[pos(i,k,n)] + y * mat[pos(i,k+1,n)];
               if (k != en - 1) {
                  p += z * mat[pos(i,k+2,n)];
                  mat[pos(i,k+2,n)] -= p*r;
               }
               mat[pos(i,k+1,n)] -= p*q;
               mat[pos(i,k,n)] -= p;
            }
            if (job) {             /* accumulate transformations */
               for (i = low; i <= hi; i++) {
                  p = x * vr[pos(i,k,n)] + y * vr[pos(i,k+1,n)];
                  if (k != en - 1) {
                     p += z * vr[pos(i,k+2,n)];
                     vr[pos(i,k+2,n)] -= p*r;
                  }
                  vr[pos(i,k+1,n)] -= p*q;
                  vr[pos(i,k,n)] -= p;
               }
            }
         }
      }
   }

   if (!job) return(0);
   if (norm != 0) {
       /* back substitute to find vectors of upper triangular form */
      for (en = n-1; en >= 0; en--) {
         p = valr[en];
         if ((q = vali[en]) < 0) {            /* complexx vector */
            m = en - 1;
            if (fabs(mat[pos(en,en-1,n)]) > fabs(mat[pos(en-1,en,n)])) {
               mat[pos(en-1,en-1,n)] = q / mat[pos(en,en-1,n)];
               mat[pos(en-1,en,n)] = (p - mat[pos(en,en,n)]) /
                     mat[pos(en,en-1,n)];
            }
            else {
               v = cdiv(myc(0.0,-mat[pos(en-1,en,n)]),
                    myc(mat[pos(en-1,en-1,n)]-p,q));
               mat[pos(en-1,en-1,n)] = v.re;
               mat[pos(en-1,en,n)] = v.im;
            }
            mat[pos(en,en-1,n)] = 0;
            mat[pos(en,en,n)] = 1;
            for (i = en - 2; i >= 0; i--) {
               w = mat[pos(i,i,n)] - p;
               ra = 0;
               sa = mat[pos(i,en,n)];
               for (j = m; j < en; j++) {
                  ra += mat[pos(i,j,n)] * mat[pos(j,en-1,n)];
                  sa += mat[pos(i,j,n)] * mat[pos(j,en,n)];
               }
               if (vali[i] < 0) {
                  z = w;
                  r = ra;
                  s = sa;
               }
               else {
                  m = i;
                  if (vali[i] == 0) {
                     v = cdiv(myc(-ra,-sa),myc(w,q));
                     mat[pos(i,en-1,n)] = v.re;
                     mat[pos(i,en,n)] = v.im;
                  }
                  else {                      /* solve complexx equations */
                     x = mat[pos(i,i+1,n)];
                     y = mat[pos(i+1,i,n)];
                     v.re = (valr[i]- p)*(valr[i]-p) + vali[i]*vali[i] - q*q;
                     v.im = (valr[i] - p)*2*q;
                     if ((fabs(v.re) + fabs(v.im)) == 0) {
                        v.re = eps * norm * (fabs(w) +
                                fabs(q) + fabs(x) + fabs(y) + fabs(z));
                     }
                     v = cdiv(myc(x*r-z*ra+q*sa,x*s-z*sa-q*ra),v);
                     mat[pos(i,en-1,n)] = v.re;
                     mat[pos(i,en,n)] = v.im;
                     if (fabs(x) > fabs(z) + fabs(q)) {
                        mat[pos(i+1,en-1,n)] = 
                             (-ra - w * mat[pos(i,en-1,n)] +
                             q * mat[pos(i,en,n)]) / x;
                        mat[pos(i+1,en,n)] = (-sa - w * mat[pos(i,en,n)] -
                             q * mat[pos(i,en-1,n)]) / x;
                     }
                     else {
                        v = cdiv(myc(-r-y*mat[pos(i,en-1,n)],
                             -s-y*mat[pos(i,en,n)]),myc(z,q));
                        mat[pos(i+1,en-1,n)] = v.re;
                        mat[pos(i+1,en,n)] = v.im;
                     }
                  }
               }
            }
         }
         else if (q == 0) {                             /* real vector */
            m = en;
            mat[pos(en,en,n)] = 1;
            for (i = en - 1; i >= 0; i--) {
               w = mat[pos(i,i,n)] - p;
               r = mat[pos(i,en,n)];
               for (j = m; j < en; j++) {
                  r += mat[pos(i,j,n)] * mat[pos(j,en,n)];
               }
               if (vali[i] < 0) {
                  z = w;
                  s = r;
               }
               else {
                  m = i;
                  if (vali[i] == 0) {
                     if ((t = w) == 0) t = eps * norm;
                     mat[pos(i,en,n)] = -r / t;
                  }
                  else {            /* solve real equations */
                     x = mat[pos(i,i+1,n)];
                     y = mat[pos(i+1,i,n)];
                     q = (valr[i] - p) * (valr[i] - p) + vali[i]*vali[i];
                     t = (x * s - z * r) / q;
                     mat[pos(i,en,n)] = t;
                     if (fabs(x) <= fabs(z)) {
                        mat[pos(i+1,en,n)] = (-s - y * t) / z;
                     }
                     else {
                        mat[pos(i+1,en,n)] = (-r - w * t) / x;
                     }
                  }
               }
            }
         }
      }
             /* vectors of isolated roots */
      for (i = 0; i < n; i++) {
         if (i < low || i > hi) {
            for (j = i; j < n; j++) {
               vr[pos(i,j,n)] = mat[pos(i,j,n)];
            }
         }
      }
       /* multiply by transformation matrix */

      for (j = n-1; j >= low; j--) {
         m = min(j,hi);
         for (i = low; i <= hi; i++) {
            for (z = 0,k = low; k <= m; k++) {
               z += vr[pos(i,k,n)] * mat[pos(k,j,n)];
            }
            vr[pos(i,j,n)] = z;
         }
      }
   }
    /* rearrange complexx eigenvectors */
   for (j = 0; j < n; j++) {
      if (vali[j] != 0) {
         for (i = 0; i < n; i++) {
            vi[pos(i,j,n)] = vr[pos(i,j+1,n)];
            vr[pos(i,j+1,n)] = vr[pos(i,j,n)];
            vi[pos(i,j+1,n)] = -vi[pos(i,j,n)];
         }
         j++;
      }
   }
   return(0);
}
void NetEvolvePlugin::SetSeed(unsigned long seed)
{
	int N = 624;
    /* setting initial seeds to mt[N] using         */
    /* the generator Line 25 of Table 1 in          */
    /* [KNUTH 1981, The Art of Computer Programming */
    /*    Vol. 2 (2nd Ed.), pp102]                  */
    mt[0]= seed & 0xffffffff;
    for (mti=1; mti<N; mti++)
        mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}
//void SetSeed (int seed)
//{
//   z_rndu = 170*(seed%178) + 137;
//   w_rndu=seed;
//}


//#ifdef FAST_RANDOM_NUMBER

//double rndu (void)
//{
//   w_rndu *= 127773;
//   return ldexp((double)w_rndu, -32);
//}

//#else 

//double rndu (void)
//{
//   static int x_rndu=11, y_rndu=23;
//   double r;

//   x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
//   y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
//   z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
//   if (x_rndu<0) x_rndu+=30269;
//   if (y_rndu<0) y_rndu+=30307;
//   if (z_rndu<0) z_rndu+=30323;
//   r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
//   return (r-(int)r);
//}

//#endif
double NetEvolvePlugin::rndu()
{
	int N=624, M=397;
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if sgenrand() has not been called, */
            SetSeed(4357); /* a default initial seed is used   */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    return ( (double)y / (unsigned long)0xffffffff ); /* reals */
    /* return y; */ /* for integer generation */
}





double NetEvolvePlugin::rndgamma (double s)
{
	double	r=0.0;
	
	if (s <= 0.0)      
		puts ("Error gamma..");
	else if (s < 1.0)  
		r = rndgamma1 (s);
	else if (s > 1.0)  
		r = rndgamma2 (s);
	else           
		r =- log(rndu());
	return (r);
}


double NetEvolvePlugin::rndgamma1 (double s)
{

	double			r, x=0.0, small=1e-37, w;
	static double	a, p, uf, ss=10.0, d;
	
	if (s!=ss) 
		{
		a  = 1.0-s;
		p  = a/(a+s*exp(-a));
		uf = p*pow(small/a,s);
		d  = a*log(a);
		ss = s;
		}
	for (;;) 
		{
		r = rndu();
		if (r > p)        
			x = a-log((1.0-r)/(1.0-p)), w=a*log(x)-d;
		else if (r>uf)  
			x = a*pow(r/p,1/s), w=x;
		else            
			return (0.0);
		r = rndu();
		if (1.0-r <= w && r > 0.0)
			if (r*(w+1.0) >= 1.0 || -log(r) <= w)  
				continue;
		break;
		}
	return (x);
}


double NetEvolvePlugin::rndgamma2 (double s)
{

	double			r ,d, f, g, x;
	static double	b, h, ss=0;
	
	if (s!=ss) 
		{
		b  = s-1.0;
		h  = sqrt(3.0*s-0.75);
		ss = s;
		}
	for (;;) 
		{
		r = rndu();
		g = r-r*r;
		f = (r-0.5)*h/sqrt(g);
		x = b+f;
		if (x <= 0.0) 
			continue;
		r = rndu();
		d = 64*r*r*g*g*g;
		if (d*x < x-2.0*f*f || log(d) < 2*(b*log(x/b)-f))  
			break;
		}
	return (x);
}

double NetEvolvePlugin::LnGamma (double alpha)
{
/* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.  
   Stirling's formula is used for the central polynomial part of the procedure.
   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double x=alpha, f=0, z;

   if (x<7) {
      f=1;  z=x-1;
      while (++z<7)  f*=z;
      x=z;   f=-log(f);
   }
   z = 1/(x*x);
   return  f + (x-0.5)*log(x) - x + .918938533204673 
	  + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
	       +.083333333333333)/x;  
}

double NetEvolvePlugin::IncompleteGamma (double x, double alpha, double ln_gamma_alpha)
{
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper 
	   limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   double p=alpha, g=ln_gamma_alpha;
   double accurate=1e-8, overflow=1e30;
   double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];

   if (x==0) return (0);
   if (x<0 || p<=0) return (-1);

   factor=exp(p*log(x)-x-g);   
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;

   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a=1-p;   b=a+x+1;  term=0;
   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
   gin=pn[2]/pn[3];
 l32:
   a++;  b+=2;  term++;   an=a*term;
   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
   if (pn[5] == 0) goto l35;
   rn=pn[4]/pn[5];   dif=fabs(gin-rn);
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:
   gin=rn;
 l35:
   for (i=0; i<4; i++) pn[i]=pn[i+2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i]/=overflow;
   goto l32;
 l42:
   gin=1-factor*gin;

 l50:
   return (gin);
}


/* functions concerning the CDF and percentage points of the gamma and
   Chi2 distribution
*/
double NetEvolvePlugin::PointNormal (double prob)
{
/* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
   returns (-9999) if in error
   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
   Applied Statistics 22: 96-97 (AS70)

   Newer methods:
     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
       normal distribution.  37: 477-484.
     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage 
       points of the normal distribution.  26: 118-121.

*/
   double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
   double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
   double b2=.531103462366, b3=.103537752850, b4=.0038560700634;
   double y, z=0, p=prob, p1;

   p1 = (p<0.5 ? p : 1-p);
   if (p1<1e-20) return (-9999);

   y = sqrt (log(1/(p1*p1)));   
   z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
   return (p<0.5 ? -z : z);
}


double NetEvolvePlugin::PointChi2 (double prob, double v)
{
/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
   returns -1 if in error.   0.000002<prob<0.999998
   RATNEST FORTRAN by
       Best DJ & Roberts DE (1975) The percentage points of the 
       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
   Converted into C by Ziheng Yang, Oct. 1993.
*/
   double e=.5e-6, aa=.6931471805, p=prob, g;
   double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;

   if (p<.000002 || p>.999998 || v<=0) return (-1);

   g = LnGamma (v/2);
   xx=v/2;   c=xx-1;
   if (v >= -1.24*log(p)) goto l1;

   ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
   if (ch-e<0) return (ch);
   goto l4;
l1:
   if (v>.32) goto l3;
   ch=0.4;   a=log(1-p);
l2:
   q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
   t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
   ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
   if (fabs(q/ch-1)-.01 <= 0) goto l4;
   else                       goto l2;
  
l3: 
   x=PointNormal (p);
   p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
   if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
l4:
   q=ch;   p1=.5*ch;
   if ((t=IncompleteGamma (p1, xx, g))<0) {
      return (-1);
   }
   p2=p-t;
   t=p2*exp(xx*aa+g+p1-c*log(ch));   
   b=t/ch;  a=0.5*t-b*c;

   s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
   s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
   s3=(210+a*(462+a*(707+932*a)))/2520;
   s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
   s5=(84+264*a+c*(175+606*a))/2520;
   s6=(120+c*(346+127*c))/5040;
   ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
   if (fabs(q/ch-1) > e) goto l4;

   return (ch);
}


#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))

int NetEvolvePlugin::DiscreteGamma (double freqK[], double rK[], 
    double alfa, double beta, int K, int median)
{
/* discretization of gamma distribution with equal proportions in each 
   category
*/
   int i;
   double gap05=1.0/(2.0*K), t, factor=alfa/beta*K, lnga1;

   if (median) {
      for (i=0; i<K; i++) rK[i]=PointGamma((i*2.0+1)*gap05, alfa, beta);
      for (i=0,t=0; i<K; i++) t+=rK[i];
      for (i=0; i<K; i++)     rK[i]*=factor/t;
   }
   else {
      lnga1=LnGamma(alfa+1);
      for (i=0; i<K-1; i++)
	 freqK[i]=PointGamma((i+1.0)/K, alfa, beta);
      for (i=0; i<K-1; i++)
	 freqK[i]=IncompleteGamma(freqK[i]*beta, alfa+1, lnga1);
      rK[0] = freqK[0]*factor;
      rK[K-1] = (1-freqK[K-2])*factor;
      for (i=1; i<K-1; i++)  rK[i] = (freqK[i]-freqK[i-1])*factor;
   }
   for (i=0; i<K; i++) freqK[i]=1.0/K;

   return (0);
}


double NetEvolvePlugin::zeroinA(double ax, double bx)//, double (*f)(double x))
{
  double a,b,c;				/* Abscissae	*/
  double fa;				/* f(a)			*/
  double fb;				/* f(b)			*/
  double fc;				/* f(c)			*/

  a = ax;  b = bx;  fa = epifunc(a);  fb = epifunc(b);
  c = a;   fc = fa;

  for(;;)		/* Main iteration loop	*/
  {
    double prev_step = b-a;	/* Distance from the last but one*/
							/* to the last approximation	*/
    double p;      			/* Interpolation step is calcu- */
    double q;      			/* lated in the form p/q; divi- */
  							/* sion operations is delayed   */
 							/* until the last moment		*/
    double new_step;      	/* Step at this iteration       */
   
    if( fabs(fc) < fabs(fb) )
    {                         		/* Swap data for b to be the 	*/
	a = b;  b = c;  c = a;          /* best approximation			*/
	fa=fb;  fb=fc;  fc=fa;
    }
    new_step = (c-b)/2;

    if( fabs(new_step) <= TOLERANCE || fb == (double)0 )
      return b;				/* Acceptable approx. is found	*/

    			/* Decide if the interpolation can be tried	*/
    if( fabs(prev_step) >= TOLERANCE	/* If prev_step was large enough*/
	&& fabs(fa) > fabs(fb) )	/* and was in true direction,	*/
    {					/* Interpolatiom may be tried	*/
	register double t1,cb,t2;
	cb = c-b;
	if( a==c )			/* If we have only two distinct	*/
	{				/* points linear interpolation 	*/
	  t1 = fb/fa;			/* can only be applied		*/
	  p = cb*t1;
	  q = 1.0 - t1;
 	}
	else				/* Quadric inverse interpolation*/
	{
	  q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
	  p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
	  q = (q-1.0) * (t1-1.0) * (t2-1.0);
	}
	if( p>(double)0 )	/* p was calculated with the op-*/
	  q = -q;			/* posite sign; make p positive	*/
	else				/* and assign possible minus to	*/
	  p = -p;			/* q				*/

	if( p < (0.75*cb*q-fabs(TOLERANCE*q)/2)	/* If b+p/q falls in [b,c]*/
	    && p < fabs(prev_step*q/2) )	/* and isn't too large	*/
	  new_step = p/q;					/* it is accepted	*/
    }

    if( fabs(new_step) < TOLERANCE )	/* Adjust the step to be not less*/
      if( new_step > (double)0 )	/* than tolerance		*/
	new_step = TOLERANCE;
      else
	new_step = -TOLERANCE;

    a = b;  fa = fb;			/* Save the previous approx.	*/ 
    b += new_step;  fb = epifunc(b);	/* Do step to a new approxim.	*/
    if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) )
    {                 			/* Adjust c for it to have a sign*/
      c = a;  fc = fa;                  /* opposite to that of b	*/
    }
  }

}
double NetEvolvePlugin::zeroinB(double ax, double bx)//, double (*f)(double x))
{
  double a,b,c;				/* Abscissae	*/
  double fa;				/* f(a)			*/
  double fb;				/* f(b)			*/
  double fc;				/* f(c)			*/

  a = ax;  b = bx;  fa = epiSubfunc(a);  fb = epiSubfunc(b);
  c = a;   fc = fa;

  for(;;)		/* Main iteration loop	*/
  {
    double prev_step = b-a;	/* Distance from the last but one*/
							/* to the last approximation	*/
    double p;      			/* Interpolation step is calcu- */
    double q;      			/* lated in the form p/q; divi- */
  							/* sion operations is delayed   */
 							/* until the last moment		*/
    double new_step;      	/* Step at this iteration       */
   
    if( fabs(fc) < fabs(fb) )
    {                         		/* Swap data for b to be the 	*/
	a = b;  b = c;  c = a;          /* best approximation			*/
	fa=fb;  fb=fc;  fc=fa;
    }
    new_step = (c-b)/2;

    if( fabs(new_step) <= TOLERANCE || fb == (double)0 )
      return b;				/* Acceptable approx. is found	*/

    			/* Decide if the interpolation can be tried	*/
    if( fabs(prev_step) >= TOLERANCE	/* If prev_step was large enough*/
	&& fabs(fa) > fabs(fb) )	/* and was in true direction,	*/
    {					/* Interpolatiom may be tried	*/
	register double t1,cb,t2;
	cb = c-b;
	if( a==c )			/* If we have only two distinct	*/
	{				/* points linear interpolation 	*/
	  t1 = fb/fa;			/* can only be applied		*/
	  p = cb*t1;
	  q = 1.0 - t1;
 	}
	else				/* Quadric inverse interpolation*/
	{
	  q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
	  p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
	  q = (q-1.0) * (t1-1.0) * (t2-1.0);
	}
	if( p>(double)0 )	/* p was calculated with the op-*/
	  q = -q;			/* posite sign; make p positive	*/
	else				/* and assign possible minus to	*/
	  p = -p;			/* q				*/

	if( p < (0.75*cb*q-fabs(TOLERANCE*q)/2)	/* If b+p/q falls in [b,c]*/
	    && p < fabs(prev_step*q/2) )	/* and isn't too large	*/
	  new_step = p/q;					/* it is accepted	*/
    }

    if( fabs(new_step) < TOLERANCE )	/* Adjust the step to be not less*/
      if( new_step > (double)0 )	/* than tolerance		*/
	new_step = TOLERANCE;
      else
	new_step = -TOLERANCE;

    a = b;  fa = fb;			/* Save the previous approx.	*/ 
    b += new_step;  fb = epiSubfunc(b);	/* Do step to a new approxim.	*/
    if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) )
    {                 			/* Adjust c for it to have a sign*/
      c = a;  fc = fa;                  /* opposite to that of b	*/
    }
  }

}

PluginProxy<NetEvolvePlugin> NetEvolvePluginProxy = PluginProxy<NetEvolvePlugin>("NetEvolve", PluginManager::getInstance());

//#undef TOLERANCE




//#undef BIG_NUMBER
//#undef BRACKET

//#undef BRACKET


//#undef BIG_NUMBER
