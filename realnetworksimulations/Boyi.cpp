
/*
Revised and updated by Feng Fu 
for gun control project
17 Jan 2023
*/

/**
1. prepare the network data file that will be imported 
2. determine model paramaters to simulate and compile the code and execute.
**/



#include "UndirectedGraph.h"
#include <math.h>

/** payoff matrix for gun adoptions */
const char C=0;
const char D=1;
double CD[2][2] = {
	-1 , 0.6 , //  -\delta_s  \delta_e
	-0.6 , -0.1  //   -\delta_e  -\delta_n     
};

double pa = 0.01; //provocation rate
double gs = -0.05; //solo payoff for owning a gun


////////// simulateon update methods //////////
#define Syn_Yes // synchronous update(Yes) vs asynchronous (No)
#define Method_2 // strategy update methods(1, 2, 3, 4): replicator, Fermi, 
double Method2_k = 1; // the temperature parameter of method 2 Fermi function
long T0=2000, T1=3000; // average window
long InitialTimes=2000; // number of initializations
#define OutEveryTime_No // output time series(Yes/No)
////////// simulateon update methods  //////////

class Boyi : public UndirectedGraph
{
public:
	Boyi(istream &in);
	~Boyi();

	/** inititialization */
	void initialize(double pc);


	inline void evolve(void);


	inline long getTotalCNodes(void) const { return TotalCNodes; }
	
private:
	char *State0, *State1; // status
	double *Payoff, *AveragePayoff; // payoff

	long TotalCNodes;

	/** synchronous update vs asynchronous update*/
	void evolveSynYes(void);
	void evolveSynNo(void);

	/** strategy update rules */
	inline char getNewState(long node);
	char getNewState1(long node);
	char getNewState2(long node);
	char getNewState3(long node);
	char getNewState4(long node);

	inline void calPayoffs(void); // caclulate Payoff
	inline void calAveragePayoffs(void); // caculate Average Payoff
};

Boyi::Boyi(istream &in)
{
	if(in==NULL) handleError(__FILE__, __LINE__, "in==NULL");
	input(in);

	// unconnected graphs
	if(!isConnected()) handleError(__FILE__, __LINE__, "The Graph is Not Connected!");

	State0=new char[N];
	if(State0==NULL) handleError(__FILE__, __LINE__, "State0==NULL");
	State1=new char[N];
	if(State1==NULL) handleError(__FILE__, __LINE__, "State1==NULL");

	Payoff=new double[N];
	if(Payoff==NULL) handleError(__FILE__, __LINE__, "Payoff==NULL");

	AveragePayoff=new double[N];
	if(AveragePayoff==NULL) handleError(__FILE__, __LINE__, "AveragePayoff==NULL");
}

Boyi::~Boyi()
{
	if(N<MinNodes) return;

	delete[] State0; State0=NULL;
	delete[] State1; State1=NULL;
	delete[] Payoff; Payoff=NULL;
	delete[] AveragePayoff; AveragePayoff=NULL;

}

/** Inititialization */
void Boyi::initialize(double pc)
{
	long * node=new long [N];
	if(node==NULL) handleError(__FILE__, __LINE__, "node==NULL");
	long i, nC2;
	for(i=0;i<N;i++) {
		State0[i]=C;
		node[i]=i;
	}
	TotalCNodes=N;
	nC2=long(N*pc);
	if(nC2>N) nC2=N;
	else if(nC2<0) nC2=0;
	
	while(TotalCNodes!=nC2) {
		i=getRand(TotalCNodes);
		State0[node[i]]=D;
		node[i]=node[--TotalCNodes];
	}

	calPayoffs();
	
	#ifdef Method_3
	calAveragePayoffs();
	#endif
	
	delete[] node;
}

/**evolutionary steps**/
inline void Boyi::evolve(void)
{
	#ifdef Syn_Yes
	evolveSynYes();
	#else
	evolveSynNo();
	#endif

	#ifdef Method_3
	calAveragePayoffs();
	#endif
}

/** update new states of nodes */
inline char Boyi::getNewState(long node)
{
	#ifdef Method_1
	return getNewState1(node);
	#else
		#ifdef Method_2
		return getNewState2(node);
		#else
			#ifdef Method_3
			return getNewState3(node);
			#else
				#ifdef Method_4
				return getNewState4(node);
				#else
					handleError(__FILE__, __LINE__, "Unknown getNewState Method");
					return -1;
				#endif  // 4
			#endif // 3
		#endif // 2
	#endif // 1
}

/** synchronous updating */
void Boyi::evolveSynYes(void)
{
	long node;


	for(node=0;node<N;node++)  State1[node]=getNewState(node);
	char *sTemp=State0; State0=State1; State1=sTemp;

	TotalCNodes=0;
	for(node=0;node<N;node++)  if(State0[node]==C) TotalCNodes++;

	calPayoffs();
}

/** asynchronous updating */
void Boyi::evolveSynNo(void)
{
	long n1=getRand(N), n2;
	char s0=State0[n1];
	char s1=getNewState(n1);
	if(s0!=s1) {
		State0[n1]=s1;
		TotalCNodes += s0-s1; // C=0, D=1
		long i;
		for(i=0;i<OutDegree[n1];i++) {
			n2=OutNeighbors[n1][i];
			Payoff[n1] += CD[s1][State0[n2]]-CD[s0][State0[n2]];
			Payoff[n2] += CD[State0[n2]][s1]-CD[State0[n2]][s0];
		}
	}
}

char Boyi::getNewState1(long node)
{
	long n1=node, n2, k1, k2;
	k1=OutDegree[n1];
	if(k1<=0) return State0[n1];
	
	n2=getRandomOutNeighbor(n1);
	k2=OutDegree[n2];
	if(k2<=0) return State0[n2];

	double p=(Payoff[n2]/k2-Payoff[n1]/k1)/(CD[1][0]-CD[0][1]);
	if(getRand01()<=p) return State0[n2];
	else return State0[n1];
}

char Boyi::getNewState2(long node)
{
	long n1=node, n2;
	n2=getRandomOutNeighbor(n1);
	if(n2<0) return State0[n1];

	double p=1/(1+exp((Payoff[n1]-Payoff[n2])/Method2_k));
	if(getRand01()<=p) return State0[n2];
	else return State0[n1];
}

char Boyi::getNewState3(long node)
{
	long n2=node, i;
	for(i=0;i<OutDegree[node];i++)
		if(AveragePayoff[OutNeighbors[node][i]]>AveragePayoff[n2]) n2=OutNeighbors[node][i];
	return State0[n2];
}

char Boyi::getNewState4(long node)
{
	long n1=node, n2, km;
	
	n2=getRandomOutNeighbor(n1);
	if(n2<0) return State0[n1];
	
	if(OutDegree[n1]>OutDegree[n2]) km=OutDegree[n1];
	else km=OutDegree[n2];

	double p=(Payoff[n2]-Payoff[n1])/km/(CD[1][0]-CD[0][1]);
	if(getRand01()<=p) return State0[n2];
	else return State0[n1];
}

/**  calculate Payoff */
inline void Boyi::calPayoffs(void)
{
	long n1, i;
	for(n1=0;n1<N;n1++)
	{
		Payoff[n1]=0;
		for(i=0;i<OutDegree[n1];i++)
			Payoff[n1]+=pa*CD[State0[n1]][State0[OutNeighbors[n1][i]]];
		Payoff[n1] = Payoff[n1]/OutDegree[n1];
		
		if(State0[n1]==C) Payoff[n1]+=gs;
	}
		

}

/**  calculate Average Payoff */
inline void Boyi::calAveragePayoffs(void)
{
	for(long node=0;node<N;node++)
		AveragePayoff[node]=Payoff[node]/OutDegree[node];
}

int main(int argc, char ** argv)
{
	time_t timeNow;
	srand((unsigned) time(&timeNow));

	outputEvent("Begin");

	const char *graphFnType="graph%02d.txt";
	char graphFn[15]; // graph00.txt, graph01.txt ...
	const char *rowTFnType="rowT%02d_b%03d_it%02d.txt";
	char rowTFn[25]; // rowT00_b100_it00.txt
	const char *rowBFnType="rowB%02d.txt";
	char rowBFn[15]; // rowB00.txt, rowB01.txt ...
	
	ifstream inGraph;
	ofstream outRowT, outRowBGraph, outRowBAll;

	Boyi *boyi;
	long N;

	double averageC, averageCInitial;
	long it, i, t, gn, gmax;
	double b;

	cout<<"please prepare graph00.txt, graph01.txt..."<<endl;
	cout<<"input the last index of the graphs (0~99): "; cin>>gmax;
	if(gmax<0 || gmax>99) handleError(__FILE__, __LINE__, "invalid graph index");
	gmax++;

	outRowBAll.open("rowBAll.txt");
	if(outRowBAll==NULL) handleError(__FILE__, __LINE__, "rowBAll.txt");
	outRowBAll<<"graph\tb\trow"<<endl;

	for(gn=0;gn<gmax;gn++)
	{
		sprintf(graphFn, graphFnType,gn);
		inGraph.open(graphFn);
		boyi=new Boyi(inGraph);
		inGraph.close();
		N=boyi->getTotalNodes();

		sprintf(rowBFn, rowBFnType, gn);
		outRowBGraph.open(rowBFn);
		if(outRowBGraph==NULL) handleError(__FILE__, __LINE__, rowBFn);
		outRowBGraph<<"b\tit\trow"<<endl;
		
		for(pa=0;pa<0.4;pa+=0.02) ////////// range of model parmaters to be simulated //////////
		{ 
			averageCInitial=0;

			for(it=0;it<InitialTimes;it++)
			{ 

				#ifdef OutEveryTime_Yes
				{
					sprintf(rowTFn, rowTFnType, gn, (int)(100*b),it);
					outRowT.open(rowTFn);
					if(outRowT==NULL) handleError(__FILE__, __LINE__, rowTFn);
					outRowT<<"time\trow"<<endl; 
				}
				#endif

				boyi->initialize(0.5); 

				#ifdef OutEveryTime_Yes
				{
					outRowT<<"0\t"<<boyi->getTotalCNodes()/(N+Zero)<<endl;
				}
				#endif
						
				i=T0+1;
				for(t=1;t<i;t++) {
					boyi->evolve();
					#ifdef OutEveryTime_Yes
					{
						outRowT<<t<<"\t"<<boyi->getTotalCNodes()/(N+Zero)<<endl;
					}
					#endif
				}

				averageC=0;
				i=T0+T1+1;
				for(t=T0+1;t<i;t++)
				{
					boyi->evolve();
					averageC+=boyi->getTotalCNodes();
					#ifdef OutEveryTime_Yes
					{
						outRowT<<t<<"\t"<<boyi->getTotalCNodes()/(N+Zero)<<endl;
					}
					#endif
				}

				#ifdef OutEveryTime_Yes
				{
					outRowT.close();
				}
				#endif

				// output data file
				outRowBGraph<<b<<"\t"<<it<<"\t"<<averageC/N/T1<<endl;
				averageCInitial+=averageC;
			} 

			// output data file
			outRowBGraph<<pa<<"\t-1\t"<<averageCInitial/N/T1/InitialTimes<<endl;
			outRowBAll<<gn<<"\t"<<pa<<"\t"<<averageCInitial/N/T1/InitialTimes<<endl;
			outputEvent("After pa=", pa);
		} 

		delete boyi;
		outRowBGraph.close();
		
		outputEvent("++ After Graph ", gn);
	} 

	outRowBAll.close();
	return 0;
}
