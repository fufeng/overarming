
#ifndef Graph_h_
#define Graph_h_

#include "PublicFunction.h"

const long MinNodes=3; 
const long EnlargeStep=5 ; 

class Graph
{
public:
	Graph();
	~Graph();

	void clean(void); 
	inline bool isValid(void) const; 
	inline bool isDirected(void) const; 
	
	inline long getTotalNodes(void) const; 
	virtual long getTotalEdges(void) const =0; 

	long getOutDegree(long node) const; 
	long getMaxOutDegree(void) const; 
	double getAverageOutDegree(void) const; 
	double getOutDegreeVariance(void) const; 

	virtual long getInDegree(long node) const =0; 
	virtual long getMaxInDegree(void) const =0; 
	virtual double getAverageInDegree(void) const =0; 
	virtual double getInDegreeVariance(void) const =0; 


	virtual void getDegreeCorrelation(double r[4]) const =0; // 
	
	virtual void outputDegrees(ostream &out) const =0; 
	virtual void outputOutDegreeDistribution(ofstream &out) const =0; 
	virtual void outputInDegreeDistribution(ofstream &out) const =0; 


	virtual int getNeighborType(long from, long to) const =0;
	
	bool isOutNeighbor(long from, long to) const; 

	const long * getOutNeighbors(long node, long &k) const;
	long getRandomOutNeighbor(long node); 

	virtual const long * getInNeighbors(long node, long &k) const =0;


	void getDistances(long node, long &nc, long &md, long &sd, long *dist=NULL)  const;


	long getDistance(long from, long to)  const;
	

	void outputDistances(ostream &out, double &avgD, long &maxD) const;
	
	virtual bool isConnected(void) const =0; 
	virtual bool isStronglyConnected(void) const =0; 

	virtual bool addEdge(long from, long to) =0; 
	virtual bool delEdge(long from, long to) =0; 


	virtual double changeEdges(double p) =0;

	virtual void input(istream &in)=0; 


	virtual void output(ostream &out, long n1=-1, long n2=-1) const=0;
	
	friend istream& operator>>(istream &in, Graph &graph);
	friend ostream& operator<<(ostream &out, const Graph &graph);


	virtual void outputToNwb(ostream &out) const =0;


	virtual void outputConnectedGraph(ostream & out, long node, long n) const =0;


	virtual void outputLargestConnectedGraph(ostream &out) const =0;

protected:
	long N; 
	bool IsDirected; 


	long * OutDegree;
	long ** OutNeighbors; // OutNeighbors[i][0]<OutNeighbors[i][1]<...<OutNeighbors[i][OutDegree[i]-1]
	long * OutCapacity;

	long * InDegree; 

};

Graph::Graph()
{
	N=-1;
	OutNeighbors=NULL;
	OutCapacity=NULL;
	OutDegree=NULL;
	InDegree=NULL;
}

Graph::~Graph()
{
	clean();
}


void Graph::clean(void)
{
	if(N<MinNodes) return;

	for(long i=0;i<N;i++)
	{
		delete[] OutNeighbors[i];   OutNeighbors[i]=NULL;
	}

	delete[] OutNeighbors;  OutNeighbors=NULL;
	delete[] OutCapacity;  OutCapacity=NULL;
	delete[] OutDegree;  OutDegree=NULL;
	delete[] InDegree;  InDegree=NULL;

	N=-1;
}


inline bool Graph::isValid(void) const
{
	return N>=MinNodes;
}


inline bool Graph::isDirected(void) const
{
	return IsDirected;
}


inline long Graph::getTotalNodes(void) const
{
	return N;
}


long Graph::getOutDegree(long node) const
{
	if(!isValid() || node<0 || node>=N) return -100;
	else return OutDegree[node];
}


long Graph::getMaxOutDegree(void) const
{
	if(!isValid()) return -100;
	long m=-1;
	for(long i=0; i<N; i++) if(OutDegree[i]>m) m=OutDegree[i];
	return m;
}


double Graph::getAverageOutDegree(void) const
{
	if(!isValid()) return -100;
	double degSum=0;
	for(long i=0; i<N; i++) degSum+=OutDegree[i];
	return degSum/N;
}


double Graph::getOutDegreeVariance(void) const
{
	if(!isValid()) return -100;
	double ad=getAverageOutDegree();
	double d2=0;
	for(long i=0;i<N;i++) d2+=(OutDegree[i]-ad)*(OutDegree[i]-ad);
	return d2/N;
}

bool Graph::isOutNeighbor(long from, long to) const
{
	if(!isValid()) return false;
	if(from<0 || from>=N || to<0 || to>=N ||from==to) return false;
	return binSearch(OutNeighbors[from], OutDegree[from], to)>=0;
}


const long * Graph::getOutNeighbors(long node, long & k) const
{
	if(!isValid() || node<0 || node>=N) { k=-100; 	return NULL; 	}
	
	k=OutDegree[node];
	return OutNeighbors[node];
}


long Graph::getRandomOutNeighbor(long node)
{
	if(!isValid() || node<0 || node>=N) return -100;
	if(OutDegree[node]==0) return -1;
	return OutNeighbors[node][getRand(OutDegree[node])];
}


void Graph::getDistances(long node, long & nc, long & md, long & sd, long *dist) const
{
	if(!isValid() || node<0 || node>=N)
	{ 
		nc=md=sd= -100;
		return;
	}

	if(OutDegree[node]==0)
	{ 
		md=sd=0;
		nc=1;
		return;
	}

	bool r=true;
	if(dist==NULL)
	{ 
		r=false;
		dist=new long[N];
		if(dist==NULL) handleError(__FILE__ , __LINE__);
	}

	long *list1=new long[N];
	if(list1==NULL) handleError(__FILE__ , __LINE__);

	long *list2=new long[N];
	if(list2==NULL) handleError(__FILE__ , __LINE__);

	long i1, i2, n1, n2, len1, len2, curDist;
	long *list3;
	
	for(i1=0;i1<N;i1++) dist[i1]=N; 
	dist[node]=curDist=0;
	list1[0]=node; len1=1; len2=0;
	nc=1; sd=0;

	while(len1>0)
	{
		curDist++;
		for(i1=0;i1<len1;i1++)
		{
			n1=list1[i1];
			for(i2=0;i2<OutDegree[n1];i2++)
			{ 
				n2=OutNeighbors[n1][i2];
				if(dist[n2]==N)
				{ 
					list2[len2++]=n2;
					dist[n2]=curDist;
					sd+=curDist;
					nc++;
				}
			}
		}

		list3=list1;
		list1=list2; len1=len2;
		list2=list3; len2=0;
	}

	md=curDist-1;

	delete[] list1;
	delete[] list2;
	if(!r) { delete[] dist; dist=NULL; }
}


long Graph::getDistance(long from, long to)  const
{
	if(!isValid() || from<0 || from>=N || to<0 || to>=N)
	{ 
		return N;
	}

	if(OutDegree[from]==0) return N; 
	if(from==to) return 0; 

	bool *isVisited=new bool[N];
	if(isVisited==NULL) handleError(__FILE__ , __LINE__);

	long *list1=new long[N];
	if(list1==NULL) handleError(__FILE__ , __LINE__);

	long *list2=new long[N];
	if(list2==NULL) handleError(__FILE__ , __LINE__);

	long i1, i2, n1, n2, len1, len2, curDist;
	long *list3;
	
	for(i1=0;i1<N;i1++) isVisited[i1]=false; 
	isVisited[from]=true;
	list1[0]=from; len1=1; len2=0;
	curDist=0;

	while(len1>0)
	{
		curDist++;
		for(i1=0;i1<len1;i1++)
		{
			n1=list1[i1];
			for(i2=0;i2<OutDegree[n1];i2++)
			{ 
				n2=OutNeighbors[n1][i2];
				if(!isVisited[n2])
				{ 
					list2[len2++]=n2;
					isVisited[n2]=true;

					if(n2==to)
					{
						delete[] list1;
						delete[] list2;
						delete[] isVisited;
						return curDist;
					}
				}
			}
		}

		list3=list1;
		list1=list2; len1=len2;
		list2=list3; len2=0;
	}

	delete[] list1;
	delete[] list2;
	delete[] isVisited;
	return N;
}

void Graph::outputDistances(ostream &out, double &avgD, long &maxD) const
{
	if(!isValid() || out==NULL) { avgD=maxD=0; return; }
	
	bool *isVisited=new bool[N];
	if(isVisited==NULL) handleError(__FILE__ , __LINE__);

	long *list1=new long[N];
	if(list1==NULL) handleError(__FILE__ , __LINE__);

	long *list2=new long[N];
	if(list2==NULL) handleError(__FILE__ , __LINE__);

	double *distD=new double[N+1]; // ¾àÀëµÄ·Ö²¼
	if(distD==NULL) handleError(__FILE__ , __LINE__);

	long i1, i2, len1, len2, curDist, node, nc;
	long *list3;

	for(i1=0;i1<N;i1++) distD[i1]=0;

	for(node=0;node<N;node++)
	{
		for(i1=0;i1<N;i1++) isVisited[i1]=false; 
		isVisited[node]=true;
		list1[0]=node; len1=1; len2=0;
		curDist=0; nc=1;

		while(len1>0)
		{
			curDist++;
			for(i1=0;i1<len1;i1++)
			{
				for(i2=0;i2<OutDegree[list1[i1]];i2++)
				{ 
					if(!isVisited[OutNeighbors[list1[i1]][i2]])
					{ 
						list2[len2++]=OutNeighbors[list1[i1]][i2];
						isVisited[OutNeighbors[list1[i1]][i2]]=true;
					}
				}
			}

			distD[curDist]+=len2;
			nc+=len2; if(nc==N) break;
			list3=list1;
			list1=list2; len1=len2;
			list2=list3; len2=0;
		}

		if(node%10000==0) outputEvent("outputDistances() ", node);
	}

	avgD=0;
	maxD=0;
	double distN=0;
	for(i1=0;i1<N;i1++)
	{
		if(distD[i1]>0)
		{
			distN+=distD[i1];
			maxD=i1;
			avgD+=i1*distD[i1];
		}
	}
	avgD/=distN;
	
	out<<"generated by Graph::outputDistances()"<<endl;
	out<<"MaxDistance= "<<maxD<<"\t\t"<<"AverageDistance= "<<avgD<<endl;
	out<<"TotalPairs(Directed)= "<<distN<<endl;
	out<<"Distance\tRatio"<<endl;
	out<<"*"<<endl;
	for(i1=0;i1<N;i1++)
	{
		if(distD[i1]>0) out<<i1<<"\t"<<distD[i1]/distN<<endl;
	}
	out<<"-1\t-1"<<endl;

	delete[] isVisited;
	delete[] list1;
	delete[] list2;
	delete[] distD;

}

istream& operator>>(istream &in, Graph &graph)
{
	graph.input(in);
	return in;
}

ostream& operator<<(ostream &out, const Graph &graph)
{
	graph.output(out);
	return out;
}

#endif

