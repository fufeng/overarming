/**
Import an undirected network

*  (start)
8   (number of total nodes)
0	1 2 3 7 -1    (list of neighbor nodes with ascending order)
2	5 6 -1  (list of neighbor nodes without repetition with the lists before)
3	4 6 -1  (Line 3)
4 -1  (Line 5)
5	7 -1
-1	-1 (end)

==========

*/

#ifndef UndirectedGraph_h_
#define UndirectedGraph_h_

#include "Graph.h"

class UndirectedGraph : public Graph
{
public:
	UndirectedGraph() { IsDirected=false; }
	~UndirectedGraph() { }

	virtual long getTotalEdges(void) const; 

	virtual long getInDegree(long node) const; 
	virtual long getMaxInDegree(void) const; 
	virtual double getAverageInDegree(void) const; 
	virtual double getInDegreeVariance(void) const; 

	
	virtual void getDegreeCorrelation(double r[4]) const;
	
	double getClustering(long node) const; 

	virtual void outputDegrees(ostream &out) const; 
	virtual void outputOutDegreeDistribution(ofstream &out) const; 
	virtual void outputInDegreeDistribution(ofstream &out) const; 

	void outputDegreesAndClusterings(ostream &out) const;
	void	outputDegreeDistributionAndClusteringAverage(ofstream &out) const; 

	
	virtual int getNeighborType(long from, long to) const;


	virtual const long * getInNeighbors(long node, long &k) const;

	virtual bool isConnected(void) const; 
	virtual bool isStronglyConnected(void) const; 

	virtual bool addEdge(long from, long to); 
	virtual bool delEdge(long from, long to); 

	
	virtual double changeEdges(double p);

	virtual void input(istream &in); 


	virtual void output(ostream &out, long n1=-1, long n2=-1) const;


	virtual void outputToNwb(ostream &out) const;


	virtual void outputConnectedGraph(ostream & out, long node, long n) const;


	virtual void outputLargestConnectedGraph(ostream &out) const;


	void outputDirectedGraph(ostream &out) const;

private:
	bool changeEdges(long n1, long n2, long n3, long n4, long *chgTimes, bool check=false);
	
};


long UndirectedGraph::getTotalEdges(void) const
{
	if(!isValid()) return -100;
	long m=0;
	for(long i=0;i<N;i++) m+=OutDegree[i];
	return m/2; 
}


long UndirectedGraph::getInDegree(long node) const
{
	return getOutDegree(node);
}


long UndirectedGraph::getMaxInDegree(void) const
{
	return getMaxOutDegree();
}


double UndirectedGraph::getAverageInDegree(void) const
{
	return getAverageOutDegree();
}


double UndirectedGraph::getInDegreeVariance(void) const
{
	return getOutDegreeVariance();
}


void UndirectedGraph::getDegreeCorrelation(double r[4]) const
{
	if(!isValid())
	{
		r[0]=r[1]=r[2]=r[3]=-100;
		return;
	}
	
	double r1=0, r2=0, r3=0, k=0;
	long i, j;
	for(i=0;i<N;i++)
	{
		j=binSearch(OutNeighbors[i], OutDegree[i], i);
		j=-1-j;
		for(;j<OutDegree[i];j++) r1+=OutDegree[i]*OutDegree[OutNeighbors[i][j]];
		k=OutDegree[i]*OutDegree[i];
		r2+=k;
		r3+=k*OutDegree[i];
	}
	r2*=r2;
	double m=getTotalEdges();
	double num=4*m*r1-r2;
	double den=2*m*r3-r2;
	if(den>-Eps && den<Eps)
	{ 
		r[0]=r[1]=r[2]=r[3]=-1234;
		return;
	}
	else r[0]=r[1]=r[2]=r[3]=num/den;
}


double UndirectedGraph::getClustering(long node) const
{
	if(!isValid() || node<0 || node>=N) return -100;


	if(OutDegree[node]==0) return -2;
	if(OutDegree[node]==1) return -1;
	
	double clus=0;
	long i,j;
	for(i=0;i<OutDegree[node];i++)
	{
		for(j=i+1;j<OutDegree[node];j++)
		{
			if(isOutNeighbor(OutNeighbors[node][i], OutNeighbors[node][j])) clus+=1;
		}
	}
	return clus*2/OutDegree[node]/(OutDegree[node]-1);
}


void UndirectedGraph::outputDegrees(ostream &out) const
{
	outputDegreesAndClusterings(out);
}


void  UndirectedGraph::outputOutDegreeDistribution(ofstream &out) const
{
	outputDegreeDistributionAndClusteringAverage(out);
}


void  UndirectedGraph::outputInDegreeDistribution(ofstream &out) const
{
	outputDegreeDistributionAndClusteringAverage(out);
}


void UndirectedGraph::outputDegreesAndClusterings(ostream &out) const
{
	if(!isValid()) return;
	out<<"Undirected Degrees And Clusterings"<<endl;
	out<<"TotalNodes= "<<N<<" MaxDegree= "<<getMaxOutDegree()<<endl;
	out<<"Node\tDegree\tClustering"<<endl;
	out<<"*"<<endl;
	for(long i=0;i<N;i++) out<<i<<"\t"<<OutDegree[i]<<"\t"<<getClustering(i)<<endl;
	out<<"-1\t-1\t-1"<<endl;
}


void UndirectedGraph::outputDegreeDistributionAndClusteringAverage(ofstream &out) const
{
	if(!isValid()) return;
	long i, j, maxD=getMaxOutDegree()+1;
	double d;
	
	long * distr = new long[maxD]; 
	if(distr==NULL) handleError(__FILE__ , __LINE__);
	double * clus = new double[maxD]; 
	if(clus==NULL) handleError(__FILE__ , __LINE__);
	double * adn = new double[maxD]; 
	if(adn==NULL) handleError(__FILE__ , __LINE__);
	
	for(i=0;i<maxD;i++) { distr[i]=0; clus[i]=0; adn[i]=0; }
	for(i=0;i<N;i++)
	{
		distr[OutDegree[i]]++;
		clus[OutDegree[i]]+=getClustering(i);
		if(OutDegree[i]!=0)
		{
			d=OutDegree[OutNeighbors[i][0]];
			for(j=1;j<OutDegree[i];j++) d+=OutDegree[OutNeighbors[i][j]];
			adn[OutDegree[i]]+=d/OutDegree[i];
		}
	}

	out<<"Undirected Degree Distribution And Clustering Average"<<endl;
	out<<"ToalNodes= "<<N<<"\tMaxDegree= "<<maxD-1<<endl;
	out<<"Degree\tRatio\tAverageDegreeOfNeighbors\tClusteringAverage"<<endl;
	out<<"*"<<endl;
	for(i=0;i<maxD;i++)
	{
		if(distr[i]!=0)  out<<i<<"\t"<<distr[i]/(N+Zero)<<"\t"
			<<adn[i]/distr[i]<<"\t"<<clus[i]/distr[i]<<endl;
	}
	out<<"-1\t-1\t-1\t-1"<<endl;
}

int UndirectedGraph::getNeighborType(long from, long to) const
{
	if(!isValid()) return -100;
	if(from<0 || from>=N || to<0 || to>=N ||from==to) return -100;

	int i=0;
	if(binSearch(OutNeighbors[from], OutDegree[from], to)>=0) i+=3;

	return i;
}


const long * UndirectedGraph::getInNeighbors(long node, long &k) const
{
	return getOutNeighbors(node, k);
}


bool UndirectedGraph::isConnected(void) const
{
	if(!isValid()) return false;
	
	long nc, md, sd;
	getDistances(0, nc, md, sd);
	return nc==N;
}


bool UndirectedGraph::isStronglyConnected(void) const
{
	return isConnected();
}


bool UndirectedGraph::addEdge(long from, long to)
{
	if(!isValid() || from<0 || from>=N || to<0 || to>=N || from==to) return false;
	

	for(int t=0;t<2;t++)
	{
		long pos=binSearch(OutNeighbors[from], OutDegree[from], to);
		if(pos>=0) return false; 
		pos=-1-pos; 
		if(OutDegree[from]==OutCapacity[from])
		{ 
			OutCapacity[from]+=EnlargeStep;
			OutNeighbors[from]=enlarge(OutNeighbors[from], OutDegree[from], EnlargeStep, pos);
			if(OutNeighbors[from]==NULL) handleError(__FILE__ , __LINE__);
		}
		else
		{
			for(long i=OutDegree[from];i>pos;i--) OutNeighbors[from][i]=OutNeighbors[from][i-1];
		}
		
		OutNeighbors[from][pos]=to;
		OutDegree[from]++;

		swap(from, to); 
	}

	return true;
}


bool UndirectedGraph::delEdge(long from, long to)
{
	if(!isValid()) return false;
	if(from<0 || from>=N || to<0 || to>=N || from==to) return false;
	

	for(int t=0;t<2;t++)
	{
		long pos=binSearch(OutNeighbors[from], OutDegree[from], to);
		if(pos<0) return false; 
		OutDegree[from]--;
		for(long i=pos; i<OutDegree[from]; i++) OutNeighbors[from][i]=OutNeighbors[from][i+1];

		swap(from, to); 
	}

	return true;
}


double UndirectedGraph::changeEdges(double p)
{
	if(!isValid()) return -100;
	if(p<=0 || p>1) return -1;

	long toChg=long(p*getTotalEdges()+1); 
	long hasChg=0; 
	long p1, p3, n1, n2, n3, n4; 
	long * chgTimes = new long [N]; 
	if(chgTimes==NULL) handleError(__FILE__ , __LINE__);
	
	long i;
	for(i=0;i<N;i++) chgTimes[i]=0;

	int v;

	for(i=2*toChg; i>0 && hasChg<toChg; i--)
	{ 

		n1=getRand(N);
	 	if(chgTimes[n1]==OutDegree[n1]) continue;
		p1=getRand(chgTimes[n1], OutDegree[n1]);
		n2=OutNeighbors[n1][p1];

		n3=getRand(N);
		if(n3==n1 || n3==n2 || chgTimes[n3]==OutDegree[n3]) continue;
		p3=getRand(chgTimes[n3], OutDegree[n3]);
		n4=OutNeighbors[n3][p3];
		if(n4==n1 || n4==n2) continue;


		swap(OutNeighbors[n1][chgTimes[n1]], OutNeighbors[n1][p1]);
		p1=lineSearch(OutNeighbors[n2], OutDegree[n2], n1);
		swap(OutNeighbors[n2][chgTimes[n2]], OutNeighbors[n2][p1]);
		swap(OutNeighbors[n3][chgTimes[n3]], OutNeighbors[n3][p3]);
		p3=lineSearch(OutNeighbors[n4], OutDegree[n4], n3);
		swap(OutNeighbors[n4][chgTimes[n4]], OutNeighbors[n4][p3]);


		v=0;
		if(lineSearch(OutNeighbors[n1], OutDegree[n1], n3)>=0) v+=1;
		if(lineSearch(OutNeighbors[n2], OutDegree[n2], n4)>=0) v+=2;
		if(lineSearch(OutNeighbors[n1], OutDegree[n1], n4)>=0) v+=4;
		if(lineSearch(OutNeighbors[n2], OutDegree[n2], n3)>=0) v+=8;

		if(v==0)
		{ 
			if(!changeEdges(n1, n2, n3, n4, chgTimes, true))
				changeEdges(n1, n2, n4, n3, chgTimes);
		}

		else if(v<4) changeEdges(n1, n2, n4, n3, chgTimes);

		else if(v%4==0) changeEdges(n1, n2, n3, n4, chgTimes);
		else continue; 
		hasChg+=2;
	}


	long *b = new long[N];
	if(b==NULL) handleError(__FILE__ , __LINE__);
	long *c = new long[N];
	if(c==NULL) handleError(__FILE__ , __LINE__);
	for(i=0;i<N;i++)  countSort(OutNeighbors[i], OutDegree[i], b, c, N);

	delete[] b;
	delete[] c;
	delete[] chgTimes;

	return hasChg/(toChg+Zero);
}


bool UndirectedGraph::changeEdges(long n1, long n2, long n3, long n4, long *chgTimes, bool check)
{
	OutNeighbors[n1][chgTimes[n1]++]=n3;
	OutNeighbors[n2][chgTimes[n2]++]=n4;
	OutNeighbors[n3][chgTimes[n3]++]=n1;
	OutNeighbors[n4][chgTimes[n4]++]=n2;
	if(check)
	{ 
		if(!isConnected())
		{
			OutNeighbors[n1][--chgTimes[n1]]=n2;
			OutNeighbors[n2][--chgTimes[n2]]=n1;
			OutNeighbors[n3][--chgTimes[n3]]=n4;
			OutNeighbors[n4][--chgTimes[n4]]=n3;
			return false;
		}
	}
	
	return true;
	
}


void UndirectedGraph::input(istream &in)
{
	if(in==NULL) handleError(__FILE__ , __LINE__);
	
	clean();
	
	char c='a';
	while(!in.eof() && c!='*') in>>c;
	if(in.eof()) handleError(__FILE__ , __LINE__);
	
	in>>N;
 	if(N<MinNodes) handleError(__FILE__ , __LINE__);

	OutNeighbors=new long * [N];
	if(OutNeighbors==NULL) handleError(__FILE__ , __LINE__);

	OutCapacity=new long[N];
	if(OutCapacity==NULL) handleError(__FILE__ , __LINE__);

	OutDegree=new long[N];
	if(OutDegree==NULL) handleError(__FILE__ , __LINE__);

	InDegree=new long[N]; // 无向图中入度只分配空间，不使用
	if(InDegree==NULL) handleError(__FILE__ , __LINE__);

	long i, j, n1, n2, t;
	
	for(i=0;i<N;i++)
	{
		OutDegree[i]=0;
		OutCapacity[i]=3*EnlargeStep;
		OutNeighbors[i]=new long[OutCapacity[i]];
		if(OutNeighbors[i]==NULL) handleError(__FILE__ , __LINE__);
	}

	if(in.eof()) handleError(__FILE__ , __LINE__);
	in>>i;
	while(i>=0)
	{
		if(i>=N) handleError(__FILE__, __LINE__, "DirectedGraph::input()", i);
		
		if(in.eof()) handleError(__FILE__ , __LINE__);
		in>>j;
		while(j>=0)
		{
			if(j>=N) handleError(__FILE__, __LINE__, "DirectedGraph::input()", j);
			
			// 文件中的节点未按从小到大排列
			if(j<=i) handleError(__FILE__, __LINE__, "DirectedGraph::input()", i, j);
			
			n1=i; n2=j; // 先将n2添加到n1邻居表中，再将n1添加到n2邻居表中
			for(t=0;t<2;t++)
			{
				// 文件中的邻居未按从小到大排列
				if(OutDegree[n1]>0 && n2<=OutNeighbors[n1][OutDegree[n1]-1])
					handleError(__FILE__, __LINE__, "DirectedGraph::input()", i, j);

				if(OutDegree[n1]==OutCapacity[n1])
				{ // 需要扩充n1的邻居表
					OutCapacity[n1]+=EnlargeStep;
					OutNeighbors[n1]=enlarge(OutNeighbors[n1], OutDegree[n1], EnlargeStep);
				}
				OutNeighbors[n1][OutDegree[n1]++]=n2;

				swap(n1, n2);
			}
			
			if(in.eof()) handleError(__FILE__ , __LINE__);
			in>>j;
		}

		if(in.eof()) handleError(__FILE__ , __LINE__);
		in>>i;
	}

	//outputEvent("UndirectedGraph::input() -- OK");
}


void UndirectedGraph::output(ostream &out, long n1, long n2) const
{
	if(!isValid()) return;
	if(out==NULL) handleError(__FILE__ , __LINE__);
	
	if(n1<0) n1=0;
	if(n2<MinNodes || n2>N) n2=N;
	if(n2-n1<MinNodes) return;
	
	out<<"Undirected Graph"<<endl;
	out<<"*"<<endl<<N<<endl;

	long i, j;
	for(i=n1;i<n2;i++)
	{
		j=binSearch(OutNeighbors[i], OutDegree[i], i);
		j=-1-j; 
		if(j<OutDegree[i])
		{
			out<<i<<"\t";
			for(;j<OutDegree[i];j++)
				if(OutNeighbors[i][j]>=n1 && OutNeighbors[i][j]<n2) out<<OutNeighbors[i][j]<<" ";
			out<<"-1"<<endl;
		}
	}

	out<<"-1\t-1"<<endl;
}


void UndirectedGraph::outputToNwb(ostream &out) const
{
	if(!isValid()) return;
	
	out<<"*Nodes "<<N<<endl;
	out<<"*UndirectedEdges"<<endl;

	long i, j;
	for(i=0;i<N;i++)
	{
		j=binSearch(OutNeighbors[i], OutDegree[i], i);
		j=-1-j;
		if(j<OutDegree[i])
		{
		
			for(;j<OutDegree[i];j++) out<<i+1<<" "<<OutNeighbors[i][j]+1<<endl;
		}
	}
}



void UndirectedGraph::outputConnectedGraph(ostream & out, long node, long n) const
{
	if(!isValid() || n<MinNodes|| node<0 || node>=N)  return;
	if(OutDegree[node]==0)
	{
		cout<<"UndirectedGraph::outputConnectedGraph() -- Node "<<node<<" is isolated"<<endl;
		return;
	}

	long top, m, n1, n2, i ,j;
	long * chosenNodes = new long [N];
	if(chosenNodes==NULL) handleError(__FILE__ , __LINE__);
	bool * isVisited = new bool [N];
	if(isVisited==NULL) handleError(__FILE__ , __LINE__);
	for(i=0;i<N;i++) isVisited[i]=false;

	chosenNodes[0]=node; isVisited[node]=true;
	top=0; m=1;
	for(top=0;top<m;top++)
	{
		n1=chosenNodes[top];
		for(j=0;j<OutDegree[n1];j++)
		{
			n2=OutNeighbors[n1][j];
			if(!isVisited[n2])
			{
				isVisited[n2]=true;
				chosenNodes[m++]=n2;
				if(m==n)
				{ 
					top=m; break;
				}
			}
		}
	}


	long *newID=chosenNodes; 
	top=0;
	for(i=0;i<N;i++)
	{ 
		if(isVisited[i]) newID[i]=top++;
		else newID[i]=N;
	}
	
	out<<"Undirected Graph"<<endl;
	out<<"Generated by UndirectedGraph::outputConnectedGraph()"<<endl;
	out<<"*"<<endl<<m<<endl;

	for(i=0;i<N;i++)
	{
		if(newID[i]<N)
		{ 
			j=binSearch(OutNeighbors[i], OutDegree[i], i);
			j=-1-j; 
			if(j<OutDegree[i])
			{
				out<<newID[i]<<"\t";
				for(;j<OutDegree[i];j++)
				{
					if(newID[OutNeighbors[i][j]]<N)  out<<newID[OutNeighbors[i][j]]<<" ";
				}
				out<<"-1"<<endl;
			}
		}
	}

	out<<"-1\t-1"<<endl;

	delete[] chosenNodes;
	delete[] isVisited;

}


void UndirectedGraph::outputLargestConnectedGraph(ostream & out) const
{
	if(!isValid()) return;
	
	long nc, md, sd, i, j,
		theOne, 
		maxConnected, 
		leftNodes, 
		v; 

	bool *isVisited = new bool[N]; 
	for(i=0;i<N;i++) isVisited[i]=false;

	leftNodes=N;


	for(v=0;v<N;v++)
	{
		if(OutDegree[v]>0) break;
		else { isVisited[v]=true; leftNodes--; }
	}
	if(v==N)
	{
		cout<<"UndirectedGraph::outputLargestConnectedGraph() -- Every node is isolated"<<endl;
		delete[] isVisited;
		return;
	}

	long *dist1 = new long[N]; 
	if(dist1==NULL) handleError(__FILE__ , __LINE__);

	getDistances(v, nc, md, sd, dist1);
	isVisited[v]=true; leftNodes-=nc;
	if(nc==N)
	{ 
		delete[] dist1;
		delete[] isVisited;
		cout<<"UndirectedGraph::outputLargestConnectedGraph() -- The Graph itself is connected"<<endl;
		output(out);
		return;
	}
	for(i=v+1;i<N;i++) if(dist1[i]<N) isVisited[i]=true;
	theOne=v; maxConnected=nc;
	for(v++; v<N && isVisited[v]; v++) {}

	long *dist2 = new long[N]; 
	if(dist2==NULL) handleError(__FILE__ , __LINE__);

	while(v<N && leftNodes>maxConnected)
	{
		getDistances(v, nc, md, sd, dist2);
		isVisited[v]=true; leftNodes-=nc;
		
		if(nc>1)
		{
			for(i=v+1;i<N;i++) if(dist2[i]<N) isVisited[i]=true;
			if(nc>maxConnected)
			{ 
				long * temp=dist1;
				dist1=dist2; dist2=temp;
				theOne=v; maxConnected=nc;
			}
		}

		for(v++; v<N && isVisited[v]; v++) {}
	}
	

	long * newID=dist1; 
	nc=0;
	for(i=0; i<N; i++)
	{
		if(dist1[i]<N) newID[i]=nc++;
	}

	out<<"Undirected Graph"<<endl;
	out<<"Generated by UndirectedGraph::outputLargestConnectedGraph()"<<endl;
	out<<"*"<<endl<<maxConnected<<endl;

	for(i=0;i<N;i++)
	{
		if(newID[i]<N)
		{ 
			j=binSearch(OutNeighbors[i], OutDegree[i], i);
			j=-1-j; 
			if(j<OutDegree[i])
			{
				out<<newID[i]<<"\t";
				for(;j<OutDegree[i];j++)
				{ 
					out<<newID[OutNeighbors[i][j]]<<" ";
				}
				out<<"-1"<<endl;
			}
		}
	}

	out<<"-1\t-1"<<endl;

	delete[] dist1;
	delete[] dist2;
	delete[] isVisited;
	
}


void UndirectedGraph::outputDirectedGraph(ostream &out) const
{
	if(!isValid()) return;
	if(out==NULL) handleError(__FILE__ , __LINE__);
	
	out<<"Directed Graph"<<endl;
	out<<"Generated by UndirectedGraph::outputDirectedGraph()"<<endl;
	out<<"*"<<endl<<N<<endl;

	long i, j;
	for(i=0;i<N;i++)
	{
		out<<i<<"\t";
		for(j=0;j<OutDegree[i];j++)
			out<<OutNeighbors[i][j]<<" ";
		out<<"-1"<<endl;
	}

	out<<"-1\t-1"<<endl;
}

#endif
