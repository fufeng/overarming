
#ifndef PublicFunction_h_
#define PublicFunction_h_

#include<iostream>
#include<fstream>
#include<stdlib.h>
using namespace std;

const double Zero=0;
const double Eps=1e-6; // Îó²î·¶Î§


void outputEvent(const char * event, double i=-1)
{
	time_t time0;
	struct tm tmNow;
	
	time(&time0);
	tmNow=*(localtime(&time0));
	cout<< tmNow.tm_hour<<":"<<tmNow.tm_min<<":"<<tmNow.tm_sec
		<<" -- "<<event;
	if(i>-1) cout<<i<<endl;
	else cout<<endl;
}


void handleError(const char * file, int line, const char * event=NULL, long i=-12345, long j=-12345)
{
	outputEvent("Run-Time Error");
	cout<<"!!!!! File: "<<file<<"  Line: "<<line<<endl;
	if(event!=NULL)
	{
		cout<<"!!!!! "<<event;
		if(i!=-12345) cout<<"  "<<i;
		if(j!=-12345) cout<<"  "<<j;
		cout<<endl;
	}
	cout<<"!!!!! Exit !!!!!"<<endl;
	exit(1);
}


inline void swap(long &n1, long &n2)
{
	long temp=n1;  n1=n2; n2=temp;
}


inline long getRand(long n2)
{
	if(n2<=0) return -100;
	else return rand()%n2;
}


inline long getRand(long n1, long n2)
{
	if(n1>=n2) return -100;
	else return rand()%(n2-n1)+n1;
}


inline double getRand01(void)
{
	return rand()/(RAND_MAX+Zero);
}



long binSearch(const double * list, long len, double v)
{
	if(list==NULL || len<=0) return -1;
	long top=0, bottom=len-1, middle;
	while(top<bottom) {
		middle=(top+bottom)/2;
		if(list[middle]-Eps<v && v<list[middle]+Eps) return middle;
		else if(list[middle]<v) top=middle+1;
		else bottom=middle-1;
	}

	if(list[top]-Eps<v && v<list[top]+Eps) return top;
	else if(list[top]<v) return -top-2;
	else return -top-1;
}


long binSearch(const long * list, long len, long v)
{
	if(list==NULL || len<=0) return -1;
	long top=0, bottom=len-1, middle;
	while(top<bottom) {
		middle=(top+bottom)/2;
		if(list[middle]==v) return middle;
		else if(list[middle]<v) top=middle+1;
		else bottom=middle-1;
	}

	if(list[top]==v) return top;
	else if(list[top]<v) return -top-2;
	else return -top-1;
}


long lineSearch(const long * list, long len, long v, long begin=-1)
{
	if(list==NULL || len<=0) return -1;
	if(begin<0) begin=0;
	for(long i=begin;i<len;i++) if(list[i]==v) return i;
	return -1;
}


long * enlarge(long * oldList, long oldLen, long step, long pos=-1)
{
	if(oldList==NULL || oldLen<0 || step<=0 || pos>oldLen) return NULL;
	
	long * newList=new long[oldLen+step];
	if(newList==NULL) handleError(__FILE__ , __LINE__);

	long i;
	if(pos<0) pos=oldLen;
	for(i=0;i<pos;i++) newList[i]=oldList[i];
	for(i=oldLen;i>pos;i--) newList[i]=oldList[i-1];

	delete[] oldList;
	return newList;
}


bool countSort(long *list, long len, long *b, long *c, long n)
{
	if(list==NULL || b==NULL || c==NULL) return false;
	if(len<=0 || n<=0) return false;
	long i;
	for(i=0;i<n;i++) c[i]=0;
	for(i=0;i<len;i++) c[list[i]]++;
	for(i=1;i<n;i++) c[i]+=c[i-1];
	for(i=len-1;i>=0;i--)
	{
		b[c[list[i]]-1] = list[i];
		c[list[i]]--;
	}
	for(i=0;i<len;i++) list[i]=b[i];
	return true;
}


bool countSort(long *list, long len, long n)
{
	if(list==NULL || len<=0 || n<=0) return false;
	long *b = new long[len];
	if(b==NULL) handleError(__FILE__ , __LINE__);
	long *c = new long[n];
	if(c==NULL) handleError(__FILE__ , __LINE__);

	bool r=countSort(list, len, b, c, n);
	delete[] b;
	delete[] c;
	return r;
}

#endif

