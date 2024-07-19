#include<bits/stdc++.h>
using namespace std;
const int INF=0x3f3f3f3f;
typedef long long ll;
struct point
{
	double x=0,y=0;
};
struct v
{
	point st,end;
};
double dotProduct(v *v1,v *v2)
{
	v q,w;
	double result=0;
	q.st.x=0;
	q.st.y=0;
	q.end.x=v1->end.x-v1->st.x;
	q.end.y=v1->end.y-v1->st.y;
	
	w.st.x=0;
	w.st.y=0;
	w.end.x=v2->end.x-v2->st.x;
	w.end.y=v2->end.y-v2->st.y;
	result=q.end.x*w.end.x+q.end.y*w.end.y;
	return result;
	
}
double crossProduct(v* v1,v* v2)
{
	v q,w;
	double result;
	q.st.x=0;
	q.st.y=0;
	q.end.x=v1->end.x-(*v1).st.x;
	q.end.y=v1->end.y-v1->st.y;
	
	w.st.x=0;
	w.st.y=0;
	w.end.x=v2->end.x-v2->st.x;
	w.end.y=v2->end.y-v2->st.y;
	result=q.end.x*w.end.y-w.end.x*q.end.y;
	return result;
	
}
bool onsegment(point p1,point p2,point q)
{
	if ((q.x-p1.x)*(p2.y-p1.y)==(p2.x-p1.x)*(q.y-p1.y)&&
	min(p1.x,p2.x)<=q.x&&q.x<=max(p1.x,p2.x)&&
	min(p1.y,p2.y)<=q.y&&q.y<=max(p1.y,p2.y))
		return true;
	else
		return false;
}
struct tri
{
	point a,b,c;
};
bool intri(tri t,point p)
{
	v ab,ac,pa,pb,pc;
	ab.st=t.a;
	ab.end=t.b;
	ac.st=t.a;
	ac.end=t.c; 
	pa.st=p;
	pa.end=t.a;
	pb.st=p;
	pb.end=t.b;
	pc.st=p;
	pc.end=t.c;

	double Sabc=fabs(crossProduct(&ab,&ac));
	double Spab=fabs(crossProduct(&pa,&pb));
	double Spac=fabs(crossProduct(&pc,&pa));
	double Spbc=fabs(crossProduct(&pb,&pc));
	//cout<<Sabc<<endl<<Spab<<endl<<Spac<<endl<<Spbc<<endl;
	if (Spab+Spac+Spbc==Sabc)	return 1;
	else return 0;
	
}

int main(){
	ios::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);
	double a1,b1,a2,b2;
	point p1,p2,q;
	//cin>>p1.x>>p1.y>>p2.x>>p2.y>>q.x>>q.y;
	tri t;
	point Q;
	cin>>t.a.x>>t.a.y>>t.b.x>>t.b.y>>t.c.x>>t.c.y>>Q.x>>Q.y;
	intri(t,Q);
	
	
	
	return 0;
}



