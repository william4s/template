#include<iostream>
#include<cmath>
#include<cstdio>
#include<algorithm>
using namespace std;
const int INF=0x3f3f3f3f;
typedef long long ll;

#define eps 1.0e-5


struct point
{
	double x;
	double y;
} po[111];  //几个点，记得改 

ll n,m;

bool online(point p1,point p,point p2)
{
	if( p.x<=max(p1.x,p2.x) && p.x>=min(p1.x,p2.x) && p.y<=max(p1.y,p2.y) && p.y>=min(p1.y,p2.y) )
	{
		if ( fabs(((p.x-p1.x)*(p2.y-p1.y) - (p.y-p1.y)*(p2.x-p1.x)))<=eps )
			return true;
	}
	return false;
}
bool inside(point p)
{
	ll cnt = 0;
	double xinter;
	point p1,p2;
	p1 = po[0];
	for(int i=1; i<=n; i++)
	{
		p2 = po[i%n];
		if( online(p1,p,p2) )	return true;
		if( p.x<=max(p1.x,p2.x)  && p.y<=max(p1.y,p2.y) && p.y>min(p1.y,p2.y) )
		{
			if(p1.y!=p2.y)
			{
				xinter = (p.y-p1.y)*(p1.x-p2.x)/(p1.y-p2.y) + p1.x;
				if(p1.x==p2.x || p.x<=xinter)	cnt++;
			}
		}
		p1 = p2;
	}
	if(cnt%2==0)	return false;
	else			return true;
}

int main()
{
	ios::sync_with_stdio(false);
	cin.tie(0);
	cout.tie(0); 
	ll kase = 0;
	point p;
	while(cin>>n)
	{
		if (n==0) break;
		cin>>m;
		if(kase)	cout<<"\n";
	
		for(int i=0; i<n; i++)	cin>>po[i].x>>po[i].y;
		cout<<"Problem "<<++kase<<":\n";
		for(int i=0; i<m; i++)
		{
			cin>>p.x>>p.y;
			if( inside(p) )		cout<<"Within\n";
			else				cout<<"Outside\n";
		}
	}
	return 0;

}
