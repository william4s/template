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
double mul(point p1,point p2,point p0)
{
	return (p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
}
ll across(v v1,v v2)
{
	if (max(v1.st.x,v1.end.x)>=min(v2.st.x,v2.end.x)&&
	        max(v2.st.x,v2.end.x)>=min(v1.st.x,v1.end.x)&&
	        max(v1.st.y,v1.end.y)>=min(v2.st.y,v2.end.y)&&
	        mul(v2.st,v1.end,v1.st)*mul(v1.end,v2.end,v1.st)>=0&&  //修改这里与下面的>=改变端点相交、线段重合 
	        mul(v1.st,v2.end,v2.st)*mul(v2.end,v1.end,v2.st)>=0)
		return 1;
	return 0;
}
int main()
{
	ios::sync_with_stdio(false);
	cin.tie(0);
	cout.tie(0);
	v v1,v2;
	cin>>v1.st.x>>v1.st.y>>v1.end.x>>v1.end.y;
	cin>>v2.st.x>>v2.st.y>>v2.end.x>>v2.end.y;
	cout<<across(v1,v2);
	return 0;
}



