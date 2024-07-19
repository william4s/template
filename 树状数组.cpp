//============================================================================
// Author      : william4s
// Description : binary tree
// Date        : 2020.9.19
//============================================================================
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <queue>
using namespace std;
int n,m;
int ans;
int he=0;
int input[500010];
struct node
{
	int left,right;
	int num;
} tree[2000010];
void build(int left,int right,int index)
{
	he++;
	tree[index].left=left;
	tree[index].right=right;
	if(left==right)
		return ;
	int mid=(right+left)/2;
	build(left,mid,index*2);
	build(mid+1,right,index*2+1);
}
int add(int index)
{
	if(tree[index].left==tree[index].right)
	{
		//cout<<index<<" "<<input[tree[index].right]<<endl;
		tree[index].num=input[tree[index].right];
		return tree[index].num;
	}
	tree[index].num=add(index*2)+add(index*2+1);
	return tree[index].num;
}
void my_plus(int index,int dis,int k)
{
	tree[index].num+=k;
	if(tree[index].left==tree[index].right)
		return ;
	if(dis<=tree[index*2].right)
		my_plus(index*2,dis,k);
	if(dis>=tree[index*2+1].left)
		my_plus(index*2+1,dis,k);
}
void search(int index,int l,int r)
{
	//cout<<index<<" ";
	if(tree[index].left>=l && tree[index].right<=r)
	{
		ans+=tree[index].num;
		return ;
	}
	if(tree[index*2].right>=l)
		search(index*2,l,r);
	if(tree[index*2+1].left<=r)
		search(index*2+1,l,r);
}
int main()
{
	cin>>n>>m;
	for(int i=1; i<=n; i++)
		scanf("%d",&input[i]);
	build(1,n,1);
	add(1);
	for(int i=1; i<=m; i++)
	{
		int a,b,c;
		scanf("%d%d%d",&a,&b,&c);
		if(a==1)
		{
			my_plus(1,b,c);
		}
		if(a==2)
		{
			ans=0;
			search(1,b,c);
			printf("%d\n",ans);
		}
	}
}
