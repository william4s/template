#include<bits/stdc++.h>  //��ҹ�����Һ�ɫ��˫��
using namespace std;    //����ȴ������Ѱ�ҹ���
const int INF=0x3f3f3f3f;//2020.10.9 23:59
typedef long long ll;
int  pre[1050];
bool t[1050];
int Find(int x) //��ѯ����BOSS 
{
	int r=x;
	while(r!=pre[r])
		r=pre[r];

	int i=x,j;
	while(pre[i]!=r)
	{
		j=pre[i];
		pre[i]=r;
		i=j;
	}
	return r;
}

void mix(int x,int y) //join in
{
	int fx=Find(x),fy=Find(y);
	if(fx!=fy)
	{
		pre[fy]=fx;
	}
}

int main()
{
	ios::sync_with_stdio(false);
	cin.tie(0);
	cout.tie(0);
	int N,M,a,b,i,j,ans;

	while(cin>>N>>M)
	{

		for(i=1; i<=N; i++)        //��ʼ����ÿ���˵��ϼ������Լ� 
			pre[i]=i;

		for(i=1; i<=M; i++)        //������������ͨ 
		{
			cin>>a>>b;
			mix(a,b);
		}


		memset(t,0,sizeof(t));
		for(i=1; i<=N; i++)        //��Ǹ����
		{
			t[Find(i)]=1;
		}
		for(ans=0,i=1; i<=N; i++)
			if(t[i])
				ans++;

		cout<<ans-1<<"\n";
	}


	return 0;
}




