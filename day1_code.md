# P1047 [NOIP2005 普及组] 校门外的树

```c++
#include<bits/stdc++.h>
using namespace std;
#define IOS ios::sync_with_stdio(false),cin.tie(0),cout.tie(0)
int a[10005];


int main(){
//加速cin,cout的读入与输出
//如果用了，就不能用c语言读入方式，比如scanf,printf
    ios::sync_with_stdio(false),cin.tie(0),cout.tie(0);
    int l,m;
    cin>>l>>m;
    while(m--)
    {
        int u,v;
        cin>>u>>v;
        for(int i=u;i<=v;i++)
            a[i]=1;//将区域上的树拔走
    }
    int cnt=0;
    for(int i=0;i<=l;i++)
        if (a[i]==0)
            cnt++;
    cout<<cnt<<endl;


    return 0;
}
```



# P1161 开灯

```c++
#include<bits/stdc++.h>
using namespace std;
#define IOS ios::sync_with_stdio(false),cin.tie(0),cout.tie(0)
const int N=2e6+100;
bool light[N];
//0->关  1->开


int main(){
    IOS;
    int n;
    cin>>n;
    while(n--)
    {
        double a;
        int t;
        cin>>a>>t;
        for(int i=1;i<=t;i++)//将a 2*a 3*a ... t*a的灯取非操作
        {
            // int pos=floor(i*a);
            int pos=(int)(i*a);
            light[pos]=!light[pos];
        }
    }
    //遍历所有的灯，看哪个开着
    for(int i=1;i<N;i++)
        if (light[i]==true)
        {
            cout<<i<<endl;
            break;
        }

    return 0;
}
```





# P5739 【深基7.例7】计算阶乘

```c++
#include<bits/stdc++.h>
using namespace std;
#define IOS ios::sync_with_stdio(false),cin.tie(0),cout.tie(0)

int fac(int n)
{
    if (n==1) return 1;
    return n*fac(n-1);
}

int main(){
    int n;
    cin>>n;
    cout<<fac(n)<<endl;
    return 0;
}
```

