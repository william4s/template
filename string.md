## 二维字符串哈希

Hint:两次需要使用不同的进制数

1.对$n*m$的矩阵进行二维哈希

```c++
base1=1331,base2=131;
for(int i=1;i<=n;i++)
    for(int j=1;j<=m;j++)
    {
        char x;
        cin>>x;
        h[i][j]=h[i][j-1]*base1+x;
    }
for(int i=1;i<=n;i++)
    for(int j=1;j<=m;j++)
        h[i][j]=h[i][j]+h[i-1][j]*base2;

```

2.对区域$[x1,y1,x2,y2]$求哈希

```c++
ull get(int x1,int y1,int x2,int y2)
{
    return h[x2][y2]-h[x2][y1-1]*p1[y2-y1+1]-h[x1-1][y2]*p2[x2-x1+1]
      +h[x1-1][y1-1]*p1[y2-y1+1]*p2[x2-x1+1];
}
```















# PAM

```c++
#include <bits/stdc++.h>
using namespace std;
const int maxn = 300000 + 5;

namespace pam {
int sz, tot, last;
int cnt[maxn], ch[maxn][26], len[maxn], fail[maxn];
char s[maxn];

int node(int l) {  // 建立一个新节点，长度为 l
  sz++;
  memset(ch[sz], 0, sizeof(ch[sz]));
  len[sz] = l;
  fail[sz] = cnt[sz] = 0;
  return sz;
}

void clear() {  // 初始化
  sz = -1;
  last = 0;
  s[tot = 0] = '$';
  node(0);
  node(-1);
  fail[0] = 1;
}

int getfail(int x) {  // 找后缀回文
  while (s[tot - len[x] - 1] != s[tot]) x = fail[x];
  return x;
}

void insert(char c) {  // 建树
  s[++tot] = c;
  int now = getfail(last);
  if (!ch[now][c - 'a']) {
    int x = node(len[now] + 2);
    fail[x] = ch[getfail(fail[now])][c - 'a'];
    ch[now][c - 'a'] = x;
  }
  last = ch[now][c - 'a'];
  cnt[last]++;
}

long long solve() {
  long long ans = 0;
  for (int i = sz; i >= 0; i--) {
    cnt[fail[i]] += cnt[i];
  }
  for (int i = 1; i <= sz; i++) {  // 更新答案
    ans = max(ans, 1ll * len[i] * cnt[i]);
  }
  return ans;
}
}  // namespace pam

char s[maxn];

int main() {
  pam::clear();
  scanf("%s", s + 1);
  for (int i = 1; s[i]; i++) {
    pam::insert(s[i]);
  }
  printf("%lld\n", pam::solve());
  return 0;
}
```



# AC自动机

```c++
#include<bits/stdc++.h>
using namespace std;
#define endl '\n'
typedef long long ll;
typedef pair<int,int> PII;
#define IOS ios::sync_with_stdio(false),cin.tie(0),cout.tie(0)
const int N=1e6+100;
int ch[N][26],fail[N],cnt[N],idx;
char s[N];
void insert(string str){
    int p=0;
    for(int i=0;str[i];i++){
        int u=str[i]-'a';
        if (!ch[p][u]) ch[p][u]=++idx;
        p=ch[p][u];
    }
    cnt[p]++;
}
void build(){
    queue<int> q;
    for(int i=0;i<26;i++){
        if (ch[0][i])
            q.push(ch[0][i]);
    }
    while(q.size()){
        int t=q.front();
        q.pop();
        for(int i=0;i<26;i++){
            int &p=ch[t][i];
            if (!p) p=ch[fail[t]][i];
            else{
                fail[p]=ch[fail[t]][i];
                q.push(p);
            }
        }
    }
}
int main(){
#ifndef ONLINE_JUDGE
    freopen("0.in", "r", stdin);
#endif
    IOS;
    int T;
    cin>>T;
    while(T--){
        int n;
        cin>>n;
        for(int i=0;i<=idx;i++)
        {
            memset(ch[i],0,sizeof ch[i]);
            fail[i]=0;
            cnt[i]=0;
        }
        idx=0;
        while(n--){
            string str;
            cin>>str;
            insert(str);
        }
        build();
        cin>>(s+1);
        int res=0;
        for(int t=1,i=0;s[t];t++){
            i=ch[i][s[t]-'a'];
            for(int j=i;j;j=fail[j])
            {
                res+=cnt[j];
                cnt[j]=0;
            }
        } 
        cout<<res<<endl;
    }
      

    return 0;
}
```

