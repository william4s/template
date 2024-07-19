# 字典树(trie)

### 一.什么是字典树

字典树很简单。就是建立一颗多叉树，字符串能够存在里面。对于有相同公共前缀的字符串而言，他们拥有相同的“祖先”。

![img](https://img-blog.csdn.net/20180823221048359?watermark/2/text/aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3dlaXhpbl8zOTc3ODU3MA==/font/5a6L5L2T/fontsize/400/fill/I0JBQkFCMA==/dissolve/70)

我们用数组$tr[MAX N][MAX CHAR]$来保存他的子节点。比如对于6号节点来说，我们接下来要找的字符是‘n’，那么$tr[6]['n'-'a']$就等于8。我们再开一个数组，$cnt[MAXN]$来保存对于一个字符串而言，我们插入了几次。比如我们用$cnt[8]$来表示在8节点处，有这样一个字符串是存在的。可以存多个相同的字符串，那么$cnt[8]$就为字符串插入几次的次数



