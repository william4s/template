//============================================================================
// Author      : william4s
// Description : segment tree
// Date        : 2020.9.19
//============================================================================
#include <bits/stdc++.h>

using namespace std;

const int Max = 1000;
//构造树        原数组     树数组      节点       构造的范围 
void build_tree(int arr[], int tree[], int node, int start, int last)  
{
	if(start == last)
	{
		tree[node] = arr[start];
	}
	else
	{
		int mid = (start + last) / 2;
		int left_node = 2 * node + 1;
		int right_node = 2 * node + 2;
		build_tree(arr, tree, left_node, start, mid);
		build_tree(arr, tree, right_node, mid + 1, last);
		tree[node] = tree[left_node] + tree[right_node];
	}

}
//                原数组     树数组      节点                          修改的位置   修改后数值 
void update_tree(int arr[], int tree[], int node, int start, int last, int idx, int val)
{
	if(start == last)
	{
		arr[idx] = val;
		tree[node] = val;
	}
	else
	{
		int mid = (start + last) / 2;
		int left_node = 2 * node + 1;
		int right_node = 2 * node + 2;
		if(idx >= start && idx <= mid)
		{
			update_tree(arr, tree, left_node, start, mid, idx, val);
		}
		else
		{
			update_tree(arr, tree, right_node, mid + 1, last, idx, val);
		}
		tree[node] = tree[left_node] + tree[right_node];
	}
}
// 查询区间[L,R]的和 
int query_tree(int arr[], int tree[], int node, int start, int last, int L, int R)
{
	cout << "start = " << start << endl;
	cout << "end   = " << last << endl;
	cout << endl;

	if(R < start || L > last)
	{
		return 0;
	}
	else if(L <= start && last <= R)
	{
		return tree[node];
	}
	else if(start == last)
	{
		return tree[node];
	}
	else
	{
		int mid = (start + last) / 2;
		int left_node = 2 * node + 1;
		int right_node = 2 * node + 2;
		int sum_left = query_tree(arr, tree, left_node, start, mid, L, R);
		int sum_right = query_tree(arr, tree, right_node, mid+1, last, L, R);
		return sum_left + sum_right;
	}
}

int main()
{
	std::ios::sync_with_stdio(false);
	int n = 6;
	int arr[] = {1, 3, 5, 7, 9, 11};
	int tree[Max] = {0};

	//Build Tree
	build_tree(arr, tree, 0, 0, n - 1);
	for(int i = 0; i < 15; i++)
	{
		cout << "tree[" << i << "] = " << tree[i] << endl;
	}
	cout << endl << endl;

	//Update Tree
	update_tree(arr, tree, 0, 0, n - 1, 4, 6);
	for(int i = 0; i < 15; i++)
	{
		cout << "tree[" << i << "] = " << tree[i] << endl;
	}
	cout << endl << endl;

	//Query Tree
	int s = query_tree(arr, tree, 0, 0, n-1, 2, 5);
	cout << "s = " << s << endl;
	return 0;
}
