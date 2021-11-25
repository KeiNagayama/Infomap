# Infomap

## Notation

```c++
int alpha, beta;    // node id
int i, j;           // module id
vector<int> ci, cj; // module vector

double EPS=1e-16;   // epsilon used for detecting convergence

vector<vector<Flow>> adj_outflows_;  // adjacent outflows
vector<vector<Flow>> nadj_inflows_;  // adjacent inflows  normalized by total outflows
vector<vector<Flow>> nadj_outflows_; // adjacent outflows normalized by total outflows


vector<int> n2c;    // module indices for each node
```

Communityクラス
```c++
struct Community
{
	int N_;
	vector<Link> links;
	vector<vector<Flow>> adj_outflows_;
	vector<vector<Flow>> nadj_inflows_;
	vector<vector<Flow>> nadj_outflows_;
	vector<double> pagerank_;
	double tau_;
	vector<vector<int>> community_;
	vector<int> n2c_;
	CodeLength code_;
};
```

Communityコンストラクタ
(1) 初期化に使用
		- weighted links
		- tau
(2) コミュニティを指定した初期化に使用
		- weighted links
		- tau
		- com
(3) 最適コミュニティの出力に使用
		- weighted links
		- adj_outflows
		- nadj_inflows
		- nadj_outflows
		- pagerank
		- tau
		- community
		- n2c


初期化の流れ

入力: ファイル名

0. read weighted links

入力: weighted links, tau

1. get N
2. get adjacent outflows
3. get normalized adjacent outflows
4. get normalized adjacent inflows
5. get pagerank

入力: community

6. get n2c
7. get codelength


