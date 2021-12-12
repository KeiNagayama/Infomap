# Infomap

2021.12.12

## Reference

- Multilevel Compression of Random Walks on Networks Reveals Hierarchical Organization in Large Integrated Systems
- Ranking and clustering of nodes in networks with smart teleportation
- Measuring Scholarly Impact Methods and Practice, Chapter. 1: Community detection and visualization of networks with the map equation framework

## Notations

### Flows

- pagerank = $p_{\alpha}^*$
- linkFlow = $q_{\alpha \curvearrowright \beta}$
- nodeFlow = $p_{\alpha}$

### CodeLength

- module flows
  - enterFlow to j = $q_{j \curvearrowleft}$
  - exitFlow to i = $q_{i \curvearrowright}$
  - totalFlow of i = $p_{i \circlearrowright}$
- code length
  - moduleCodeLength = $q_{\curvearrowleft}H(Q)$
  - indexCodeLength of module i = $p_{i \circlearrowright}H(P^i)$

## Structs (write only member variables)

### directed link with weight

```c++
struct Link
{
	int source;
	int target;
	double weight;
};
```

### flows w.r.t pagerank with unrecorded link teleportation

```c++
struct FlowData
{
	vector<Link> links;
	vector<Link> nlinks; // links normalized by outflows
	vector<Link> linkFlows;
	vector<vector<Link>> adj_outLinkFlows; // adjacency list of outlinkFlows
	vector<vector<Link>> adj_inLinkFlows;  // adjacency list of  inlinkFlows
	vector<double> nodeFlows;
};
```

### code length

```c++
struct CodeLength
{
	double L;
	double moduleCodeLength;
	vector<double> indexCodeLength;
	vector<double> enterFlows; // for module Code
	vector<double> exitFlows;  // for index Code
	vector<double> totalFlows; // for index Code
};
```

### community

```c++
struct Community
{
	int N;
	FlowData flowdata;
	double tau;
	vector<vector<int>> community;
	vector<int> n2c; // community indices for each node
	CodeLength code;
};
```

### delta code length

```c++
struct DeltaCodeLength
{
	int targetModule;
	double delta_codeLength;
	CodeLength new_code;
};
```

## Notes

### Node probabilities

- pagerank with unrecorded link teleportation を使用．
- ノードへの転移確率を outflow の比で与えたページランクをべき乗法で計算した後，リンクに沿った遷移のみを抽出する．
- パラメータ
  - max_iter = 100
  - tau = 0.15 (teleportation ratio)

### Optimization

1. 全ノードをシャッフルし，各 1 回ずつ最適ノードを探索する．
2. 1 を 3 回行って更新がなければ収束したものとみなす．
3. 1 の最大反復数はノード数 N とした．

### Aggregation

- 自己ループを含める
