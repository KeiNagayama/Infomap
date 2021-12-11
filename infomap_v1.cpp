// 2021/12/11 ~

#include <iostream>
#include <cstdio>  //printf
#include <iomanip> //setprecision
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <tuple>
#include <random>
#include <utility>
#include <cmath>
#include <algorithm>
#include <cassert>

using namespace std;

vector<string> split(const string &s, char delim)
{
	vector<string> elems;
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim))
	{
		elems.push_back(item);
	}
	if (s.back() == delim)
		elems.push_back(string());
	return elems;
}

struct Link
{
	int source;
	int target;
	double weight;
	Link(int source, int target, int weight) : source(source), target(target), weight(weight){};
};

void print_link(vector<Link> &links, int precision = 2)
{
	if (links.size() == 0)
	{
		cout << "{}" << endl;
		return;
	}
	for (auto &link : links)
	{
		cout << link.source << " " << link.target << " " << setprecision(precision) << link.weight << endl;
	}
	return;
}

// to read directed links from a text file
vector<Link> read_links(string fname, bool directed = true, bool use_weight = true)
{
	vector<Link> links;
	ifstream ifs(fname);
	string line;
	int i, j;
	double w = 1;
	while (getline(ifs, line))
	{
		stringstream ss(line);
		ss >> i >> j;
		if (use_weight)
			ss >> w;
		links.push_back(Link(i, j, w));
		if (!directed)
			links.push_back(Link(j, i, w));
	}
	return links;
}

double sum(vector<double> &v)
{
	double s = 0;
	for (auto &x : v)
	{
		s += x;
	}
	return s;
}

void DivideVectorByScaler(vector<double> &v, double x)
{
	for (int i = 0; i < v.size(); i++)
	{
		v[i] /= x;
	}
	return;
}

vector<double> slice_vector_by_indices(vector<double> &v, vector<int> indicies)
{
	vector<double> v1;
	for (auto &id : indicies)
	{
		v1.push_back(v[id]);
	}
	return v1;
}

double l1_norm(vector<double> &x, vector<double> &y)
{
	double d = 0;
	for (int i = 0; i < x.size(); i++)
	{
		d += abs(x[i] - y[i]);
	}
	return d;
}

template <typename T>
void print_vector(vector<T> &v, char br = '\0', int precision = 2)
{
	if (v.size() == 0)
	{
		cout << "{}" << endl;
		return;
	}
	cout << "{";
	cout << fixed << setprecision(precision);
	for (int i = 0; i < v.size() - 1; i++)
	{
		cout << v[i] << ",";
	}
	cout << v.back() << "}" << br;
	return;
}

int get_N(vector<Link> &links)
{
	map<int, int> node_count;
	for (auto link : links)
	{
		node_count[link.source]++;
		node_count[link.target]++;
	}
	return node_count.size();
}

// weighted total_outflows
vector<double> get_total_outflows(vector<Link> &links, int N)
{
	vector<double> total_outflows(N, 0);
	for (auto flow : links)
	{
		total_outflows[flow.source] += flow.weight;
	}
	return total_outflows;
}

double plogp(double p, double EPS = 1e-15)
{
	return (abs(p) < EPS) ? 0 : p * log(p) / log(2);
}

double entropy(vector<double> &ps, double EPS = 1e-15)
{
	double h = 0;
	for (auto &p : ps)
	{
		h += plogp(p, EPS);
	}
	return h;
}

vector<int> get_dangling_nodes(vector<double> &total_outflow)
{
	vector<int> dangling_nodes;
	double EPS = 1e-15;
	for (int i = 0; i < total_outflow.size(); i++)
	{
		if (total_outflow[i] < EPS)
			dangling_nodes.push_back(i);
	}
	return dangling_nodes;
}

double get_dangling_sum(vector<double> &v, vector<int> &dangling_nodes)
{
	double s = 0;
	for (auto &i : dangling_nodes)
	{
		s += v[i];
	}
	return s;
}

// to calculate teleportation rate for each node
// outflows[i] = w_i^{out}
// total_outflow = W
vector<double> get_teleportation_rate(vector<double> &total_outflows, double sum_outflow)
{
	vector<double> teleportation_rate = total_outflows;
	for (auto &x : teleportation_rate)
	{
		x /= sum_outflow;
	}
	return teleportation_rate;
}

struct FlowData
{
	vector<Link> nlinks;
	vector<Link> linkFlows;
	vector<vector<Link>> adj_outLinkFlows; // adjacency list of outlinkFlows
	vector<vector<Link>> adj_inLinkFlows;  // adjacency list of  inlinkFlows
	vector<double> nodeFlows;
	// construtor
	FlowData(){};
	FlowData(
		vector<Link> &nlinks,
		vector<Link> &linkFlows,
		vector<vector<Link>> &adj_outLinkFlows,
		vector<vector<Link>> &adj_inLinkFlows,
		vector<double> &nodeFlows)
		: nlinks(nlinks),
		  linkFlows(linkFlows),
		  adj_outLinkFlows(adj_outLinkFlows),
		  adj_inLinkFlows(adj_inLinkFlows),
		  nodeFlows(nodeFlows){};
};

// to calculate pagerank with unrecorded link teleportation
vector<double> get_pagerank(int N, vector<Link> &nlinks, vector<int> &dangling_nodes, vector<double> &teleportation_rate, double tau = 0.15, int T = 100000, double EPS = 1e-6)
{
	// initialization
	vector<double> pagerank = teleportation_rate;
	vector<double> pagerank_last = teleportation_rate;

	cout << "  #dangling nodes = " << dangling_nodes.size() << endl;
	cout << "  teleportation rate = ";
	print_vector(teleportation_rate);
	cout << endl;

	// iteration by power method
	int t = 0;
	for (int t = 0; t < T; t++)
	{
		// calc teleportation
		for (int alpha = 0; alpha < N; alpha++)
		{
			double dangling_sum = get_dangling_sum(pagerank_last, dangling_nodes);
			double teleportationFlow = tau + (1.0 - tau) * dangling_sum;
			pagerank[alpha] = teleportationFlow * teleportation_rate[alpha];
		}
		// calc flow along with link
		for (auto &link : nlinks)
		{
			int alpha = link.source;
			int beta = link.target;
			double w = link.weight;
			pagerank[beta] += (1 - tau) * w * pagerank_last[alpha];
		}
		// check convergence
		if (l1_norm(pagerank, pagerank_last) < EPS)
		{
			cout << "  total steps to converge: " << t << " (converged)" << endl;
			return pagerank;
		}
		pagerank_last = pagerank;
	}
	cout << "  total steps to converge: " << T << " (not converged)" << endl;
	return pagerank;
}

// weighted outflows
FlowData get_flowdata(vector<Link> &links, int N)
{
	// to calculate outflows and its sum
	vector<double> outflows(N, 0);
	double total_outflow = 0;
	for (auto &link : links)
	{
		double w = link.weight;
		outflows[link.source] += w;
		total_outflow += w;
	}

	// to normalize links by source node outflows
	vector<Link> nlinks = links;
	for (auto &link : nlinks)
	{
		link.weight /= outflows[link.source];
	}

	// to calculate pagerank
	vector<int> dangling_nodes = get_dangling_nodes(outflows);
	vector<double> teleportation_rate = get_teleportation_rate(outflows, total_outflow);
	vector<double> pagerank = get_pagerank(N, nlinks, dangling_nodes, teleportation_rate);
	cout << "  pagerank = ";
	print_vector(pagerank, '\n');

	// to calculate final flows
	vector<Link> linkFlows = nlinks;
	vector<vector<Link>> adj_outLinkFlows(N);
	vector<vector<Link>> adj_inLinkFlows(N);
	vector<double> nodeFlows(N, 0);
	for (auto &linkFlow : linkFlows)
	{
		double flow = linkFlow.weight * pagerank[linkFlow.source];
		linkFlow.weight = flow;
		adj_outLinkFlows[linkFlow.source].push_back(linkFlow);
		adj_inLinkFlows[linkFlow.target].push_back(linkFlow);
		nodeFlows[linkFlow.target] += flow;
	}
	return FlowData(nlinks, linkFlows, adj_outLinkFlows, adj_inLinkFlows, nodeFlows);
}

struct CodeLength
{
	double L;
	double moduleCodeLength;
	vector<double> indexCodeLength;
	vector<double> enterFlows;
	vector<double> exitFlows;
	vector<double> totalFlows;
	CodeLength(){};
	CodeLength(double L, double moduleCodeLength, vector<double> &indexCodeLength, vector<double> &enterFlows, vector<double> &exitFlows, vector<double> &totalFlows)
		: L(L), moduleCodeLength(moduleCodeLength), indexCodeLength(indexCodeLength), enterFlows(enterFlows), exitFlows(exitFlows), totalFlows(totalFlows){};
};

vector<vector<int>> get_init_partition(int N)
{
	vector<vector<int>> community(N);
	for (int alpha = 0; alpha < N; alpha++)
	{
		community[alpha] = {alpha};
	}
	return community;
}

// to get community index of each node
// if a node belongs to no community, its community index is -1
// if comunity partition is well-done, -1 does not remain
vector<int> get_n2c(int N, vector<vector<int>> &community)
{
	vector<int> n2c(N, -1);
	for (int i = 0; i < community.size(); i++)
	{
		for (auto &alpha : community[i])
		{
			n2c[alpha] = i;
		}
	}
	return n2c;
}

// to calculate map equation
CodeLength get_codelength(int N, vector<Link> &linkFlows, vector<double> &nodeFlows, double tau, vector<vector<int>> &community, vector<int> &n2c)
{
	int NC = community.size();

	// to calculate exit/enterFlows of each module
	vector<double> enterFlows(NC, 0);
	vector<double> exitFlows(NC, 0);
	for (auto &linkFlow : linkFlows)
	{
		int i = n2c[linkFlow.source];
		int j = n2c[linkFlow.target];
		double w = linkFlow.weight;
		if (i != j)
		{
			enterFlows[j] += w;
			exitFlows[i] += w;
		}
	}

	// to calculate moduleCodeLength
	double moduleCodeLength = 0;
	double sum_enterFlow = 0;
	for (int i = 0; i < NC; i++)
	{
		double enterFlow = enterFlows[i];
		sum_enterFlow += enterFlow;
		moduleCodeLength += plogp(enterFlow);
	}
	moduleCodeLength -= plogp(sum_enterFlow);

	// to calculate indexCodeLength (along with totalFlows)
	vector<double> indexCodeLength(NC, 0);
	vector<double> totalFlows(NC, 0);
	for (int i = 0; i < NC; i++)
	{
		double sum_nodeFlow = 0; // sum of p_alpha
		double sum_nodeFlow_log_nodeFlow = 0;

		for (auto &alpha : community[i])
		{
			double nodeFlow = nodeFlows[alpha];
			sum_nodeFlow += nodeFlow;
			sum_nodeFlow_log_nodeFlow += plogp(nodeFlow);
		}
		double exitFlow = exitFlows[i];
		double totalFlow = exitFlow + sum_nodeFlow;
		totalFlows[i] = totalFlow;
		indexCodeLength[i] = plogp(exitFlow) + sum_nodeFlow_log_nodeFlow - plogp(totalFlow);
	}

	// to sum up codeLength
	double codeLength = moduleCodeLength + sum(indexCodeLength);

	cout << "  indexCodeLength = ";
	print_vector(indexCodeLength);
	cout << ", " << sum(indexCodeLength) << endl;
	cout << "  moduleCodeLength = " << moduleCodeLength << endl;
	cout << "  nodeFlow = ";
	print_vector(nodeFlows, '\n');
	cout << "  exitFlow = ";
	print_vector(exitFlows, '\n');
	cout << "  enterFlow = ";
	print_vector(enterFlows, '\n');
	// cout << "  linkFlows = " << endl;
	// print_link(linkFlows);

	return CodeLength(codeLength, moduleCodeLength, indexCodeLength, enterFlows, exitFlows, totalFlows);
}

struct Community
{
	int N;
	FlowData flowdata;
	double tau;
	vector<vector<int>> community;
	vector<int> n2c;
	CodeLength code;
	// Community() { this->N = 0; };
	// initialize with all separated partition
	Community(vector<Link> &weighted_links, double tau) : tau(tau)
	{
		this->N = get_N(weighted_links);
		this->flowdata = get_flowdata(weighted_links, this->N);
		this->community = get_init_partition(this->N);
		this->n2c = get_n2c(this->N, this->community);
		this->code = get_codelength(this->N, this->flowdata.linkFlows, this->flowdata.nodeFlows, tau, this->community, this->n2c);
	};
	// initialize with given partition
	Community(vector<Link> &weighted_links, double tau, vector<vector<int>> &community) : tau(tau), community(community)
	{
		this->N = get_N(weighted_links);
		this->flowdata = get_flowdata(weighted_links, this->N);
		this->n2c = get_n2c(this->N, community);
		this->code = get_codelength(this->N, this->flowdata.linkFlows, this->flowdata.nodeFlows, tau, community, this->n2c);
	};
	// bool isNull()
	// {
	// 	return this->N == 0;
	// };
};

void print_community(vector<vector<int>> &community, char br = '\0')
{
	cout << "{";
	for (int i = 0; i < community.size() - 1; i++)
	{
		print_vector(community[i], ',');
	}
	print_vector(community.back());
	cout << "}" << br;
	return;
}

// delta(plogp)
double delta_plogp(double p1, double p2)
{
	return plogp(p2) - plogp(p1);
}

struct DeltaCodeLength
{
	int targetModule;
	double delta_codeLength;
	CodeLength new_code;
	DeltaCodeLength(int targetModule, double delta_codeLength, CodeLength &new_code)
		: targetModule(targetModule), delta_codeLength(new_codeLength), new_code(new_code){};
};

// to calculate delta codelength when moving node gamma from sourceModule to targetModule
DeltaCodeLength
get_deltaCodeLength(Community &C, int gamma, int sourceModule, int targetModule)
{
	vector<int> &n2c = C.n2c;
	double delta_source_enterFlow = 0;
	double delta_source_exitFlow = 0;
	double delta_target_enterFlow = 0;
	double delta_target_exitFlow = 0;

	assert(n2c[gamma] == sourceModule);

	// outlinkflow
	for (auto &outLinkFlow : C.flowdata.adj_outLinkFlows[gamma])
	{
		int other = n2c[outLinkFlow.target];
		// old module
		if (other == sourceModule)
			delta_source_enterFlow += outLinkFlow.weight;
		else
			delta_source_exitFlow -= outLinkFlow.weight;
		// new module
		if (other == targetModule)
			delta_target_enterFlow -= outLinkFlow.weight;
		else
			delta_target_exitFlow += outLinkFlow.weight;
	}

	// inlinkflows
	for (auto &inLinkFlow : C.flowdata.adj_inLinkFlows[gamma])
	{
		int other = n2c[inLinkFlow.source];
		// old module
		if (other != sourceModule)
			delta_source_enterFlow -= inLinkFlow.weight;
		else
			delta_source_exitFlow += inLinkFlow.weight;
		// new module
		if (other != targetModule)
			delta_target_enterFlow += inLinkFlow.weight;
		else
			delta_target_exitFlow -= inLinkFlow.weight;
	}

	double nodeFlow = C.flowdata.nodeFlows[gamma];
	double delta_source_totalFlow = delta_source_exitFlow - nodeFlow;
	double delta_target_totalFlow = delta_target_exitFlow + nodeFlow;

	// sum of enter flow
	double old_sumEnter = sum(C.code.enterFlows);
	double new_sumEnter = old_sumEnter + delta_source_enterFlow + delta_target_enterFlow;
	double delta_sumEnter_log_sumEnter = plogp(new_sumEnter) - plogp(old_sumEnter);

	// souce module flows
	double old_source_enterFlow = C.code.enterFlows[sourceModule];
	double old_source_exitFlow = C.code.exitFlows[sourceModule];
	double old_source_totalFlow = C.code.totalFlows[sourceModule];
	double new_source_enterFlow = old_source_enterFlow + delta_source_enterFlow;
	double new_source_exitFlow = old_source_exitFlow + delta_source_exitFlow;
	double new_source_totalFlow = old_source_totalFlow + delta_source_totalFlow;
	double delta_source_enter_log_enter = plogp(new_source_enterFlow) - plogp(old_source_enterFlow);
	double delta_source_exit_log_exit = plogp(new_source_exitFlow) - plogp(old_source_exitFlow);
	double delta_source_total_log_total = plogp(new_source_totalFlow) - plogp(old_source_totalFlow);

	// target module flows
	double old_target_enterFlow = C.code.enterFlows[targetModule];
	double old_target_exitFlow = C.code.exitFlows[targetModule];
	double old_target_totalFlow = C.code.totalFlows[targetModule];
	double new_target_enterFlow = old_target_enterFlow + delta_target_enterFlow;
	double new_target_exitFlow = old_target_exitFlow + delta_target_exitFlow;
	double new_target_totalFlow = old_target_totalFlow + delta_target_totalFlow;
	double delta_target_enter_log_enter = plogp(new_target_enterFlow) - plogp(old_target_enterFlow);
	double delta_target_exit_log_exit = plogp(new_target_exitFlow) - plogp(old_target_exitFlow);
	double delta_target_total_log_total = plogp(new_target_totalFlow) - plogp(old_target_totalFlow);

	double sum_delta_enter_log_enter = delta_source_enter_log_enter + delta_target_enter_log_enter;

	double delta_moduleCodeLength = sum_delta_enter_log_enter - delta_sumEnter_log_sumEnter;
	double delta_source_indexCodeLength = delta_source_exit_log_exit - delta_source_total_log_total;
	double delta_target_indexCodeLength = delta_target_exit_log_exit - delta_target_total_log_total;
	double delta_indexCodeLength = delta_source_indexCodeLength + delta_target_indexCodeLength;
	double delta_codeLength = delta_moduleCodeLength + delta_indexCodeLength;

	// cout << "--------------------------" << endl;
	// cout << "delta_source_enterFlow = " << delta_source_enterFlow << endl;
	// cout << "delta_source_exitFlow  = " << delta_source_exitFlow << endl;
	// cout << "delta_target_enterFlow = " << delta_target_enterFlow << endl;
	// cout << "delta_target_exitFlow  = " << delta_target_exitFlow << endl;

	// cout << "--------------------------" << endl;
	// cout << "old_sumEnter = " << old_sumEnter << endl;
	// cout << "new_sumEnter = " << new_sumEnter << endl;
	// cout << "delta_sumEnter_log_sumEnter = " << delta_sumEnter_log_sumEnter << endl;

	// cout << "--------------------------" << endl;
	// cout << "old_source_enterFlow = " << old_source_enterFlow << endl;
	// cout << "old_source_exitFlow  = " << old_source_exitFlow << endl;
	// cout << "old_target_enterFlow = " << old_target_enterFlow << endl;
	// cout << "old_target_exitFlow  = " << old_target_exitFlow << endl;
	// cout << "new_source_enterFlow = " << new_source_enterFlow << endl;
	// cout << "new_source_exitFlow  = " << new_source_exitFlow << endl;
	// cout << "new_target_enterFlow = " << new_target_enterFlow << endl;
	// cout << "new_target_exitFlow  = " << new_target_exitFlow << endl;

	// cout << "--------------------------" << endl;
	// cout << "delta_source_enter_log_enter = " << delta_source_enter_log_enter << endl;
	// cout << "delta_source_exit_log_exit   = " << delta_source_exit_log_exit << endl;
	// cout << "delta_source_total_log_total = " << delta_source_total_log_total << endl;
	// cout << "delta_target_enter_log_enter = " << delta_target_enter_log_enter << endl;
	// cout << "delta_target_exit_log_exit   = " << delta_target_exit_log_exit << endl;
	// cout << "delta_target_total_log_total = " << delta_target_total_log_total << endl;

	// cout << "--------------------------" << endl;
	// cout << "delta_moduleCodeLength = " << delta_moduleCodeLength << endl;
	// cout << "delta_indexCodeLength = " << delta_indexCodeLength << endl;
	// cout << "delta_codeLength = " << delta_codeLength << endl;

	// cout << "--------------------------" << endl;

	double new_codeLength = C.code.L + delta_codeLength;
	double new_moduleCodeLength = C.code.moduleCodeLength + delta_moduleCodeLength;
	vector<double> new_indexCodeLength = C.code.indexCodeLength;
	new_indexCodeLength[sourceModule] += delta_source_indexCodeLength;
	new_indexCodeLength[targetModule] += delta_target_indexCodeLength;

	vector<double> new_enterFlows = C.code.enterFlows;
	vector<double> new_exitFlows = C.code.exitFlows;
	vector<double> new_totalFlows = C.code.totalFlows;
	new_enterFlows[sourceModule] += delta_source_enterFlow;
	new_exitFlows[sourceModule] += delta_source_exitFlow;
	new_totalFlows[sourceModule] += delta_source_totalFlow;
	new_enterFlows[targetModule] += delta_target_enterFlow;
	new_exitFlows[targetModule] += delta_target_exitFlow;
	new_totalFlows[targetModule] += delta_target_totalFlow;

	CodeLength new_code = CodeLength(new_codeLength, new_moduleCodeLength, new_indexCodeLength, new_enterFlows, new_exitFlows, new_totalFlows);

	return DeltaCodeLengthData(targetModule, delta_codeLength, new_code);
}

// to find a module for a node gamma to move
DeltaCodeLengthData get_optimal_target_module(Community &C, int gamma, int sourceModule)
{
	vector<int> &n2c = C.n2c_;
	int opt_targetModule = sourceModule;
	double opt_delta_L = 0;
	DeltaCodeLength &opt_deltaCodeLength;

	for (auto &outflow : C.flowdata.adj_outLinkFlows[gamma])
	{
		int targetModule = n2c[outflow.target];
		if (targetModule == sourceModule)
			continue;
		DeltaCodeLength deltaCodeLength = get_deltaCodeLength(C, gamma, sourceModule, targetModule);
		if (delta_codeLengthData.delta_codeLength > opt_delta_L)
		{
			opt_deltaCodeLength = deltaCodeLength;
		}
	}
	return opt_deltaCodeLength;
}

// to update community assignment after getting module to move gamma to
// Even if erasing gamma from module i made {}, remain null module i,
// which is removed in next process (re-indexing)
void update_community(vector<vector<int>> &community, vector<int> &n2c, int gamma, int sourceModule, int targetModule)
{
	// update community
	int remove_index = community[sourceModule].size();
	for (int index = 0; index < community[sourceModule].size(); index++)
	{
		if (community[sourceModule][index] == gamma)
			remove_index = index;
	}
	community[sourceModule].erase(community[sourceModule].begin() + remove_index);
	community[targetModule].push_back(gamma);

	// update n2c
	n2c[gamma] = targetModule;
	return;
}

// community: return new community
// n2c: overwrite
vector<vector<int>> re_index(vector<vector<int>> &community, vector<int> &n2c)
{
	// get unique community indices
	vector<int> unq(n2c);
	sort(unq.begin(), unq.end());
	unq.erase(unique(unq.begin(), unq.end()), unq.end());
	// make a comunity index mapper (old -> new)
	map<int, int> old2new;
	for (int i = 0; i < unq.size(); i++)
	{
		old2new[unq[i]] = i;
	}
	// re-index node-community assignment
	for (int i = 0; i < n2c.size(); i++)
	{
		n2c[i] = old2new[n2c[i]];
	}
	// re-index modules
	vector<vector<int>> new_community;
	for (auto &i : unq)
	{
		new_community.push_back(community[i]);
	}
	return new_community;
}

// to find optimal community at each level
Community get_optimal_community(Community &C, int seed = 1)
{
	int N = C.N;
	vector<vector<int>> community = C.community;
	vector<int> n2c = C.n2c;

	// set random generator
	mt19937 mt(seed);
	// node indices for shuffling
	vector<int> nodes(N);
	for (int alpha = 0; alpha < N; alpha++)
	{
		nodes[alpha] = alpha;
	}

	// iteration
	int T = N; // max iteration counts
	double total_delta_L = 0;
	int convergence_counter = 0;
	for (int t = 0; t < T; t++)
	{
		// cout << "t = " << t << endl;
		double all_delta_L = 0;
		shuffle(nodes.begin(), nodes.end(), mt);
		for (auto &alpha : nodes)
		{
			// optimiation
			tuple<int, double, CodeLength> opt = get_optimal_target_module(C, alpha);
			int j = get<0>(opt);
			double delta_L = get<1>(opt);
			CodeLength code = get<2>(opt);
			// printf("optimization: alpha = %d, j = %d / delta_L = %:.3f\n", alpha, j, delta_L);
			// update community assignment and codelength
			if (delta_L > 0)
			{
				// printf("optimization: alpha = %d, j = %d / delta_L = %:.3f\n", alpha, j, delta_L);
				update_community(community, n2c, alpha, j);
				C.community = community;
				C.n2c = n2c;
				C.code = code;
				// // debug
				// int precision = 3;
				// cout << "\e[0;32m  updated community = \e[0m";
				// print_community(C.community, '\n');
				// cout << "\e[0;32m  updated qs = \e[0m";
				// print_vector(C.code.qs_, '\n', precision);
				// cout << "\e[0;32m  updated ps = \e[0m";
				// print_vector(C.code.ps_, '\n', precision);
				// printf("\e[0;32m  L = \e[0m%1.3f\n", C.code.L);
			}
			all_delta_L += delta_L;
		}
		total_delta_L += all_delta_L;
		// at convergence
		if (all_delta_L == 0)
		{
			if (convergence_counter++ > 5)
			{
				// cout << "total steps to get optimal commmunity: " << t << endl;
				break;
			}
		}
	}

	// re-index community
	vector<vector<int>> new_community = re_index(community, n2c);
	cout << "n2c = ";
	print_vector(n2c);

	if (total_delta_L > 0)
		return Community(C.nadj_inflows_, C.nodeFlow, C.tau, new_community, n2c);
	else
		return Community();
}

// vector<vector<Flow>> get_agg_outflows(vector<vector<Flow>> &flows, vector<vector<int>> &community, vector<int> &n2c)
// {
// 	cout << "call get aggregated flows" << endl;
// 	int N = community.size();
// 	vector<vector<Flow>> agg_outflows(N);
// 	for (int i = 0; i < N; i++)
// 	{
// 		vector<double> agg_outflows_i_tmp(N, 0);
// 		for (auto alpha : community[i])
// 		{
// 			for (auto &flow : flows[alpha])
// 			{
// 				int j = n2c[flow.first];
// 				double w = flow.second;
// 				agg_outflows_i_tmp[j] += w;
// 			}
// 		}
// 		// normalize
// 		double s = 0;
// 		for (auto w : agg_outflows_i_tmp)
// 		{
// 			s += w;
// 		}
// 		for (int j = 0; j < N; j++)
// 		{
// 			double w = agg_outflows_i_tmp[j];
// 			if (w > 0)
// 				agg_outflows[i].push_back(make_pair(j, w / s));
// 		}
// 	}
// 	cout << "end get aggregated flows" << endl;
// 	return agg_outflows;
// }

// Community get_aggregated_community(Community C)
// {
// 	cout << "call get aggregated community" << endl;
// 	// get C elements
// 	double tau = C.tau;
// 	vector<vector<int>> community = C.community;
// 	vector<int> n2c = C.n2c;
// 	CodeLength code = C.code;

// 	int N = community.size();
// 	vector<vector<Flow>> aggregated_nadj_inflows = get_agg_outflows(C.nadj_inflows_, community, n2c);
// 	vector<double> aggregated_pagerank = code.qs_;
// 	// assing individual community to each (aggregated) node
// 	vector<vector<int>> agg_community(N);
// 	for (int i = 0; i < N; i++)
// 	{
// 		agg_community[i] = {i};
// 	}
// 	vector<int> agg_n2c = get_n2c(N, agg_community);

// 	cout << "end get aggregated community" << endl;

// 	return Community(aggregated_nadj_inflows, aggregated_pagerank, tau, agg_community, agg_n2c);
// }

// vector<Community> recursive_search(vector<Link> links, double tau)
// {
// 	vector<Community> multilevel_C;
// 	cout << "initialization" << endl;
// 	Community agg_C = initialization(links, tau);
// 	cout << "optimization" << endl;
// 	Community opt_C = get_optimal_community(agg_C);
// 	cout << "push C" << endl;
// 	multilevel_C.push_back(opt_C);
// 	while (!opt_C.isNull())
// 	{
// 		cout << "aggregation" << endl;
// 		agg_C = get_aggregated_community(opt_C);
// 		cout << "optimization" << endl;
// 		opt_C = get_optimal_community(agg_C);
// 		multilevel_C.push_back(opt_C);
// 	}
// 	// reverse(multilevel_C.begin(), multilevel_C.end());
// 	return multilevel_C;
// }

// vector<Community> find_community(vector<Link> links, double tau)
// {
// 	return recursive_search(links, tau);
// }

// void print_multilevel_community(vector<Community> &MC)
// {
// 	cout << "{";
// 	for (int i = 0; i < MC.size() - 1; i++)
// 	{
// 		print_community(MC[i].community, '\n');
// 	}
// 	print_community(MC.back().community);
// 	cout << "}" << endl;
// 	return;
// }

void test_4()
{
	string graph_name = "graphs/4.txt";
	cout << "========================================" << endl;
	cout << "test for " << graph_name << endl;

	vector<Link> links = read_links(graph_name);
	// int N = get_N(links);
	// vector<vector<Flow>> adj_flows = get_nadj_inflows(links, N);
	double tau = 0.15;
	int seed = 1;
	int precision = 3;

	vector<vector<int>> init_community = {{0}, {1}, {2}, {3}};
	// vector<vector<int>> init_community = {{0, 1}, {2, 3}};
	cout << "\e[0;32m  init community = \e[0m";
	print_community(init_community, '\n');
	Community C = Community(links, tau, init_community);
	printf("\e[0;32m  tau = \e[0m%1.3f\n", tau);
	cout << "\e[0;32m  init nodeFlow = \e[0m";
	print_vector(C.flowdata.nodeFlows, '\n', precision);
	// cout << "\e[0;32m  init qs = \e[0m";
	// print_vector(C.code.qs_, '\n', precision);
	// cout << "\e[0;32m  init ps = \e[0m";
	// print_vector(C.code.ps_, '\n', precision);
	double L0 = C.code.L;
	double dL0 = get_deltaCodeLength(C, 3, 3, 2);
	printf("\e[0;32m  L0 = \e[0m%1.3f\n", -L0);
	printf("\e[0;32m  dL0 = \e[0m%1.3f\n", dL0);
	// printf("\e[0;32m  L = \e[0m%1.3f\n", -C.code.L / log(2));

	init_community = {{0}, {1}, {2, 3}};
	cout << "\e[0;32m  init community = \e[0m";
	print_community(init_community, '\n');
	C = Community(links, tau, init_community);
	double L1 = C.code.L;
	printf("\e[0;32m  L1 = \e[0m%1.3f\n", -L1);

	init_community = {{0, 1, 2}, {3}};
	cout << "\e[0;32m  init community = \e[0m";
	print_community(init_community, '\n');
	C = Community(links, tau, init_community);
	double L2 = C.code.L;
	double dL2 = get_deltaCodeLength(C, 3, 1, 0);
	printf("\e[0;32m  L2 = \e[0m%1.3f\n", -L2);
	printf("\e[0;32m  dL2 = \e[0m%1.3f\n", dL2);

	init_community = {{0, 1, 2, 3}};
	cout << "\e[0;32m  init community = \e[0m";
	print_community(init_community, '\n');
	C = Community(links, tau, init_community);
	double L3 = C.code.L;
	printf("\e[0;32m  L3 = \e[0m%1.3f\n", -L3);

	// ==========================================================
	// test for get optimal community
	// ==========================================================
	// cout << "========================================" << endl;
	// cout << "test for get optimal community" << endl;

	// // get optimal community
	// cout << "========================================" << endl;
	// cout << "[info]" << endl;
	// Community opt_C = get_optimal_community(C, seed);
	// cout << "\e[0;32m  seed = \e[0m" << seed << endl;
	// cout << "\e[0;32m  detected community = \e[0m";
	// print_community(opt_C.community, '\n');
	// cout << "\e[0;32m  L = \e[0m";
	// printf("%1.3f", -opt_C.code.L);

	// cout << "========================================" << endl;
	// cout << "test for find multi-level community" << endl;
	// cout << "========================================" << endl;
	// cout << "[info]" << endl;

	// vector<Community> MC = find_community(links, tau);
	// cout << "\e[0;32m  detected community = \e[0m";
	// print_multilevel_community(MC);

	return;
}

int main(int argc, char *argv[])
{
	// test_L();
	// test_3(); // all test passed!! (for one-layer search)
	test_4();
	// test_27();
	// test_3_1();
	// test_97();
	// test_lfr();
	// test_entropy();
	// test_plogp();
	return 0;
}
