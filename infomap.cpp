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

using Link = tuple<int, int, double>;
using Flow = pair<int, double>;

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
		links.push_back(make_tuple(i, j, w));
		if (!directed)
			links.push_back(make_tuple(j, i, w));
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
		int i = get<0>(link);
		int j = get<1>(link);
		node_count[i]++;
		node_count[j]++;
	}
	return node_count.size();
}

// weighted total_outflows
vector<double> get_total_outflows(vector<Link> &links, int N)
{
	vector<double> total_outflows(N, 0);
	for (auto flow : links)
	{
		int i = get<0>(flow);
		double w = get<2>(flow);
		total_outflows[i] += w;
	}
	return total_outflows;
}

// weighted adjacency list
vector<vector<Flow>> get_adj_outflows(int N, vector<Link> &links)
{
	vector<vector<Flow>> adj_outflows(N);
	for (auto link : links)
	{
		int i = get<0>(link);
		int j = get<1>(link);
		double w = get<2>(link);
		adj_outflows[i].push_back(make_pair(j, w));
	}
	return adj_outflows;
}

// weighted adjacency list
vector<vector<Flow>> get_adj_inflows(vector<Link> &links, int N)
{
	vector<vector<Flow>> adj_inflows(N);
	for (auto flow : links)
	{
		int i = get<0>(flow);
		int j = get<1>(flow);
		double w = get<2>(flow);
		adj_inflows[j].push_back(make_pair(i, w));
	}
	return adj_inflows;
}

// inflow list normalized by total_outflows
vector<vector<Flow>> get_nadj_inflows(vector<Link> &links, int N)
{
	vector<vector<Flow>> adj_inflows = get_adj_inflows(links, N);
	vector<vector<Flow>> nadj_inflows = adj_inflows;
	vector<double> total_outflows = get_total_outflows(links, N);
	for (int j = 0; j < nadj_inflows.size(); j++)
	{
		for (int k = 0; k < nadj_inflows[j].size(); k++)
		{
			int i = nadj_inflows[j][k].first;
			double w = nadj_inflows[j][k].second;
			nadj_inflows[j][k] = make_pair(i, w / total_outflows[i]);
		}
	}
	return nadj_inflows;
}

// outflow list normalized by total_outflows
vector<vector<Flow>> inflows2outflows(vector<vector<Flow>> &inflows)
{
	vector<vector<Flow>> outflows(inflows.size());
	for (int j = 0; j < inflows.size(); j++)
	{
		for (auto flow : inflows[j])
		{
			int i = flow.first;
			double w = flow.second;
			outflows[i].push_back(make_pair(j, w));
		}
	}
	return outflows;
}

double plogp(double p, double EPS = 1e-15)
{
	return (abs(p) < EPS) ? 0 : p * log(p);
}

double plog2p(double p, double EPS = 1e-15)
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

struct FlowSet
{
	vector<vector<Flow>> adj_outflows_;
	vector<vector<Flow>> adj_inflows_;
	vector<vector<Flow>> nadj_outflows_;
	vector<vector<Flow>> nadj_inflows_;
	vector<int> dangling_nodes_;
	vector<double> teleportation_rate_;
	// construtor
	FlowSet(
		vector<vector<Flow>> &adj_outflows,
		vector<vector<Flow>> &adj_inflows,
		vector<vector<Flow>> &nadj_outflows,
		vector<vector<Flow>> &nadj_inflows,
		vector<int> &dangling_nodes,
		vector<double> &teleportation_rate) : adj_outflows_(adj_outflows), adj_inflows_(adj_inflows), nadj_outflows_(nadj_outflows), nadj_inflows_(nadj_inflows), dangling_nodes_(dangling_nodes), teleportation_rate_(teleportation_rate){};
};

// weighted total_outflows
FlowSet get_flowset(vector<Link> &links, int N)
{
	vector<double> total_outflows(N, 0);
	double sum_outflow = 0;
	vector<vector<Flow>> adj_outflows(N);
	vector<vector<Flow>> adj_inflows(N);
	for (auto &link : links)
	{
		int i = get<0>(link);
		int j = get<1>(link);
		double w = get<2>(link);
		total_outflows[i] += w;
		sum_outflow += w;
		adj_outflows[i].push_back(make_pair(j, w));
		adj_inflows[j].push_back(make_pair(i, w));
	}
	// normalize outflows by total outflows
	vector<vector<Flow>> nadj_outflows = adj_outflows;
	for (int i = 0; i < N; i++)
	{
		double s = total_outflows[i];
		for (auto &outflow : nadj_outflows[i])
		{
			outflow.second /= s;
		}
	}
	// normalize inflows by total outflows
	vector<vector<Flow>> nadj_inflows = adj_inflows;
	for (int j = 0; j < N; j++)
	{
		for (auto &inflow : nadj_inflows[j])
		{
			int i = inflow.first;
			double s = total_outflows[i];
			inflow.second /= s;
		}
	}
	vector<int> dangling_nodes = get_dangling_nodes(total_outflows);
	vector<double> teleportation_rate = get_teleportation_rate(total_outflows, sum_outflow);
	return FlowSet(adj_outflows, adj_inflows, nadj_outflows, nadj_inflows, dangling_nodes, teleportation_rate);
}

// to calculate pagerank with unrecorded link teleportation
vector<double> get_pagerank(int N, vector<vector<Flow>> &nadj_inflows, vector<int> &dangling_nodes, vector<double> &teleportation_rate, double tau = 0.15, int T = 100000, double EPS = 1e-6)
{
	// initialization
	vector<double> pagerank(N, 1.0 / N);
	vector<double> pagerank_last(N, 1.0 / N);

	cout << "dangling nodes: " << dangling_nodes.size() << endl;

	// iteration by power method
	int t = 0;
	for (int t = 0; t < T; t++)
	{
		for (int j = 0; j < N; j++)
		{
			double x = 0;
			for (auto &flow : nadj_inflows[j])
			{
				int i = flow.first;
				double w = flow.second;
				x += w * pagerank_last[i];
			}
			double dangling_sum = get_dangling_sum(pagerank_last, dangling_nodes);
			double teleportationFlow = tau + (1.0 - tau) * dangling_sum;
			pagerank[j] = teleportationFlow * teleportation_rate[j] + (1.0 - tau) * x;
		}
		if (l1_norm(pagerank, pagerank_last) < EPS)
		{
			// damp one step (without teleportation)
			for (int j = 0; j < N; j++)
			{
				double x = 0;
				for (auto &flow : nadj_inflows[j])
				{
					int i = flow.first;
					double w = flow.second;
					x += w * pagerank_last[i];
				}
				pagerank[j] = (1.0 - tau) * x;
			}
			DivideVectorByScaler(pagerank, sum(pagerank));
			cout << "total steps to converge: " << t << endl;
			return pagerank;
		}
		pagerank_last = pagerank;
		// cout << "t = " << t << endl;
		// print_vector(pagerank, '\n', 3);
	}
	cout << "total steps to converge: " << T << " (not converged)" << endl;
	return pagerank;
}

// // to calculate pagerank
// vector<double> get_pagerank(int N, vector<vector<Flow>> &nadj_inflows, double tau = 0.15, int T = 100000, double EPS = 1e-15)
// {
// 	// initialization
// 	vector<double> pagerank(N, 1.0 / N);
// 	vector<double> pagerank_tmp(N, 1.0 / N);

// 	// iteration by power method
// 	int t = 0;
// 	while (t++ < T)
// 	{
// 		for (int j = 0; j < N; j++)
// 		{
// 			double x = 0;
// 			for (auto &flow : nadj_inflows[j])
// 			{
// 				int i = flow.first;
// 				double w = flow.second;
// 				x += w * pagerank[i];
// 			}
// 			pagerank_tmp[j] = tau / N + (1.0 - tau) * x;
// 		}
// 		if (l1_norm(pagerank, pagerank_tmp) < EPS)
// 		{
// 			// damp one step (without teleportation)
// 			for (int j = 0; j < N; j++)
// 			{
// 				double x = 0;
// 				for (auto &flow : nadj_inflows[j])
// 				{
// 					int i = flow.first;
// 					double w = flow.second;
// 					x += w * pagerank_tmp[i];
// 				}
// 				pagerank[j] = (1.0 - tau) * x;
// 			}
// 			DivideVectorByScaler(pagerank, sum(pagerank));
// 			break;
// 		}
// 		pagerank = pagerank_tmp;
// 	}
// 	return pagerank;
// }

struct CodeLength
{
	double L_;
	double q_;
	vector<double> qs_;
	vector<double> ps_;
	CodeLength(){};
	CodeLength(double L, double &q, vector<double> &qs, vector<double> &ps) : L_(L), q_(q), qs_(qs), ps_(ps){};
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
CodeLength get_code_length_0(vector<vector<Flow>> &nadj_inflows, vector<double> &pagerank, double tau, vector<vector<int>> &community, vector<int> &n2c)
{
	int N = pagerank.size();
	int NC = community.size();
	vector<vector<Flow>> nadj_outflows = inflows2outflows(nadj_inflows);
	vector<double> qs(NC, 0);
	vector<double> ps(NC, 0);

	double q_enter = 0;
	for (int i = 0; i < NC; i++)
	{
		int Ni = community[i].size();
		vector<double> pr_i = slice_vector_by_indices(pagerank, community[i]);
		double qi_tel = sum(pr_i) * (N - Ni) / N;
		double qi_adj = 0;
		for (auto &alpha : community[i])
		{
			for (auto &flow : nadj_outflows[alpha])
			{
				if (n2c[flow.first] != i)
					qi_adj += flow.second * pagerank[alpha];
			}
		}
		qs[i] = (1 - tau) * qi_adj;
		// qs[i] = tau * qi_tel + (1 - tau) * qi_adj;
		ps[i] = qs[i] + sum(pr_i);

		for (auto &alpha : community[i])
		{
			for (auto &flow : nadj_inflows[alpha])
			{
				if (n2c[flow.first] != i)
					q_enter += flow.second * pagerank[alpha];
			}
		}
	}
	// cout << "===================================" << endl;
	// cout << "test for get code length" << endl;
	// cout << "===================================" << endl;
	double q = sum(qs);
	double EPS = 1e-15;
	int precision = 4;
	// cout << "  pr = ";
	// print_vector(pagerank, '\n', precision);
	// cout << "  qs = ";
	// print_vector(qs, '\n', precision);
	// cout << "  ps = ";
	// print_vector(ps, '\n', precision);
	// cout << "H(Q) = " << entropy(qs, EPS) << endl;
	// cout << "H(P) = " << entropy(ps, EPS) << endl;
	// cout << "H(q) = " << plogp(q) << endl;
	// cout << "H(PR) = " << entropy(pagerank, EPS) << endl;
	// cout << "===================================" << endl;
	double L = 2 * entropy(qs, EPS) - entropy(ps, EPS) - plogp(q, EPS) + entropy(pagerank, EPS);

	return CodeLength(L, q, qs, ps);
}

// to calculate map equation
CodeLength get_code_length(vector<vector<Flow>> &nadj_inflows, vector<double> &pagerank, double tau, vector<vector<int>> &community, vector<int> &n2c)
{
	int N = pagerank.size();
	int NC = community.size();
	vector<vector<Flow>> nadj_outflows = inflows2outflows(nadj_inflows);
	vector<double> qs(NC, 0);
	vector<double> ps(NC, 0);

	vector<double> indexCodeLength(NC, 0);
	double moduleCodeLength = 0;

	// to calculate indexCodeLength
	for (int i = 0; i < NC; i++)
	{
		double exitFlow = 0;  // q_i_exit
		double totalFlow = 0; // p_i_circle

		double sum_flow = 0; // sum of p_alpha
		double sum_flow_log_flow = 0;

		for (auto &alpha : community[i])
		{
			for (auto &outflow : nadj_outflows[alpha])
			{
				int beta = outflow.first;
				if (n2c[beta] != i)
				{
					exitFlow += (1 - tau) * outflow.second * pagerank[alpha];
				}
			}
			double flow = pagerank[alpha];
			sum_flow += flow;
			sum_flow_log_flow += plog2p(flow);
		}
		totalFlow = exitFlow + sum_flow;
		indexCodeLength[i] = plog2p(exitFlow) + sum_flow_log_flow - plog2p(totalFlow);
	}

	// to calculate moduleCodeLength
	double sum_enterFlow = 0;
	for (int i = 0; i < NC; i++)
	{
		double enterFlow = 0; // q_i_enter

		for (auto &beta : community[i])
		{
			for (auto &inflow : nadj_inflows[beta])
			{
				int alpha = inflow.first;
				if (n2c[alpha] != i)
				{
					enterFlow += (1 - tau) * inflow.second * pagerank[alpha];
				}
			}
		}
		sum_enterFlow = enterFlow;
		moduleCodeLength += plog2p(enterFlow);
	}
	moduleCodeLength -= plog2p(sum_enterFlow);

	double codeLength = moduleCodeLength + sum(indexCodeLength);

	cout << "  indexCodeLength = ";
	print_vector(indexCodeLength, '\n');
	cout << "  moduleCodeLength = " << moduleCodeLength << endl;

	double q = sum(qs);

	return CodeLength(codeLength, q, qs, ps);
}

struct Community
{
	int N_;
	// vector<vector<Flow>> adj_outflows_;
	vector<vector<Flow>> nadj_inflows_;
	vector<vector<Flow>> nadj_outflows_;
	// vector<int> dangling_nodes_;
	vector<double> pagerank_;
	double tau_;
	vector<vector<int>> community_;
	vector<int> n2c_;
	CodeLength code_;
	Community() { this->N_ = 0; };
	// initialize with all separated partition
	Community(vector<Link> &weighted_links, double tau) : tau_(tau)
	{
		this->N_ = get_N(weighted_links);
		FlowSet flowset = get_flowset(weighted_links, this->N_);
		this->nadj_inflows_ = flowset.nadj_inflows_;
		this->nadj_outflows_ = flowset.nadj_outflows_;
		// this->dangling_nodes_ = flowset.dangling_nodes_;
		this->pagerank_ = get_pagerank(this->N_, this->nadj_inflows_, flowset.dangling_nodes_, flowset.teleportation_rate_, tau);
		this->community_ = get_init_partition(this->N_);
		this->n2c_ = get_n2c(this->N_, this->community_);
		this->code_ = get_code_length(this->nadj_inflows_, this->pagerank_, tau, this->community_, this->n2c_);
	}
	// initialize with given partition
	Community(vector<Link> &weighted_links, double tau, vector<vector<int>> &community) : tau_(tau), community_(community)
	{
		this->N_ = get_N(weighted_links);
		FlowSet flowset = get_flowset(weighted_links, this->N_);
		this->nadj_inflows_ = flowset.nadj_inflows_;
		this->nadj_outflows_ = flowset.nadj_outflows_;
		// this->dangling_nodes_ = flowset.dangling_nodes_;
		this->pagerank_ = get_pagerank(this->N_, this->nadj_inflows_, flowset.dangling_nodes_, flowset.teleportation_rate_, tau);
		this->n2c_ = get_n2c(this->N_, community);
		this->code_ = get_code_length(this->nadj_inflows_, this->pagerank_, tau, community, this->n2c_);
	}
	// do not use
	Community(vector<vector<Flow>> &nadj_inflows, vector<double> &pagerank, double tau, vector<vector<int>> &community, vector<int> &n2c) : nadj_inflows_(nadj_inflows), pagerank_(pagerank), tau_(tau), community_(community), n2c_(n2c)
	{
		this->N_ = pagerank.size();
		this->nadj_outflows_ = inflows2outflows(nadj_inflows_);
		this->code_ = get_code_length(nadj_inflows, pagerank, tau, community, n2c);
	};
	bool isNull()
	{
		return this->N_ == 0;
	}
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
	double EPS = 1e-15;
	return plogp(p2, EPS) - plogp(p1, EPS);
}

// return deltaL and updated CodeLength
pair<double, CodeLength> get_delta_L(Community &C, int alpha_0, int i_id, int j_id)
{
	vector<vector<Flow>> nadj_inflows = C.nadj_inflows_;
	vector<vector<Flow>> nadj_outflows = C.nadj_outflows_;

	vector<double> pagerank = C.pagerank_;
	double tau = C.tau_;
	vector<int> n2c = C.n2c_;
	CodeLength code = C.code_;
	vector<vector<int>> community = C.community_;
	vector<int> i = community[i_id];
	vector<int> j = community[j_id];

	int N = pagerank.size();
	int Ni = i.size();
	int Nj = j.size();

	// code length values
	double q = code.q_;
	double qi = code.qs_[i_id];
	double qj = code.qs_[j_id];
	double pi = code.ps_[i_id];
	double pj = code.ps_[j_id];

	// pageranks in module
	double p_alpha_0 = pagerank[alpha_0];

	// ========================================
	// delta_qi = delta_qi_tel + delta_qi_adj
	// ========================================
	// delta_qi_tel
	vector<double> ps_i = slice_vector_by_indices(pagerank, i);
	double delta_qi_tel = 1.0 / N * sum(ps_i) - (N - Ni + 1.0) / N * p_alpha_0;
	// delta_qi_adj
	double delta_qi_adj = 0;
	// 1st term of delta_qi_adj
	if (Ni > 1)
	{
		for (auto &inflow : nadj_inflows[alpha_0])
		{
			int alpha = inflow.first;
			if (n2c[alpha] == i_id && alpha != alpha_0)
				delta_qi_adj += inflow.second * pagerank[alpha];
		}
	}

	// 2nd term of delta_qi_adj
	for (auto &outflow : nadj_outflows[alpha_0])
	{
		if (n2c[outflow.first] != i_id)
			delta_qi_adj -= outflow.second * p_alpha_0;
	}

	// delta_qi
	double delta_qi = (1 - tau) * delta_qi_adj;
	// double delta_qi = tau * delta_qi_tel + (1 - tau) * delta_qi_adj;

	// ========================================
	// delta_qj = delta_qj_tel + delta_qj_adj
	// ========================================
	// delta_qj_tel
	vector<double> ps_j = slice_vector_by_indices(pagerank, j);
	double delta_qj_tel = -1.0 / N * sum(ps_j) + (N - Nj - 1.0) / N * p_alpha_0;

	// delta_qj_adj
	double delta_qj_adj = 0;
	// 1st term of delta_qj_adj
	for (auto &inflow : nadj_inflows[alpha_0])
	{
		int alpha = inflow.first;
		if (n2c[alpha] == j_id)
			delta_qj_adj -= inflow.second * pagerank[alpha];
	}

	// 2nd term of delta_qj_adj
	for (auto &outflow : nadj_outflows[alpha_0])
	{
		int beta = outflow.first;
		if (n2c[beta] != j_id || beta == alpha_0)
			delta_qj_adj += outflow.second * p_alpha_0;
	}
	// delta_qj
	double delta_qj = (1 - tau) * delta_qj_adj;
	// double delta_qj = tau * delta_qj_tel + (1 - tau) * delta_qj_adj;

	// ========================================
	// delta_q
	// ========================================
	double delta_q = delta_qi + delta_qj;

	// delta_pi, delta_pj
	double delta_pi = delta_qi - p_alpha_0;
	double delta_pj = delta_qj + p_alpha_0;

	// cout << "delta_qi = " << delta_qi << endl;
	// cout << "delta_qj = " << delta_qj << endl;
	// cout << "delta_pi = " << delta_pi << endl;
	// cout << "delta_pj = " << delta_pj << endl;

	// delta_L
	double qi2 = qi + delta_qi;
	double qj2 = qj + delta_qj;
	double q2 = q + delta_q;
	double pi2 = pi + delta_pi;
	double pj2 = pj + delta_pj;

	if (Ni == 1)
	{
		qi2 = 0;
		pi2 = 0;
	}
	if (Nj + 1 == N)
	{
		qj2 = 0;
	}

	double delta_L = 2 * (delta_plogp(qi, qi2) + delta_plogp(qj, qj2)) - (delta_plogp(pi, pi2) + delta_plogp(pj, pj2)) - delta_plogp(q, q2);
	// update CodeLength
	code.L_ += delta_L;
	code.q_ = q2;
	code.qs_[i_id] = qi2;
	code.qs_[j_id] = qj2;
	code.ps_[i_id] = pi2;
	code.ps_[j_id] = pj2;

	return make_pair(delta_L, code);
}

// initialization
Community initialization(vector<Link> &links, double tau)
{
	return Community(links, tau);
}

// to find a module for a node alpha to move
tuple<int, double, CodeLength> get_optimal_target_module(Community &C, int alpha)
{
	vector<int> n2c = C.n2c_;
	int i = n2c[alpha];
	int opt_j = i;
	double opt_delta_L = 0;
	CodeLength updated_code;

	for (auto &outflow : C.nadj_outflows_[alpha])
	{
		int j = n2c[outflow.first];
		pair<double, CodeLength> pr = get_delta_L(C, alpha, i, j);
		double delta_L = (i == j) ? 0 : pr.first;
		updated_code = pr.second;
		if (delta_L > opt_delta_L)
		{
			opt_j = j;
			opt_delta_L = delta_L;
		}
	}
	return make_tuple(opt_j, opt_delta_L, updated_code);
}

// to update community assignment after getting module to move alpha to
// Even if erasing alpha from module i made {}, remain null module i,
// which is removed in next process (re-indexing)
void update_community(vector<vector<int>> &community, vector<int> &n2c, int alpha, int j)
{
	// move alpha from module i to j
	int i = n2c[alpha];
	int alpha_index;
	for (int index = 0; index < community[i].size(); index++)
	{
		if (community[i][index] == alpha)
			alpha_index = index;
	}
	community[i].erase(community[i].begin() + alpha_index);
	community[j].push_back(alpha);
	// add alpha to module j
	n2c[alpha] = j;
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
	int N = C.N_;
	vector<vector<int>> community = C.community_;
	vector<int> n2c = C.n2c_;

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
				C.community_ = community;
				C.n2c_ = n2c;
				C.code_ = code;
				// // debug
				// int precision = 3;
				// cout << "\e[0;32m  updated community = \e[0m";
				// print_community(C.community_, '\n');
				// cout << "\e[0;32m  updated qs = \e[0m";
				// print_vector(C.code_.qs_, '\n', precision);
				// cout << "\e[0;32m  updated ps = \e[0m";
				// print_vector(C.code_.ps_, '\n', precision);
				// printf("\e[0;32m  L = \e[0m%1.3f\n", C.code_.L_);
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
		return Community(C.nadj_inflows_, C.pagerank_, C.tau_, new_community, n2c);
	else
		return Community();
}

vector<vector<Flow>> get_agg_outflows(vector<vector<Flow>> &flows, vector<vector<int>> &community, vector<int> &n2c)
{
	cout << "call get aggregated flows" << endl;
	int N = community.size();
	vector<vector<Flow>> agg_outflows(N);
	for (int i = 0; i < N; i++)
	{
		vector<double> agg_outflows_i_tmp(N, 0);
		for (auto alpha : community[i])
		{
			for (auto &flow : flows[alpha])
			{
				int j = n2c[flow.first];
				double w = flow.second;
				agg_outflows_i_tmp[j] += w;
			}
		}
		// normalize
		double s = 0;
		for (auto w : agg_outflows_i_tmp)
		{
			s += w;
		}
		for (int j = 0; j < N; j++)
		{
			double w = agg_outflows_i_tmp[j];
			if (w > 0)
				agg_outflows[i].push_back(make_pair(j, w / s));
		}
	}
	cout << "end get aggregated flows" << endl;
	return agg_outflows;
}

Community get_aggregated_community(Community C)
{
	cout << "call get aggregated community" << endl;
	// get C elements
	double tau = C.tau_;
	vector<vector<int>> community = C.community_;
	vector<int> n2c = C.n2c_;
	CodeLength code = C.code_;

	int N = community.size();
	vector<vector<Flow>> aggregated_nadj_inflows = get_agg_outflows(C.nadj_inflows_, community, n2c);
	vector<double> aggregated_pagerank = code.qs_;
	// assing individual community to each (aggregated) node
	vector<vector<int>> agg_community(N);
	for (int i = 0; i < N; i++)
	{
		agg_community[i] = {i};
	}
	vector<int> agg_n2c = get_n2c(N, agg_community);

	cout << "end get aggregated community" << endl;

	return Community(aggregated_nadj_inflows, aggregated_pagerank, tau, agg_community, agg_n2c);
}

vector<Community> recursive_search(vector<Link> links, double tau)
{
	vector<Community> multilevel_C;
	cout << "initialization" << endl;
	Community agg_C = initialization(links, tau);
	cout << "optimization" << endl;
	Community opt_C = get_optimal_community(agg_C);
	cout << "push C" << endl;
	multilevel_C.push_back(opt_C);
	while (!opt_C.isNull())
	{
		cout << "aggregation" << endl;
		agg_C = get_aggregated_community(opt_C);
		cout << "optimization" << endl;
		opt_C = get_optimal_community(agg_C);
		multilevel_C.push_back(opt_C);
	}
	// reverse(multilevel_C.begin(), multilevel_C.end());
	return multilevel_C;
}

vector<Community> find_community(vector<Link> links, double tau)
{
	return recursive_search(links, tau);
}

void print_multilevel_community(vector<Community> &MC)
{
	cout << "{";
	for (int i = 0; i < MC.size() - 1; i++)
	{
		print_community(MC[i].community_, '\n');
	}
	print_community(MC.back().community_);
	cout << "}" << endl;
	return;
}

void test_agg_flow()
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
	Community C = Community(links, tau, init_community);

	vector<vector<Flow>> agg_flows = get_agg_outflows(C.nadj_outflows_, C.community_, C.n2c_);

	return;
}

void test_3()
{
	string graph_name = "graphs/3.txt";
	cout << "========================================" << endl;
	cout << "test for " << graph_name << endl;

	vector<Link> links = read_links(graph_name);
	// int N = get_N(links);
	// vector<vector<Flow>> adj_flows = get_nadj_inflows(links, N);
	double tau = 0.15;

	// ==========================================================
	// test for delta L
	// ==========================================================
	cout << "========================================" << endl;
	cout << "test for delta L" << endl;

	// community 0
	cout << "========================================" << endl;
	vector<vector<int>> community0 = {{0}, {1}, {2}};
	cout << "[info]" << endl;
	cout << "  community0: ";
	print_community(community0, '\n');

	Community C0 = Community(links, tau, community0);

	int precision = 4;
	cout << "  pr = ";
	print_vector(C0.pagerank_, '\n', precision);
	cout << "  qs = ";
	print_vector(C0.code_.qs_, '\n', precision);
	cout << "  ps = ";
	print_vector(C0.code_.ps_, '\n', precision);
	double L0 = C0.code_.L_;
	double dL1_ = get_delta_L(C0, 1, 1, 2).first;
	double L1_ = get_delta_L(C0, 1, 1, 2).second.L_;
	cout << "  L  = " << L0 << endl;
	cout << "  L1_  = " << L1_ << endl;
	cout << "  dL1_  = " << dL1_ << endl;
	cout << "  test for get_delta_L >>";
	assert(L1_ == L0 + dL1_);
	cout << "\e[0;32m passed!! \e[0m" << endl;

	// community 1: moved node 1 from module 1 to 2
	cout << "========================================" << endl;
	vector<vector<int>> community1 = {{0}, {1, 2}};
	cout << "[info]" << endl;
	cout << "  community1: ";
	print_community(community1, '\n');
	Community C1 = Community(links, tau, community1);
	cout << "  pr = ";
	print_vector(C1.pagerank_, '\n', precision);
	cout << "  qs = ";
	print_vector(C1.code_.qs_, '\n', precision);
	cout << "  ps = ";
	print_vector(C1.code_.ps_, '\n', precision);
	// test delta L
	double L1 = C1.code_.L_;
	double dL2_ = get_delta_L(C1, 0, 0, 1).first;
	double L2_ = get_delta_L(C1, 0, 0, 1).second.L_;
	cout << "  L1  = " << L1 << endl;
	cout << "  L2_  = " << L2_ << endl;
	cout << "  dL2_  = " << dL2_ << endl;
	cout << "  test for get_delta_L >>";
	assert(abs(L1 - L1_) < 1e-15);
	assert(abs(L2_ - (L1 + dL2_)) < 1e-15);
	cout << "\e[0;32m passed!! \e[0m" << endl;

	// community 2: moved node 0 from module 0 to 1
	cout << "========================================" << endl;
	vector<vector<int>> community2 = {{0, 1, 2}};
	cout << "[info]" << endl;
	cout << "  community2: ";
	print_community(community2, '\n');
	Community C2 = Community(links, tau, community2);
	cout << "  pr = ";
	print_vector(C2.pagerank_, '\n');
	cout << "  qs = ";
	print_vector(C2.code_.qs_, '\n');
	cout << "  ps = ";
	print_vector(C2.code_.ps_, '\n');
	// test delta L
	double L2 = C2.code_.L_;
	cout << "  L  = " << L2 << endl;
	cout << "  test for get_delta_L >>";
	assert(L2 == L2_);
	cout << "\e[0;32m passed!! \e[0m" << endl;

	// ==========================================================
	// test for optimal target module
	// ==========================================================
	cout << "========================================" << endl;
	cout << "test for optimal target module" << endl;

	// move 1 in C0
	cout << "========================================" << endl;
	cout << "[info]" << endl;
	int alpha0 = 1;
	tuple<int, double, CodeLength> target0 = get_optimal_target_module(C0, alpha0);
	int j0 = get<0>(target0);
	printf("  move %d to module %d\n", alpha0, j0);
	cout << "  test for target module >>";
	assert(j0 != 1);
	cout << "\e[0;32m passed!! \e[0m" << endl;

	// move 0 in C1
	cout << "========================================" << endl;
	cout << "[info]" << endl;
	int alpha1 = 0;
	tuple<int, double, CodeLength> target1 = get_optimal_target_module(C0, alpha1);
	int j1 = get<0>(target1);
	printf("  move %d to module %d\n", alpha1, j1);
	cout << "  test for target module >>";
	assert(j1 == 1);
	cout << "\e[0;32m passed!! \e[0m" << endl;

	// ==========================================================
	// test for update community
	// ==========================================================
	cout << "========================================" << endl;
	cout << "test for update community" << endl;

	// udpate after "move 1 in C0"
	cout << "========================================" << endl;
	cout << "[info]" << endl;
	printf("  move %d to module %d\n", alpha0, 2);
	vector<vector<int>> c0 = C0.community_;
	vector<int> n2c0 = C0.n2c_;
	update_community(c0, n2c0, alpha0, 2);
	cout << "  test for update community >>";
	assert(c0[0].size() == 1);
	assert(c0[1].size() == 0);
	assert(c0[2].size() == 2);
	cout << "\e[0;32m passed!! \e[0m" << endl;

	// udpate after "move 0 in C1"
	cout << "========================================" << endl;
	cout << "[info]" << endl;
	printf("  move %d to module %d\n", alpha1, 1);
	vector<vector<int>> c1 = C1.community_;
	vector<int> n2c1 = C1.n2c_;
	update_community(c1, n2c1, alpha1, 1);
	cout << "  test for update community >>";
	assert(c1[1].size() == 3);
	cout << "\e[0;32m passed!! \e[0m" << endl;

	// ==========================================================
	// test for get optimal community
	// ==========================================================
	cout << "========================================" << endl;
	cout << "test for get optimal community" << endl;

	// get optimal community
	cout << "========================================" << endl;
	cout << "[info]" << endl;
	Community opt_C = get_optimal_community(C0);
	cout << "\e[0;32m  detected community = \e[0m";
	print_community(opt_C.community_, '\n');
}

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
	cout << "\e[0;32m  init pagerank = \e[0m";
	print_vector(C.pagerank_, '\n', precision);
	cout << "\e[0;32m  init qs = \e[0m";
	print_vector(C.code_.qs_, '\n', precision);
	cout << "\e[0;32m  init ps = \e[0m";
	print_vector(C.code_.ps_, '\n', precision);
	printf("\e[0;32m  L = \e[0m%1.3f\n", -C.code_.L_);
	// printf("\e[0;32m  L = \e[0m%1.3f\n", -C.code_.L_ / log(2));

	init_community = {{0}, {1}, {2, 3}};
	cout << "\e[0;32m  init community = \e[0m";
	print_community(init_community, '\n');
	C = Community(links, tau, init_community);
	printf("  L = %1.3f\n", -C.code_.L_);

	init_community = {{0, 1, 2}, {3}};
	cout << "\e[0;32m  init community = \e[0m";
	print_community(init_community, '\n');
	C = Community(links, tau, init_community);
	printf("  L = %1.3f\n", -C.code_.L_);

	init_community = {{0, 1, 2, 3}};
	cout << "\e[0;32m  init community = \e[0m";
	print_community(init_community, '\n');
	C = Community(links, tau, init_community);
	printf("  L = %1.3f\n", -C.code_.L_);

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
	// print_community(opt_C.community_, '\n');
	// cout << "\e[0;32m  L = \e[0m";
	// printf("%1.3f", -opt_C.code_.L_);

	// cout << "========================================" << endl;
	// cout << "test for find multi-level community" << endl;
	// cout << "========================================" << endl;
	// cout << "[info]" << endl;

	// vector<Community> MC = find_community(links, tau);
	// cout << "\e[0;32m  detected community = \e[0m";
	// print_multilevel_community(MC);

	return;
}

void test_27()
{
	string graph_name = "graphs/27.txt";
	cout << "========================================" << endl;
	cout << "test for " << graph_name << endl;

	vector<Link> links = read_links(graph_name);
	double tau = 0.15;

	// for (int seed = 0; seed < 10; seed++)
	// {
	// 	cout << "seed = " << seed << ", ";
	// 	vector<vector<int>> init_community(N);
	// 	for (int i = 0; i < N; i++)
	// 	{
	// 		init_community[i] = {i};
	// 	}

	// 	Community C = Community(links, tau, init_community);

	// 	// ==========================================================
	// 	// test for get optimal community
	// 	// ==========================================================
	// 	// cout << "========================================" << endl;
	// 	// cout << "test for get optimal community" << endl;

	// 	// get optimal community
	// 	// cout << "========================================" << endl;
	// 	// cout << "[info]" << endl;
	// 	Community opt_C = get_optimal_community(C, seed);
	// 	printf(", L = %1.3f\n", -opt_C.code_.L_ / log(2));
	// 	// cout << "\e[0;32m  detected community = \e[0m";
	// 	// print_community(opt_C.community_, '\n');
	// }

	// Fig.1-A
	vector<vector<int>> community = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}, {9, 10, 11}, {15, 16, 17}, {18, 19, 20}, {21, 22, 23}, {12, 13, 14}, {24, 25, 26}};
	Community C = Community(links, tau, community);
	int precision = 4;
	cout << "pagerank = ";
	print_vector(C.pagerank_, '\n', precision);
	cout << "qs = ";
	print_vector(C.code_.qs_, '\n', precision);
	cout << "ps = ";
	print_vector(C.code_.ps_, '\n', precision);
	cout << "n2c = ";
	print_vector(C.n2c_, '\n', precision);
	printf(", L = %1.3f\n", -C.code_.L_ / log(2));
}

void test_3_1()
{
	string graph_name = "graphs/3.txt";
	cout << "========================================" << endl;
	cout << "test for " << graph_name << endl;

	vector<Link> links = read_links(graph_name);
	double tau = 0;
	cout << "tau = " << tau << endl;

	vector<vector<int>> community = {{0, 1, 2}};
	Community C = Community(links, tau, community);
	int precision = 4;
	cout << "pagerank = ";
	print_vector(C.pagerank_, '\n', precision);
	cout << "qs = ";
	print_vector(C.code_.qs_, '\n', precision);
	cout << "ps = ";
	print_vector(C.code_.ps_, '\n', precision);
	cout << "n2c = ";
	print_vector(C.n2c_, '\n', precision);
	printf("L = %1.3f\n", -C.code_.L_ / log(2));
}

// n2c to community
void test_27_2()
{
	vector<int> n2c = {18, 4, 4, 4, 4, 4, 6, 6, 4, 6, 6, 6, 6, 6, 6, 6, 10, 6, 6, 6, 6, 6, 10, 4, 4, 10, 10, 4, 4, 10, 4, 10, 10, 10, 4, 4, 4, 4, 4, 2, 3, 0, 1, 0, 0, 3, 3, 1, 1, 1, 3, 0, 0, 0, 1, 3, 1, 9, 1, 1, 1, 3, 3, 3, 1, 2, 2, 0, 0, 9, 9, 9, 0, 9, 2, 2, 0, 9, 0, 0, 0, 0, 0, 1, 0, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 8, 14, 1, 0, 0, 8, 2, 0, 0, 0, 0, 0, 0, 14, 2, 0, 0, 0, 8, 0, 13, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 8, 0, 0, 2, 0, 0, 11, 11, 0, 0, 11, 0, 8, 0, 0, 13, 11, 0, 0, 0, 0, 0, 0, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 7, 0, 17, 17, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 5, 5, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 5, 12, 12, 12, 12, 5, 12};
	vector<vector<int>> community(19);
	for (int i = 0; i < n2c.size(); i++)
	{
		community[n2c[i]].push_back(i);
	}
	print_community(community);
}

void test_97()
{
	string graph_name = "graphs/USAir97_edges_reindexed.txt";
	cout << "========================================" << endl;
	cout << "test for " << graph_name << endl;

	bool directed = false;
	vector<Link> links = read_links(graph_name, directed);
	int N = get_N(links);
	// vector<vector<Flow>> adj_flows = get_nadj_inflows(links, N);
	double tau = 0.15;

	int seed = 1;
	cout << "seed = " << seed << ", ";

	vector<vector<int>> init_community(N);
	for (int i = 0; i < N; i++)
	{
		init_community[i] = {i};
	}

	Community C = Community(links, tau, init_community);
	printf(", init L = %1.3f\n", -C.code_.L_ / log(2));
	Community opt_C = get_optimal_community(C, seed);
	printf(", L = %1.3f\n", -opt_C.code_.L_ / log(2));

	for (int seed = 0; seed < 10; seed++)
	{
		cout << "seed = " << seed << ", ";
		vector<vector<int>> init_community(N);
		for (int i = 0; i < N; i++)
		{
			init_community[i] = {i};
		}

		Community C = Community(links, tau);
		Community opt_C = get_optimal_community(C, seed);
		printf(", L = %1.3f\n", -opt_C.code_.L_ / log(2));
		// cout << "\e[0;32m  detected community = \e[0m";
		// print_community(opt_C.community_, '\n');
	}

	vector<vector<int>> community = {{41, 43, 44, 51, 52, 53, 67, 68, 72, 76, 78, 79, 80, 81, 82, 84, 88, 89, 90, 91, 92, 93, 94, 95, 96, 98, 99, 100, 101, 105, 106, 109, 110, 111, 112, 113, 114, 117, 118, 119, 121, 123, 124, 125, 126, 127, 128, 129, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 146, 147, 149, 150, 152, 153, 156, 157, 159, 161, 162, 165, 166, 167, 168, 169, 170, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 229, 230, 231, 232, 233, 234, 235, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 268, 273, 274, 275, 276, 277, 278, 279, 280, 282, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 320, 321, 322, 323, 324, 325}, {42, 47, 48, 49, 54, 56, 58, 59, 60, 64, 83, 85, 86, 97, 104}, {39, 65, 66, 74, 75, 87, 108, 116, 151}, {40, 45, 46, 50, 55, 61, 62, 63}, {1, 2, 3, 4, 5, 8, 23, 24, 27, 28, 30, 34, 35, 36, 37, 38}, {312, 313, 314, 315, 316, 317, 318, 319, 326, 331}, {6, 7, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21}, {266, 267, 271, 272, 281, 283}, {102, 107, 120, 130, 145, 148, 160}, {57, 69, 70, 71, 73, 77}, {16, 22, 25, 26, 29, 31, 32, 33}, {154, 155, 158, 164, 171}, {327, 328, 329, 330, 332}, {122, 163}, {103, 115}, {195, 207}, {228, 236}, {269, 270}, {0}};

	C = Community(links, tau, community);
	cout << "n2c = ";
	print_vector(C.n2c_);
	printf(", L = %1.3f\n", -C.code_.L_ / log(2));
}

void test_lfr()
{
	string graph_name = "graphs/lfr250.txt";
	cout << "========================================" << endl;
	cout << "test for " << graph_name << endl;

	bool directed = false;
	bool use_weight = false;
	vector<Link> links = read_links(graph_name, directed, use_weight);
	int N = get_N(links);

	// vector<vector<Flow>> adj_flows = get_nadj_inflows(links, N);
	double tau = 0.15;

	int seed = 1;
	cout << "seed = " << seed << ", ";

	vector<vector<int>> init_community(N);
	for (int i = 0; i < N; i++)
	{
		init_community[i] = {i};
	}

	Community C = Community(links, tau);
	printf(", init L = %1.3f\n", C.code_.L_ / log(2));
	Community opt_C = get_optimal_community(C, seed);
	printf(", L = %1.3f\n", opt_C.code_.L_ / log(2));
	print_community(opt_C.community_, '\n');
	cout << "NC = " << opt_C.community_.size() << endl;

	vector<vector<int>> community = {{1, 7, 13, 14, 20, 23, 24, 27, 28, 29, 31, 33, 37, 38, 39, 40, 46, 54, 55, 58, 62, 69, 71, 72, 75, 76, 78, 80, 83, 86, 90, 91, 96, 102, 103, 106, 107, 109, 113, 114, 116, 117, 120, 125, 128, 132, 135, 152, 153, 160, 165, 166, 168, 169, 177, 181, 182, 188, 192, 194, 202, 204, 205, 217, 228, 235, 240, 243, 244, 247, 248, 249}, {0, 2, 3, 4, 5, 6, 9, 10, 11, 15, 17, 25, 26, 30, 32, 35, 36, 41, 42, 43, 44, 45, 47, 48, 49, 50, 51, 52, 53, 56, 57, 59, 60, 63, 64, 66, 67, 68, 73, 74, 77, 79, 81, 82, 85, 87, 88, 89, 92, 93, 94, 95, 97, 98, 99, 100, 101, 104, 105, 108, 110, 111, 112, 115, 118, 119, 121, 122, 123, 124, 126, 127, 129, 131, 133, 134, 136, 137, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 154, 155, 156, 157, 158, 159, 161, 162, 163, 164, 167, 170, 171, 172, 173, 175, 176, 178, 179, 184, 185, 187, 189, 190, 191, 193, 195, 196, 197, 198, 199, 200, 201, 203, 206, 207, 208, 210, 212, 213, 214, 215, 216, 218, 219, 220, 221, 222, 223, 224, 225, 226, 229, 230, 231, 232, 233, 234, 236, 237, 238, 239, 241, 242, 245, 246}, {8, 12, 16, 18, 19, 21, 22, 34, 61, 65, 70, 84, 130, 138, 174, 180, 183, 186, 209, 211, 227}};

	C = Community(links, tau, community);
	cout << "n2c = ";
	print_vector(C.n2c_);
	printf(", ground L = %1.3f\n", C.code_.L_ / log(2));
}

void test_entropy()
{
	double EPS = 1e-6;
	vector<double> ps = {0.5, 0.5};
	vector<double> ps2 = {0.5, 0.25, 0.25};
	cout << "  test for entropy >>";
	assert(abs(entropy(ps) - log(2)) < EPS);
	assert(abs(entropy(ps2) - 1.5 * log(2)) < EPS);
	cout << "\e[0;32m passed!! \e[0m" << endl;
}

void test_plogp()
{
	double x1 = plogp(0.5);
	cout << x1 << " " << log(2) / 2 << endl;
	double d1 = delta_plogp(0.5, 0.25);
	cout << d1 << endl;
}

void test_L()
{
	double q0 = 0.2182;
	double q1 = 0.2068;
	double p0 = 0.4750;
	double p1 = 0.9500;
	double q = q0 + q1;
	double pa0 = 0.2565;
	double pa1 = 0.4865;
	double pa2 = 0.2565;
	double L = 2 * (plogp(q0) + plogp(q1)) - (plogp(p0) + plogp(p1)) - plogp(q) + (plogp(pa0) + plogp(pa1) + plogp(pa2));
	cout << L << endl;

	// vector<double> pr = {0.2568, 0.4865, 0.2568};
	// vector<double> qs = {0.2182, 0.4135, 0.2182};
	// vector<double> ps = {0.4750, 0.9000, 0.4750};

	vector<double> pr = {0.2565, 0.4865, 0.2565};
	vector<double> qs = {0.2182, 0.2068};
	vector<double> ps = {0.4750, 0.9500};

	double L2 = 2 * entropy(qs) - entropy(ps) - plogp(sum(qs)) + entropy(pr);
	cout << L2 << endl;

	cout << "qs: " << plogp(q0) + plogp(q1) << " " << entropy(qs) << endl;
	cout << "ps: " << plogp(p0) + plogp(p1) << " " << entropy(ps) << endl;
	cout << "q: " << plogp(q) << " " << plogp(sum(qs)) << endl;
	cout << "pr: " << plogp(pa0) + plogp(pa1) + plogp(pa2) << " " << entropy(pr) << endl;
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
