#include <bits/stdc++.h>
#include <fstream>
using namespace std;
typedef int type;

const type INF = 2000000000;
double timer;
inline void timer_start() { timer = clock(); }
inline double timer_end() { double t = (clock() - timer) / CLOCKS_PER_SEC; return t; }

long long int mem = 0;

// A simplified queue that actually implemented by an array, use Q(size) to initialize, where size is the size of the array. 
struct Queue {
	type* node; type head, tail;
	Queue() {}
	Queue(type size) { node = (type*)malloc(size * sizeof(type)); head = tail = 0; }
	void alloc(type sz) { head = tail = 0; node = (type*)malloc(sz * sizeof(type)); }
	~Queue() { free(node); }
	bool empty() { return head == tail; }
	type pop() { return node[head++]; }
	void push(type x) { node[tail++] = x; }
	void clear() { head = tail = 0; }
};
// A simplified queue that also implemented by an array, use set(size) to initialize, use reset(size) to reinitialize. 
struct Set {
	type* nodes;
	bool* in;
	type size = -1;
	Set() {}
	Set(type sz) { size = 0; nodes = (type*)malloc(sz * sizeof(type)); in = (bool*)malloc(sz * sizeof(type)); memset(in, 0, sz * sizeof(type)); }
	void alloc(type sz) { size = 0; nodes = (type*)malloc(sz * sizeof(type)); in = (bool*)malloc(sz * sizeof(type)); memset(in, 0, sz * sizeof(type)); }
	void freememory() { free(nodes), free(in); }
	void erase(Set& s) { for (type i = 0; i < s.size; i++) in[s.nodes[i]] = false; }
	void insert(type x) { nodes[size++] = x; in[x] = true; }
	void insert(Set& s) { for (type i = 0; i < s.size; i++) insert(s.nodes[i]); }
	void clear() { for (type i = 0; i < size; i++) in[nodes[i]] = false; size = 0; }
	~Set() { free(nodes); free(in); }
};
// Edge structure, t1 and t2 are two endpoints, "to" is the endpoint the edge points to. 
struct Edge { type t1, t2, to; bool exists; };

vector<Set> RSet;
vector<Set> ASet;

set<type> totalupdateedges;

// Most computation is in the Graph structure. 
struct Graph {
	Graph() {}
	void read_edge(char* address); // Read the graph file, input the address of the file. 
	void read_R(char* address); // Read the R nodes set. 
	void read_A(char* address); // Read the A nodes set. 
	void get_densest(); // Get the nearly densest anchored subgraph using re-orientation network. 
	void generateUpdateEdges(type updatenum);
	set<type> updateedges;

	Edge* e; // Edges
	type n, m; // Number of nodes and edges. 
	type* d; // Indegree plus weight of nodes. 
	type* undeg; // Undirected degree in G.
	type** adj; // The adjacent edges of nodes. 
	type* adj_length; // The length of the adjacent list of every single node, i.e., doubled degree of nodes. 

	Set R, A; // R and A nodes set. 
	type volR;
	Set nownodes, allnodes;
	Set NearDensest; // The nearly densest anchored subgraph we obtain. 

	void init_max_d();
	type max_d, max_d_number; bool max_d_changed;
	inline void decrease_d(type u); // In the dynamic algorithms, if decrease a node's d, should invoke this function.
	inline void increase_d(type u); // In the dynamic algorithms, if increase a node's d, should invoke this function.
	void ProcessIncVertex(type u);
	void ProcessDecVertex(type u, type now_max_d);
	inline void reverse_edge(type eid);

	Queue Q; // A queue used in BFS. 

	type test_value, pivot; // The test value used in the re-orientation network, pivot = test_value - 1. 
	type alpha;
	void ReTest(type test_value);
	inline bool density_ge_test_value();
	Set ReTestSubgraph; // The dense subgraph that ReTest has found. 
	type ReTestSubgraphEdge;
	type* dist; // dist array denote "distance", used in BFS to store the distance between source and nodes. Besides, dist array also additionally replace the "visit" array, where dist[u] = -1 if u has not been visited. 
	type* cur; // used in DinicDFS, functioning as a well-known optimizing method of Dinic algorithm. 
	type* p; // p array denote "parent", used in DFS to reverse paths, storing the edge index of the BFS search tree. Besides, p array also additionally replace the "visit" array, where p[u] = -1 if u has not been visited. 
	Set s_node; // Used in re-orientation, contain all nodes directly connect to source (i.e. d[u] < pivot). 
	bool DinicBFS(); // Return whether source s can reach sink t. 
	bool DinicDFS(type); // Return whether this recursion has found an augment path. 

	void ImproveDelete(type eid);
	void ImproveInsert(type eid);

	void output(char* output_file_address);
	void output();
	void display_all_edges();

	double mybd = 0.0;
	double getlbfromSetR();

	type get_max_d(); // A simple function that returns the max indegree. 
	bool check(); // Test whether the result graph is right, if so, output True. 

	~Graph() {
		free(e); free(d); free(adj_length); free(cur); free(dist); free(p); free(undeg);

		for (int i = 0; i < n; i++) {
			free(adj[i]);
		}
		free(adj);

	}
};


Graph* G;

void Graph::read_edge(char* address) {
	// Because the input file does not tell the algorithm how many edges there are, so the "e" array is incremented step by step, the constant value below defines the size of this "step". 
	const type ARRAY_SIZE_INCREMENT = 500000;
	type Emax = ARRAY_SIZE_INCREMENT; // The size of the "e" array.

	FILE* file = fopen(address, "r");

	n = m = 0;
	e = (Edge*)malloc(Emax * sizeof(Edge));

	// Repeatedly read a line (an edge), then convert this edge to two edges (a multiple edge). 
	char line[200];
	while (fgets(line, 200, file)) {
		type i = 0;
		e[m].t1 = e[m].t2 = 0;
		while (line[i] < '0' || line[i] > '9') i++;
		while (line[i] >= '0' && line[i] <= '9') e[m].t1 = e[m].t1 * 10 + line[i] - '0', i++;
		while (line[i] < '0' || line[i] >'9') i++;
		while (line[i] >= '0' && line[i] <= '9') e[m].t2 = e[m].t2 * 10 + line[i] - '0', i++;
		n = max(n, max(e[m].t1, e[m].t2));
		e[m].to = e[m].t1; e[m].exists = true;
		m++;
		e[m].t2 = e[m - 1].t1, e[m].t1 = e[m - 1].t2;
		e[m].to = e[m].t1; e[m].exists = true;
		m++;
		if (m + 1 >= Emax) {
			Emax += ARRAY_SIZE_INCREMENT; // Increment the size of the "e" array. 
			e = (Edge*)realloc(e, Emax * sizeof(Edge));
		}
	}
	fclose(file);

	n++;

	// Allocate memory for these arrays. 
	e = (Edge*)realloc(e, m * sizeof(Edge));
	d = (type*)malloc(n * sizeof(type)); memset(d, 0, n * sizeof(type));
	dist = (type*)malloc(n * sizeof(type));
	cur = (type*)malloc(n * sizeof(type));
	p = (type*)malloc(n * sizeof(type));
	adj = (type**)malloc(n * sizeof(type*));
	adj_length = (type*)malloc(n * sizeof(type));
	undeg = (type*)malloc(n * sizeof(type));

	//mem of d, dist, cur, p

	// Allocate memory for our array-implemented sets and queue. 
	R.alloc(n); A.alloc(n); NearDensest.alloc(n); ReTestSubgraph.alloc(n); s_node.alloc(n); Q.alloc(n); nownodes.alloc(n); allnodes.alloc(n);

	//mem of ReTestSubgraph, s_node, Q
	// mem+=3*n* sizeof(type);

	// Construct the arrays "adj" and "adj_length". 
	type* now_edge_num = (type*)malloc(n * sizeof(type));
	for (type i = 0; i < n; i++) adj_length[i] = now_edge_num[i] = 0;
	for (type i = 0; i < m; i++) { adj_length[e[i].t1]++, adj_length[e[i].t2]++; }
	for (type i = 0; i < n; i++) adj[i] = (type*)malloc(adj_length[i] * sizeof(type));
	for (type i = 0; i < m; i++) {
		type t1 = e[i].t1, t2 = e[i].t2;
		adj[t1][now_edge_num[t1]++] = i;
		adj[t2][now_edge_num[t2]++] = i;
	}
	free(now_edge_num);

	for (type i = 0; i < n; i++) undeg[i] = adj_length[i] / 2;

	// printf("n = %d, m = %d, ", n, m / 2);
	return;
}

void Graph::read_R(char* address) {
	FILE* file = fopen(address, "r");

	// Repeatedly read a line (a node in R). 
	char line[200];
	while (fgets(line, 200, file)) {
		type i = 0, tem = 0;
		while (line[i] < '0' || line[i] > '9') i++;
		while (line[i] >= '0' && line[i] <= '9') tem = tem * 10 + line[i] - '0', i++;
		R.insert(tem);
	}
	fclose(file);

	printf("|R| = %d, ", R.size);
}

char buff[4194304];

vector<bool> RSetflag;
type Rsettruecnt = 0;
void read_RSet(char* address, type num) {
	type querytimes = 0;
	fstream finA(address, std::ios::in);
	finA >> querytimes;
	// cout<< querytimes<<endl;
	RSet.resize(querytimes);
	RSetflag.resize(querytimes, 1);
	for (int i = 0; i < querytimes; i++) {
		RSet[i].alloc(num);
	}
	// cout<<"emmmm"<<endl;
	for (int i = 0; i < querytimes; i++)
	{
		finA >> buff;
		// std::cout<<buff<<"\n";
		char delims[] = ",";
		char* res = strtok(buff, delims);
		if (res) {
			// std::cout<<"res="<<res<<std::endl;
			// std::cout<<"res="<<stoi(res)<<std::endl;
			RSet[i].insert(std::stoi(res));
		}
		while (res = strtok(NULL, delims)) {
			// std::cout<<"res="<<res<<std::endl;
			RSet[i].insert(std::stoi(res));
		}
		// std::sort(RSet[i].begin(), RSet[i].end());
		// std::cout<<"check R"<<i<<std::endl;
		// at least one edge check
		bool flag = false;
		for (int j = 0; j < RSet[i].size; j++)
		{
			for (int k = 0; k < RSet[i].size; k++)
			{
				if (j < k)
				{
					type u = RSet[i].nodes[j];
					type v = RSet[i].nodes[k];
					for (int j = 0; j < G->adj_length[u]; j++) {
						Edge& ne = G->e[G->adj[u][j]];
						if ((ne.t1 == u && ne.t2 == v) || (ne.t1 == v && ne.t2 == u)) {
							flag = true;
							break;
						}
					}
				}
				if (flag) {
					break;
				}
			}
			if (flag) {
				break;
			}
		}
		if (!flag) {
			RSetflag[i] = 0;
			Rsettruecnt++;
		}
		// std::cout<<"R"<<i<<"\t"<<flag<<std::endl;
		// break;
	}

	finA.close();
}


void read_ASet(char* address, type num) {

	type querytimes = 0;
	fstream finA(address, std::ios::in);
	finA >> querytimes;
	// cout<< querytimes<<endl;
	ASet.resize(querytimes);
	for (int i = 0; i < querytimes; i++) {
		ASet[i].alloc(num);
	}
	for (int i = 0; i < querytimes; i++) {
		finA >> buff;
		// std::cout<<buff<<"\n";
		char delims[] = ",";
		char* res = strtok(buff, delims);
		if (res) {
			// std::cout<<"res="<<res<<std::endl;
			ASet[i].insert(std::stoi(res));
		}
		while (res = strtok(NULL, delims)) {
			// std::cout<<"res="<<res<<std::endl;
			ASet[i].insert(std::stoi(res));
		}
		// std::sort(RSet[i].begin(), RSet[i].end());
		// std::cout<<"output A"<<i<<std::endl;
		// for(int j=0; j<ASet[i].size; j++){
		// 	std::cout<<ASet[i].nodes[j]<<"\t";
		// }
		// std::cout<<"\n";
	}
	finA.close();

	// FILE* file = fopen(address, "r");
	// // Repeatedly read a line (a node in R). 
	// char line[200];
	// fgets(line, 200, file);
	// type i = 0, querytimes = 0;
	// type querytimes=0;
	// while (line[i] < '0' || line[i] > '9') i++;
	// while (line[i] >= '0' && line[i] <= '9') querytimes = querytimes * 10 + line[i] - '0', i++;


	// while (fgets(line, 200, file)) {
	// 	type i = 0, tem = 0;
	// 	while (line[i] < '0' || line[i] > '9') i++;
	// 	while (line[i] >= '0' && line[i] <= '9') tem = tem * 10 + line[i] - '0', i++;
	// 	RSet.insert(tem);
	// }
	// fclose(file);
	// printf("|R| = %d, ", R.size);
}


void read_updateedges(char* address) {

	type querytimes = 0;
	fstream finA(address, std::ios::in);
	finA >> querytimes;
	// cout<< querytimes<<endl;
	for (int i = 0; i < querytimes; i++) {
		type eid;
		finA >> eid;
		totalupdateedges.insert(eid);
	}
	finA.close();
}



void Graph::read_A(char* address) {
	FILE* file = fopen(address, "r");

	// Repeatedly read a line (a node in A). 
	char line[200];
	while (fgets(line, 200, file)) {
		type i = 0, tem = 0;
		while (line[i] < '0' || line[i] > '9') i++;
		while (line[i] >= '0' && line[i] <= '9') tem = tem * 10 + line[i] - '0', i++;
		A.insert(tem);
	}
	fclose(file);

	printf("|A| = %d\n", A.size);
}

void Graph::get_densest() {
	for (int i = 0; i < R.size; i++) {
		int u = R.nodes[i];
		for (int j = 0; j < adj_length[u]; j++) {
			Edge& ne = e[adj[u][j]];
			int v = ne.t1 == u ? ne.t2 : ne.t1;
			if (R.in[v])
				printf("%d %d\n", u, v);
		}
	}

	// Initiate the d array and edges in R set. 
	volR = 0;
	for (type i = 0; i < R.size; i++) {
		type u = R.nodes[i];
		volR += undeg[u];
		allnodes.insert(u);
		if (A.in[u]) d[u] = undeg[u] + INF;
		else d[u] = undeg[u];
	}

	type du, dl;
	// First get the binary search range du and dl, the answer is in (dl, du]. 
	dl = du = 0;
	for (type i = 0; i < R.size; i++) du = du >= adj_length[R.nodes[i]] / 2 ? du : adj_length[R.nodes[i]] / 2;

	dl = (type)ceil(getlbfromSetR()) - 1;

	// Start binary search. 
	while (dl < du) {
		type dm = (dl + du + 1) / 2;
		ReTest(dm);
		if (density_ge_test_value()) dl = dm;
		else du = dm - 1;
	}
	type ans = dl;

	// Now we have the ans, which equals to the round-up value of the nearly densest anchored subgraph. 
	// Next invoke ReTest(ans) once and the ReTestSubgraph is the NearDensest we want. 
	ReTest(ans);
	NearDensest.insert(ReTestSubgraph);

	alpha = ans;

	// cout<<ans<<endl;

	return;
}

void Graph::generateUpdateEdges(type updatenum) {
	srand(0);
	while (updateedges.size() < updatenum) {
		type ueid = (rand() % (m));
		if (ueid % 2 == 1)
			ueid--;
		updateedges.insert(ueid);
	}
}

inline bool Graph::density_ge_test_value() {
	return ReTestSubgraphEdge > (test_value - 1) * ReTestSubgraph.size;
}

void Graph::init_max_d() {
	for (type i = 0; i < n; i++)
		if (undeg[i] >= volR)
			d[i] -= INF;
	max_d = 0;
	for (type i = 0; i < allnodes.size; i++) {
		type u = allnodes.nodes[i];
		// if (A.in[u]) continue;
		int nowd = A.in[u] ? d[u] - INF : d[u];
		if (nowd > max_d)
			max_d = nowd, max_d_number = 0;
		if (nowd == max_d)
			max_d_number++;
	}
}

inline void Graph::decrease_d(type u) {
	if (d[u] == max_d)
		max_d_number--, max_d_changed = true;
	d[u]--;
	if (max_d_number == 0) {
		max_d = 0;
		for (type i = 0; i < allnodes.size; i++) {
			u = allnodes.nodes[i];
			if (A.in[u]) continue;
			if (d[u] > max_d)
				max_d = d[u], max_d_number = 0;
			if (d[u] == max_d)
				max_d_number++;
		}
	}
}

inline void Graph::increase_d(type u) {
	if (!allnodes.in[u])
		allnodes.insert(u);
	d[u]++;
	if (A.in[u])
		return;
	if (d[u] > max_d)
		max_d = d[u], max_d_number = 1, max_d_changed = true;
	else if (d[u] == max_d)
		max_d_number++, max_d_changed = true;
}

inline void Graph::reverse_edge(type eid) {
	type from = e[eid].t1 == e[eid].to ? e[eid].t2 : e[eid].t1;
	increase_d(from);
	decrease_d(e[eid].to);
	e[eid].to = from;
}

void Graph::ProcessIncVertex(type u) {
	if (d[u] == max_d) {
		Q.clear();
		for (type i = 0; i > allnodes.size; i++) {
			type x = allnodes.nodes[i];
			p[x] = -1;
		}
		Q.push(u); p[u] = -2;
		type tar = -1;
		while (!Q.empty()) {
			type x = Q.pop();
			for (type j = 0; j < adj_length[x]; j++) {
				Edge& ne = e[adj[x][j]];
				if (!ne.exists) continue;
				if (ne.to != x) continue;
				type from = ne.t1 == x ? ne.t2 : ne.t1;
				if (d[from] <= max_d - 2) {
					p[from] = adj[x][j]; tar = from; break;
				}
				if (p[from] != -1) continue;
				p[from] = adj[x][j];
				Q.push(from);
			}
			if (tar != -1)
				break;
		}
		if (tar != -1) {
			while (tar != u) {
				type tem = e[p[tar]].to;
				reverse_edge(p[tar]);
				tar = tem;
			}
		}
	}
}

void Graph::ProcessDecVertex(type u, type now_max_d) {
	if (d[u] >= max_d - 2) {
		Q.clear();
		for (type i = 0; i < allnodes.size; i++) {
			type x = allnodes.nodes[i];
			p[x] = -1;
		}
		Q.push(u); p[u] = -2;
		type tar = -1;
		while (!Q.empty()) {
			type x = Q.pop();
			for (type j = 0; j < adj_length[x]; j++) {
				Edge& ne = e[adj[x][j]];
				if (!ne.exists) continue;
				if (ne.to == x) continue;
				if (d[ne.to] == max_d) {
					p[ne.to] = adj[x][j]; tar = ne.to; break;
				}
				if (p[ne.to] != -1) continue;
				p[ne.to] = adj[x][j];
				Q.push(ne.to);
			}
			if (tar != -1)
				break;
		}
		if (tar != -1) {
			while (tar != u) {
				reverse_edge(p[tar]);
				tar = e[p[tar]].to;
			}
		}
	}

	if (max_d < now_max_d)
		ReTest(now_max_d - 1);
}

void Graph::ReTest(type test_value_) {
	test_value = test_value_; pivot = test_value - 1;
	// cout << pivot << endl;

	if (test_value == 1) {
		ReTestSubgraph.clear();
		ReTestSubgraph.insert(R);
		goto compute;
	}

	nownodes.clear();
	for (type i = 0; i < allnodes.size; i++) {
		type u = allnodes.nodes[i];
		if (d[u] >= pivot && undeg[u] < volR)
			nownodes.insert(u);
	}

	while (DinicBFS()) {
		for (type i = 0; i < nownodes.size; i++) {
			type u = nownodes.nodes[i];
			p[u] = -1; cur[u] = 0;
		}
		for (type i = 0; i < s_node.size; i++) {
			DinicDFS(s_node.nodes[i]);
		}
	}

	// Now, the ReTestSubgraph equals to all nodes reachable to a node u, where d[u] >= pivot + 1. Next we compute it. 
	ReTestSubgraph.clear(); Q.clear();
	for (type i = 0; i < nownodes.size; i++) {
		type u = nownodes.nodes[i];
		if (d[u] >= pivot + 1 && undeg[u] < volR)
			dist[u] = 1, Q.push(u), ReTestSubgraph.insert(u);
		else
			dist[u] = 0;
	}

	while (!Q.empty()) { // Start from nodes with d[] > pivot, traversing reversely. 
		type u = Q.pop();
		for (type i = 0; i < adj_length[u]; i++) {
			Edge& ne = e[adj[u][i]];
			if (!ne.exists) continue;
			if (ne.to != u) continue;
			type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
			if (dist[from]) continue;
			dist[from] = 1; Q.push(from); ReTestSubgraph.insert(from);
		}
	}

	// Compute the number of edges in ReTestSubgraph, for computing its density.
compute:
	type Edgenumber = 0, penalty = 0;
	for (type i = 0; i < ReTestSubgraph.size; i++) {
		type u = ReTestSubgraph.nodes[i];
		if (!R.in[u]) penalty += undeg[u];
		for (type j = 0; j < adj_length[u]; j++) {
			Edge& ne = e[adj[u][j]];
			if (!ne.exists) continue;
			type v = ne.t1 == u ? ne.t2 : ne.t1;
			if (!ReTestSubgraph.in[v]) continue;
			Edgenumber++;
		}
	}
	ReTestSubgraphEdge = Edgenumber / 2 - penalty;

	return;
}

bool Graph::DinicBFS() {
	// Obtain the s_node set. 
	s_node.clear();
	for (type i = 0; i < nownodes.size; i++) if (d[nownodes.nodes[i]] > pivot && undeg[nownodes.nodes[i]] < volR) s_node.insert(nownodes.nodes[i]);

	type dist_t = INF; // The distance from source s to sink t. 

	// Initialize queue and some arrays. 
	Q.clear();
	for (type i = 0; i < nownodes.size; i++) {
		dist[nownodes.nodes[i]] = -1;
	}
	for (type i = 0; i < s_node.size; i++) dist[s_node.nodes[i]] = 1, Q.push(s_node.nodes[i]);

	// Start traversal. 
	bool break_loop = false;
	while (!Q.empty()) {
		type u = Q.pop();
		for (type i = 0; i < adj_length[u]; i++) {
			Edge& ne = e[adj[u][i]];
			if (!ne.exists) continue;
			if (ne.to != u) continue;
			type from = ne.t1 == u ? ne.t2 : ne.t1;
			if (d[from] < pivot || undeg[from] >= volR) {
				dist_t = dist[u] + 2; break_loop = true; break;
			}
			if (dist[from] != -1) continue;
			dist[from] = dist[u] + 1;
			Q.push(from);
		}
		if (break_loop) break;
	}
	return dist_t != INF;
}

bool Graph::DinicDFS(type u)
{
	if (d[u] < pivot || undeg[u] >= volR) {
		d[u]++; d[e[p[u]].to]--; e[p[u]].to = u;
		if (undeg[u] < volR && !allnodes.in[u])
			allnodes.insert(u);
		if (d[u] == pivot && undeg[u] < volR)
			nownodes.insert(u);
		return true;
	}
	for (type& i = cur[u]; i < adj_length[u]; i++) {
		Edge& ne = e[adj[u][i]];
		if (!ne.exists) continue;
		if (ne.to != u) continue;
		type from = ne.t1 == u ? ne.t2 : ne.t1;
		if ((dist[from] != dist[u] + 1) && !(d[from] < pivot || undeg[from] >= volR)) continue;
		p[from] = adj[u][i];
		if (DinicDFS(from)) {
			if (s_node.in[u]) {
				if (d[u] == pivot && undeg[u] < volR) return true; // This means arc (s, u) has been saturated, so we can return. 
				continue; // This means arc (s, u) has not been saturated, so we continue to search for augment path from u. 
			}
			d[u]++; d[e[p[u]].to]--; e[p[u]].to = u;
			return true;
		}
	}
	return false;
}

void Graph::ImproveDelete(type eid) {
	max_d_changed = false;

	type u = e[eid].t1, v = e[eid].t2;

	undeg[u]--, undeg[v]--;

	if (!R.in[u]) {
		increase_d(u);
		ProcessIncVertex(u);
	}
	if (!R.in[v]) {
		increase_d(v);
		ProcessIncVertex(v);
	}

	type x, y;
	if (e[eid].to == e[eid].t1) x = e[eid].t2, y = e[eid].t1;
	else x = e[eid].t1, y = e[eid].t2;
	e[eid].exists = false;
	decrease_d(y);
	ProcessDecVertex(y, max_d);

	eid++;
	if (e[eid].to == e[eid].t1) x = e[eid].t2, y = e[eid].t1;
	else x = e[eid].t1, y = e[eid].t2;
	e[eid].exists = false;
	decrease_d(y);
	ProcessDecVertex(y, max_d);

	if (max_d_changed) {
		NearDensest.clear(); Q.clear();
		for (type i = 0; i < allnodes.size; i++) {
			type t = allnodes.nodes[i];
			if (d[t] >= max_d)
				dist[t] = 1, Q.push(t), NearDensest.insert(t);
			else
				dist[t] = 0;
		}
		while (!Q.empty()) {
			type t = Q.pop();
			for (type j = 0; j < adj_length[t]; j++) {
				Edge& ne = e[adj[t][j]];
				if (!ne.exists) continue;
				if (ne.to != t) continue;
				type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
				if (dist[from]) continue;
				dist[from] = 1; Q.push(from); NearDensest.insert(from);
			}
		}
	}

	return;
}

void Graph::ImproveInsert(type eid) {
	max_d_changed = false;
	type u = e[eid].t1, v = e[eid].t2;

	undeg[u]++, undeg[v]++;

	if (!R.in[u]) {
		decrease_d(u);
		ProcessDecVertex(u, max_d);
	}
	if (!R.in[v]) {
		decrease_d(v);
		ProcessDecVertex(v, max_d);
	}

	type x, y;
	if (d[u] >= d[v]) x = u, y = v;
	else x = v, y = u;
	e[eid].exists = true; e[eid].to = y;
	increase_d(y);
	ProcessIncVertex(y);

	eid++;
	if (d[u] >= d[v]) x = u, y = v;
	else x = v, y = u;
	e[eid].exists = true; e[eid].to = y;
	increase_d(y);
	ProcessIncVertex(y);

	if (max_d_changed) {
		NearDensest.clear(); Q.clear();
		for (type i = 0; i < allnodes.size; i++) {
			type t = allnodes.nodes[i];
			if (d[t] >= max_d)
				dist[t] = 1, Q.push(t), NearDensest.insert(t);
			else
				dist[t] = 0;
		}
		while (!Q.empty()) {
			type t = Q.pop();
			for (type j = 0; j < adj_length[t]; j++) {
				Edge& ne = e[adj[t][j]];
				if (!ne.exists) continue;
				if (ne.to != t) continue;
				type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
				if (dist[from]) continue;
				dist[from] = 1; Q.push(from); NearDensest.insert(from);
			}
		}
	}

	return;
}

double Graph::getlbfromSetR() {
	type sn = R.size;
	type sm = 0;
	for (type i = 0; i < R.size; i++) {
		type u = R.nodes[i];
		for (type j = 0; j < adj_length[u]; j++) {
			Edge& ne = e[adj[u][j]];
			if (!ne.exists) continue;
			if (R.in[ne.t1] && R.in[ne.t2])
				sm++;
		}
	}
	sm /= 2;
	// printf("sn:%ld, sm:%ld\n", sn, sm);
	return (double)sm / sn;
}

void Graph::output(char* output_file_address) {
	// Customize what you want to output. 
	const bool OUTPUT_ANCHORED_DENSITY = true;
	const bool OUTPUT_NODES = true;
	const bool OUTPUT_D_ARRAY = false;
	const bool OUTPUT_EDGES = false;

	FILE* file = fopen(output_file_address, "w");

	if (OUTPUT_ANCHORED_DENSITY) {
		fprintf(file, "Anchored Subgraph: \n");
		printf("Anchored Subgraph: \n");
		fprintf(file, "Number of Nodes: %d. \n", NearDensest.size);
		printf("Number of Nodes: %d. \n", NearDensest.size);
		type numerator = 0;
		for (type i = 0; i < NearDensest.size; i++)
			numerator += A.in[NearDensest.nodes[i]] ? d[NearDensest.nodes[i]] - INF : d[NearDensest.nodes[i]];
		type edge_cnt = numerator;
		for (type i = 0; i < NearDensest.size; i++)
			if (!R.in[NearDensest.nodes[i]])
				edge_cnt += adj_length[NearDensest.nodes[i]] / 2;
		fprintf(file, "Number of Edges: %d * 2 = %d. \n", edge_cnt / 2, edge_cnt);
		printf("Number of Edges: %d * 2 = %d. \n", edge_cnt / 2, edge_cnt);
		fprintf(file, "Edges Minus Weight: %d. \n", numerator);
		printf("Edges Minus Weight: %d. \n", numerator);
		fprintf(file, "R-density:%lf\n", 1.0 * numerator / NearDensest.size);
		printf("R-density:%lf\n", 1.0 * numerator / NearDensest.size);
	}

	if (OUTPUT_NODES) {
		fprintf(file, "Near Densest Anchored Subgraph: ");
		printf("Near Densest Anchored Subgraph: ");
		sort(NearDensest.nodes, NearDensest.nodes + NearDensest.size);
		for (type i = 0; i < NearDensest.size; i++)	fprintf(file, "%d ", NearDensest.nodes[i]), printf("%d ", NearDensest.nodes[i]);
		fprintf(file, "\n");
		printf("\n");
	}

	if (OUTPUT_D_ARRAY) {
		fprintf(file, "D array: ");
		for (type i = 0; i < NearDensest.size; i++) fprintf(file, "d[%d] = %d\n", NearDensest.nodes[i], d[NearDensest.nodes[i]]);
	}

	if (OUTPUT_EDGES) {
		fprintf(file, "Edges: \n");
		for (type i = 0; i < NearDensest.size; i++) {
			type u = NearDensest.nodes[i];
			for (type j = 0; j < adj_length[u]; j++) {
				Edge& ne = e[adj[u][j]];
				if (ne.to != u) continue;
				type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
				fprintf(file, "%d %d\n", from, u);
			}
		}
	}
}

vector<double> history;
void Graph::output() {
	type Ecnt = 0, penalty = 0, cut = 0, vol_R_cap_S = 0, vol_S_cap_nR = 0, S_cap_R = 0;

	sort(NearDensest.nodes, NearDensest.nodes + NearDensest.size);
	printf("Densest Subgraph: ");
	for (type i = 0; i < NearDensest.size; i++) {
		type u = NearDensest.nodes[i];

		printf("%d ", u);

		if (!R.in[u]) {
			penalty += undeg[u];
			vol_S_cap_nR += undeg[u];
		}

		if (R.in[u]) {
			vol_R_cap_S += undeg[u];
			S_cap_R++;
		}

		for (type j = 0; j < adj_length[u]; j++) {
			Edge& ne = e[adj[u][j]];
			if (!ne.exists) continue;
			type v = ne.t1 == u ? ne.t2 : ne.t1;
			if (NearDensest.in[v])
				Ecnt++;
			if (!NearDensest.in[v])
				cut++;
		}
	}
	printf("\n");
	Ecnt /= 4;
	cut /= 2;

	// printf("|S|: %d\n", NearDensest.size);
	// printf("Density: %lf\n", 1.0 * Ecnt / NearDensest.size);
	// printf("Max_d: %d\n", max_d);
	printf("R-subgraph density: (2 * %d - %d) / %d = %lf\n", Ecnt, penalty, NearDensest.size, 1.0 * (2 * Ecnt - penalty) / NearDensest.size);
	// printf("Conductance: %d / %d = %lf\n", cut, min(NearDensest.size, n - NearDensest.size), 1.0 * cut / min(NearDensest.size, n - NearDensest.size));
	// assert(vol_R_cap_S != vol_S_cap_nR); printf("Local conductance: %d / (%d - %d) = %lf\n", cut, vol_R_cap_S, vol_S_cap_nR, 1.0 * cut / (vol_R_cap_S - vol_S_cap_nR));
	// printf("F1-score: 2 * %d / (%d + %d) = %lf\n", S_cap_R, NearDensest.size, R.size, 2.0 * S_cap_R / (NearDensest.size + R.size));
	history.push_back(1.0 * (2 * Ecnt - penalty) / NearDensest.size);

	long long memory = m * sizeof(Edge) + n * sizeof(type) + m * sizeof(type) + n * sizeof(type) + 7 * n * (sizeof(type) + sizeof(bool)) + n * sizeof(type) + 3 * n * sizeof(type);
	//                 e                  d                  adj                adj_length         7 * Set                                 Q                  dist+cur+p
	// printf("Memory: %lf MB = %lf GB\n", 1.0 * memory / 1024 / 1024 + 1.0 * memory / 1024 / 1024 / 1024);
}

void Graph::display_all_edges() {
	printf("All edges:------------------------\n");
	for (type u = 0; u < n; u++)
		for (type j = 0; j < adj_length[u]; j++) {
			Edge& ne = e[adj[u][j]];
			if (!ne.exists) continue;
			if (ne.to != u) continue;
			type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
			printf("%d %d\n", from, ne.to);
		}
	printf("----------------------------------\n");
}
bool Graph::check() {
	if (NearDensest.size == R.size)
		return true;
	// This function check the correctness of NearDensest, denoted as Ans for simplicity. 
	// 1. Check all edge between Ans and (V - Ans) should point to (V - Ans);
	// 2. Check the sum of d[] in Ans should equal to (2|E(Ans)| - \sum degree) (the definition of anchored density);
	// 3. Note that this function does not check whether Ans is the densest subgraph, if you want to check, please check if the result is the same as using the global algorithm.

	// 1. Check all edge between Ans and (V - Ans) should point to (V - Ans);
	for (type i = 0; i < NearDensest.size; i++) {
		type u = NearDensest.nodes[i];
		for (type j = 0; j < adj_length[u]; j++) {
			Edge& ne = e[adj[u][j]];
			if (!ne.exists) continue;
			type v = ne.t1 == u ? ne.t2 : ne.t1;
			if (NearDensest.in[v]) continue;
			if (!(ne.to == v))
				return false;
		}
	}

	// cout<<"ok1"<<endl;
	// 2. Check the sum of d[] in Ans should equal to (2|E(Ans)| - \sum degree) (the definition of anchored density);
	type sum_d = 0, numerator = 0; // numerator is (2|E(Ans)| - \sum degree)
	for (type i = 0; i < NearDensest.size; i++) sum_d += A.in[NearDensest.nodes[i]] ? d[NearDensest.nodes[i]] - INF : d[NearDensest.nodes[i]];

	// cout<<"sum_d"<<sum_d<<endl;
	for (type i = 0; i < NearDensest.size; i++) { // numerator += 2|E(Ans)|
		type u = NearDensest.nodes[i];
		for (type j = 0; j < adj_length[u]; j++) {
			Edge& ne = e[adj[u][j]];
			if (!ne.exists) continue;
			type v = ne.t1 == u ? ne.t2 : ne.t1;
			if (!NearDensest.in[v]) continue;
			numerator++;
		}
	}
	if (numerator % 4 != 0) // An edge is counted four time, so we need to divide it back then times 2. 
		return false;
	// cout<<"ok2"<<endl;
	numerator /= 2;
	// cout<<"numerator="<<numerator<<endl;
	for (type i = 0; i < NearDensest.size; i++) {
		if (!R.in[NearDensest.nodes[i]]) {
			numerator -= undeg[NearDensest.nodes[i]]; // numerator -= \sum degree
		}
	}
	// cout<< numerator<<endl;
	// cout<< sum_d<<endl;
	if (numerator != sum_d)
		return false;
	// cout<<"ok3"<<endl;
	return true;
}

int main(int argc, char** argv) {
	char adjacent_list_address[100], R_set_address[100], A_set_address[100], updateedges_address[100];
	if (argc == 5) {
		strcpy(adjacent_list_address, argv[1]); strcpy(R_set_address, argv[2]); strcpy(A_set_address, argv[3]); strcpy(updateedges_address, argv[4]);
	}
	else {
		printf("Usage:\n");
		printf("./main adjacent_list_address R_set_address A_set_address updateedges_address\n");
		return 0;
	}
	double runtime;

	printf("----------Now processing graph: %s----------\n", adjacent_list_address);


	G = new Graph();
	timer_start();
	G->read_edge(adjacent_list_address);
	runtime = timer_end();
	// printf("Reading edges runtime:\t\t%.2lf seconds\n", runtime);

	// G.read_R(R_set_address);
	// G.read_A(A_set_address);

	read_RSet(R_set_address, G->n);
	read_ASet(A_set_address, G->n);

	read_updateedges(updateedges_address);

	delete G;

	double total_insert_time = 0.0;
	double total_delete_time = 0.0;
	int realcnt = 0;

	for (type i = 0; i < RSet.size(); i++) {
		if (RSetflag[i] == 0) {
			continue;
		}
		G = new Graph();
		G->read_edge(adjacent_list_address);
		for (type j = 0; j < RSet[i].size; j++) {
			G->R.insert(RSet[i].nodes[j]);
		}
		for (type j = 0; j < ASet[i].size; j++) {
			G->A.insert(ASet[i].nodes[j]);
		}
		double total_insert_singlequery_time = 0.0;
		double total_delete_singlequery_time = 0.0;
		G->get_densest(); // Set alpha as the round-up value of current density.
		G->output();
		assert(G->check());

		// G->generateUpdateEdges(1000);
		// ???? i ??????????e[i]??e[i+1]???????????????update edges??????????eid?????
		G->init_max_d();
		cout<<"G->max_D=" <<G->max_d<<endl;
		cout<<"G->alpha=" <<G->alpha<<endl;
		if(G->max_d<3){
			delete G;
			continue;
		}
		realcnt++;
		vector<type> temedges;
		for (type eid : totalupdateedges) {
			G->updateedges.insert(eid);
		}

		for (type eid : G->updateedges)
			temedges.push_back(eid);

		for (type p = 0; p < temedges.size(); p++) {
			type eid = temedges[p];
			// eid = 20744;
			printf("----------------\nNow update edge (%d, %d) %d\n", G->e[eid].t1, G->e[eid].t2, eid);
			timer_start();
			G->ImproveDelete(eid);
			runtime = timer_end();
			total_delete_time += runtime;
			total_delete_singlequery_time += runtime;
			// G->output();
			if (G->check() == false) {
				printf("!!!!!!!!!!!!!!!the %d %d query time check error\n", i, p);
			}
			printf("Update delete time:\t%.6f seconds\n", runtime);

			timer_start();
			G->ImproveInsert(eid);
			runtime = timer_end();
			total_insert_time += runtime;
			total_insert_singlequery_time += runtime;
			// G->output();
			if (G->check() == false) {
				printf("!!!!!!!!!!!!!!!the %d %d query time check error\n", i, p);
			}
			printf("Update insert time:\t%.6f seconds\n", runtime);

			// printf("query %d cur total insert time:\t%.6f seconds\n", i, total_insert_singlequery_time);
			// printf("query %d cur total delete time:\t%.6f seconds\n", i, total_delete_singlequery_time);
			printf("query %d cur average insert time:\t%.6f seconds\n", i, total_insert_time / (p + 1));
			printf("query %d cur average delete time:\t%.6f seconds\n", i, total_delete_time / (p + 1));
		}

		// G->output();
		if (G->check() == false) {
			printf("!!!!!!!!!!!!!!!the %d query time check error\n", i);
			delete G;
			continue;
		}

		printf("query %d total insert time:\t%.6f seconds\n", i, total_insert_singlequery_time);
		printf("query %d total delete time:\t%.6f seconds\n", i, total_delete_singlequery_time);
		printf("query %d average insert time:\t%.6f seconds\n", i, total_insert_time / realcnt);
		printf("query %d average delete time:\t%.6f seconds\n", i, total_delete_time / realcnt);
		delete G;
		// break;
	}
	printf("total insert time:\t%.6f seconds\n", total_insert_time);

	printf("total delete time:\t%.6f seconds\n", total_delete_time);

	printf("average insert time:\t%.6f seconds\n", total_insert_time / (RSet.size() - Rsettruecnt));

	printf("average delete time:\t%.6f seconds\n", total_delete_time / (RSet.size() - Rsettruecnt));

	return 0;
}

