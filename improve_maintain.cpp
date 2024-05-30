// Improve dynamic algorithms ImproveIns and ImproveDel to maintain AADS
// Usage: ./main <graph_dataset_address> <R_set_address> <A_set_address>
#include <bits/stdc++.h>
using namespace std;
const int INF = 1000000000;

struct Timer {
	double start_time, end_time;
	void start() { start_time = clock(); }
	void end() { end_time = clock(); }
	double time() { return (end_time - start_time) / CLOCKS_PER_SEC; }
};
template <class T>
struct Set {
	T* nodes; bool* in; int size = -1;
	Set() {}
	Set(int sz) { size = 0; nodes = (T*)malloc(sz * sizeof(T)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); }
	void alloc(int sz) { size = 0; nodes = (T*)malloc(sz * sizeof(int)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); }
	void insert(T x) { nodes[size++] = x; in[x] = true; }
	void clear() { for (int i = 0; i < size; i++) in[nodes[i]] = false; size = 0; }
	~Set() { free(nodes), free(in); }
	void operator=(Set& S) { clear(); for (int i = 0; i < S.size; i++) insert(S.nodes[i]); }
};
template <class T>
struct Map {
	int* nodes; bool* in; int size = -1; T* value;
	Map() {}
	Map(int sz) { size = 0; nodes = (int*)malloc(sz * sizeof(int)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); value = (T*)malloc(sz * sizeof(T)); memset(value, 0, sz * sizeof(T)); }
	void alloc(int sz) { size = 0; nodes = (int*)malloc(sz * sizeof(int)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); value = (T*)malloc(sz * sizeof(T)); memset(value, 0, sz * sizeof(T)); }
	void freememory() { free(nodes), free(in), free(value); }
	void clear() { for (int i = 0; i < size; i++) in[nodes[i]] = false, value[nodes[i]] = 0; size = 0; }
	T& operator[](int x) { if (!in[x]) nodes[size++] = x, in[x] = true; return value[x]; }
	~Map() { free(nodes), free(in), free(value); }
};
template <class T>
struct Queue {
	T* nodes; int head, tail;
	Queue() {}
	Queue(int size) { nodes = (T*)malloc(size * sizeof(T)); head = tail = 0; }
	void alloc(int sz) { head = tail = 0; nodes = (int*)malloc(sz * sizeof(int)); }
	~Queue() { free(nodes); }
	bool empty() { return head == tail; }
	int pop() { return nodes[head++]; }
	void push(T x) { nodes[tail++] = x; }
	void clear() { head = tail = 0; }
};

inline int read_number(FILE* in) {
	int x = 0; char ch = 0; while (ch < '0' || ch > '9') ch = fgetc(in); while (ch >= '0' && ch <= '9') { x = x * 10 + (ch - '0'); ch = fgetc(in); } return x;
}
inline void check(bool flag, const char* message) {
	if (!flag) {
		printf("!!!!! CHECK ERROR !!!!!\n");
		printf("Error message: %s\n", message);
		assert(0);
	}
}

struct Edge { int t1, t2, to; bool exists; };
struct Graph {
	int n, m;
	Edge* e;
	Set<int> R, A;
	int* undeg, * bounty, * adj_length;
	int** adj;
	Set<int> V_A;
	void read_graph_from_dataset(FILE* dataset_file);
	bool read_R_and_A_set(FILE* R_file, FILE* A_file);

	void initialize();

	double allquery_ave_delete_time=0.0;
	double allquery_ave_insert_time=0.0;


	Set<int> ReTestSubgraph; int ReTestSubgraphBounty;
	void get_AADS();

	int test_value, pivot, alpha;
	Queue<int> Q;
	Map<int> parent, cur, dist; Set<int> vis;
	void ReTest(int test_value_);
	bool DinicBFS();
	bool DinicDFS(int u);

	void get_update_edge(); set<int> edges_to_be_updated;
	void maintain();
	void edge_delete(int eid);
	void edge_insert(int eid);
	void IncBounty(int x);
	void DecBounty(int x, int eid);
	int count_bounty_in_set(Set<int>& S);

	void read_update_edge(FILE* dataset_file);

	Set<int> S_alpha, S_alpha1; int S_alpha_bounty, S_alpha1_bounty;
	void check_correctness();
	void output_AADS();
	void display();

	void reset();

	~Graph() {
		free(e), free(undeg), free(bounty), free(adj_length);
		for (int i = 0; i < n; i++)
			free(adj[i]);
		free(adj);
	}
};
void Graph::reset() {
	for (int i = 0; i < m; i++)
		e[i].to = e[i].t2, e[i].exists = true;
	memset(bounty, 0, n * sizeof(int));
	R.clear(), A.clear(), V_A.clear();
	ReTestSubgraph.clear(), Q.clear();
	parent.clear(), cur.clear(), dist.clear(), vis.clear();
	// edges_to_be_updated.clear();
	S_alpha.clear(), S_alpha1.clear();
}

void Graph::read_graph_from_dataset(FILE* dataset_file) {
	const int ARRAY_SIZE_INCREMENT = 500000;
	int Emax = ARRAY_SIZE_INCREMENT;

	n = m = 0;
	e = (Edge*)malloc(Emax * sizeof(Edge));

	char line[1000];
	while (fgets(line, 200, dataset_file)) {
		int i = 0;
		e[m].t1 = e[m].t2 = 0;
		while (line[i] < '0' || line[i] > '9') i++;
		while (line[i] >= '0' && line[i] <= '9') e[m].t1 = e[m].t1 * 10 + line[i] - '0', i++;
		while (line[i] < '0' || line[i] >'9') i++;
		while (line[i] >= '0' && line[i] <= '9') e[m].t2 = e[m].t2 * 10 + line[i] - '0', i++;
		n = max(n, max(e[m].t1, e[m].t2)); e[m].exists = true;
		m++;
		e[m].t2 = e[m - 1].t1, e[m].t1 = e[m - 1].t2; e[m].exists = true;
		m++;
		if (m + 1 >= Emax) {
			Emax += ARRAY_SIZE_INCREMENT; // Increment the size of the "e" array. 
			e = (Edge*)realloc(e, Emax * sizeof(Edge));
		}
	}
	fclose(dataset_file);

	n++;

	e = (Edge*)realloc(e, m * sizeof(Edge));
	bounty = (int*)malloc(n * sizeof(int)); memset(bounty, 0, n * sizeof(int));
	undeg = (int*)malloc(n * sizeof(int)); memset(undeg, 0, n * sizeof(int));
	adj = (int**)malloc(n * sizeof(int*));
	adj_length = (int*)malloc(n * sizeof(int));
	undeg = (int*)malloc(n * sizeof(int));

	R.alloc(n); A.alloc(n); ReTestSubgraph.alloc(n); Q.alloc(n); parent.alloc(n); vis.alloc(n); V_A.alloc(n); cur.alloc(n); dist.alloc(n); S_alpha.alloc(n); S_alpha.alloc(n); S_alpha1.alloc(n);

	int* now_edge_num = (int*)malloc(n * sizeof(int));
	for (int i = 0; i < n; i++) adj_length[i] = now_edge_num[i] = 0;
	for (int i = 0; i < m; i++) { adj_length[e[i].t1]++, adj_length[e[i].t2]++; }
	for (int i = 0; i < n; i++) adj[i] = (int*)malloc(adj_length[i] * sizeof(int));
	for (int i = 0; i < m; i++) {
		int t1 = e[i].t1, t2 = e[i].t2;
		adj[t1][now_edge_num[t1]++] = i;
		adj[t2][now_edge_num[t2]++] = i;
	}
	free(now_edge_num);

	for (int i = 0; i < n; i++) undeg[i] = adj_length[i] / 2;

	return;
}
char buf[10000];
bool Graph::read_R_and_A_set(FILE* R_file, FILE* A_file) {
	fgets(buf, 10000, R_file); int len = strlen(buf);
	int pointer = 0;
	while (pointer < len) {
		int tem = 0;
		while ((buf[pointer] < '0' || buf[pointer] > '9') && pointer < len) pointer++;
		if (pointer >= len) break;
		while (buf[pointer] >= '0' && buf[pointer] <= '9') tem = tem * 10 + buf[pointer] - '0', pointer++;
		R.insert(tem);
	}

	buf[0] = '\0';
	fgets(buf, 10000, A_file); len = strlen(buf);
	pointer = 0;
	while (pointer < len) {
		int tem = 0;
		while ((buf[pointer] < '0' || buf[pointer] > '9') && pointer < len) pointer++;
		if (pointer >= len) break;
		while (buf[pointer] >= '0' && buf[pointer] <= '9') tem = tem * 10 + buf[pointer] - '0', pointer++;
		A.insert(tem);
	}

	// Insert many nodes to A, to test the correctness of our algorithms.
	/*for (int u = 1; u <= 10000; u += 200) {
		R.insert(u), A.insert(u);
	}*/

	int edges_in_R = 0;
	for (int i = 0; i < R.size; i++) {
		int u = R.nodes[i];
		for (int j = 0; j < adj_length[u]; j++) {
			Edge& ne = e[adj[u][j]];
			int v = ne.t1 == u ? ne.t2 : ne.t1;
			if (!R.in[v]) continue;
			if (ne.t1 == u)
				edges_in_R++;
		}
	}
	printf("- %-30s: %d\n", "Number of Edges in R", edges_in_R);
	if(edges_in_R<=0){
		return false;
	}
	return true;

	// check(edges_in_R > 0, "E(R) is empty!");
}
void Graph::initialize() {
	for (int i = 0; i < R.size; i++) {
		int u = R.nodes[i];
		V_A.insert(u);
		bounty[u] = undeg[u];
		if (A.in[u])
			bounty[u] += INF;
	}
}
void Graph::get_AADS() {
	int alpha_l = 0, alpha_u = 0; for (int i = 0; i < R.size; i++) alpha_u = max(alpha_u, undeg[R.nodes[i]]);

	// cout<<"alpha_l = " << alpha_l<<endl;
	// cout<<"alpha_u = " << alpha_u<<endl;

	while (alpha_l < alpha_u) {
		int alpha_m = (alpha_l + alpha_u + 1) / 2;
		// cout<<"alpha_m = " << alpha_m<<endl;
		ReTest(alpha_m);
		// cout<<"ReTestSubgraphBounty = " << ReTestSubgraphBounty <<endl;
		if (ReTestSubgraphBounty > (alpha_m - 1) * ReTestSubgraph.size) alpha_l = alpha_m;
		else alpha_u = alpha_m - 1;

		// cout<<"alpha_l = " << alpha_l<<endl;
		// cout<<"alpha_u = " << alpha_u<<endl;
	}
	alpha = alpha_l;
	// cout<<"alpha = " << alpha<<endl;
	if (alpha == 1) {
		printf("Alpha = 1 is too low, contracdicting the assumption of the paper!\n");
		return;
	}

	ReTest(alpha);
	for (int i = 0; i < ReTestSubgraph.size; i++) S_alpha.insert(ReTestSubgraph.nodes[i]);

	return;
}
void Graph::ReTest(int test_value_) {
	test_value = test_value_; pivot = test_value - 1;

	if (test_value == 1) {
		ReTestSubgraph.clear(); for (int i = 0; i < R.size; i++) ReTestSubgraph.insert(R.nodes[i]);
		goto compute;
	}

	while (DinicBFS()) {
		parent.clear(); cur.clear();
		for (int i = 0; i < V_A.size; i++) {
			int u = V_A.nodes[i];
			if (bounty[u] > pivot)
				parent[u] = -2, cur[u] = 0, DinicDFS(u);
		}
	}

	ReTestSubgraph.clear(); Q.clear(); vis.clear();
	for (int i = 0; i < V_A.size; i++) {
		int u = V_A.nodes[i];
		if (bounty[u] > pivot)
			vis.insert(u), Q.push(u), ReTestSubgraph.insert(u);
	}

	while (!Q.empty()) {
		int u = Q.pop();
		for (int i = 0; i < adj_length[u]; i++) {
			Edge& ne = e[adj[u][i]];
			if (!ne.exists) continue;
			if (ne.to != u) continue;
			int from = ne.to == ne.t1 ? ne.t2 : ne.t1;
			if (vis.in[from]) continue;
			vis.insert(from); Q.push(from); ReTestSubgraph.insert(from);
		}
	}

compute:
	int Edgenumber = 0, penalty = 0;
	for (int i = 0; i < ReTestSubgraph.size; i++) {
		int u = ReTestSubgraph.nodes[i];
		if (!R.in[u]) penalty += undeg[u];
		for (int j = 0; j < adj_length[u]; j++) {
			Edge& ne = e[adj[u][j]];
			if (!ne.exists) continue;
			int v = ne.t1 == u ? ne.t2 : ne.t1;
			if (!ReTestSubgraph.in[v]) continue;
			Edgenumber++;
		}
	}
	// cout<<"Edgenumber = " << Edgenumber <<endl;
	// cout<<"penalty = " << penalty <<endl;
	ReTestSubgraphBounty = Edgenumber / 2 - penalty;
	// cout<<"ReTestSubgraphBounty = " << ReTestSubgraphBounty <<endl;
	return;
}
bool Graph::DinicBFS() {
	int dist_t = INF;

	Q.clear(); dist.clear();
	for (int i = 0; i < V_A.size; i++) {
		int u = V_A.nodes[i];
		if (bounty[u] > pivot)
			dist[u] = 1, Q.push(u);
	}

	bool break_loop = false;
	while (!Q.empty()) {
		int u = Q.pop();
		for (int i = 0; i < adj_length[u]; i++) {
			Edge& ne = e[adj[u][i]];
			if (!ne.exists) continue;
			if (ne.to != u) continue;
			int from = ne.t1 == u ? ne.t2 : ne.t1;
			if (bounty[from] < pivot) {
				dist_t = dist[u] + 2; break_loop = true; break;
			}
			if (dist.in[from]) continue;
			dist[from] = dist[u] + 1;
			Q.push(from);
		}
		if (break_loop) break;
	}
	return dist_t != INF;
}
bool Graph::DinicDFS(int u)
{
	if (bounty[u] < pivot) {
		bounty[u]++; bounty[e[parent[u]].to]--; e[parent[u]].to = u;
		if (!V_A.in[u])
			V_A.insert(u);
		return true;
	}
	for (int& i = cur[u]; i < adj_length[u]; i++) {
		Edge& ne = e[adj[u][i]];
		if (!ne.exists) continue;
		if (ne.to != u) continue;
		int from = ne.t1 == u ? ne.t2 : ne.t1;
		if ((!dist.in[from] || dist[from] != dist[u] + 1) && !(bounty[from] < pivot)) continue;
		parent[from] = adj[u][i], cur[from] = 0;
		if (DinicDFS(from)) {
			if (parent[u] == -2) {
				if (bounty[u] == pivot) return true;
				continue;
			}
			bounty[u]++; bounty[e[parent[u]].to]--; e[parent[u]].to = u;
			if (!V_A.in[u])
				V_A.insert(u);
			return true;
		}
	}
	return false;
}
void Graph::get_update_edge() {
	edges_to_be_updated.clear();
	int updated_number = 10000;
	while (edges_to_be_updated.size() < updated_number) {
		edges_to_be_updated.insert((rand() % (m / 2)) * 2);
	}
}

void Graph::read_update_edge(FILE* dataset_file){
	edges_to_be_updated.clear();
	int updated_number = 10000;
	char line[1000];
	fgets(line, 200, dataset_file);
	
	while (fgets(line, 200, dataset_file)) {
		int i = 0;
		int eid=0;
		while (line[i] < '0' || line[i] > '9') i++;
		while (line[i] >= '0' && line[i] <= '9') eid = eid * 10 + line[i] - '0', i++;	
		edges_to_be_updated.insert(eid);
		// cout<<eid<<endl;
	}
	fclose(dataset_file);

}



void Graph::maintain() {
	ReTest(alpha); S_alpha = ReTestSubgraph; S_alpha_bounty = ReTestSubgraphBounty;
	ReTest(alpha + 1); S_alpha1 = ReTestSubgraph; S_alpha1_bounty = ReTestSubgraphBounty;

	//if (true) { //DEBUG
	//	printf("Original\n");
	//	for (int i = 0; i < S_alpha1.size; i++) {
	//		int u = S_alpha1.nodes[i];
	//		printf("u = %d, bounty[u] = %d\n", u, bounty[u]);
	//	}
	//}

	Timer insert_timer, delete_timer; double insert_total_time = 0.0, delete_total_time = 0.0;
	int count_updated_edges = 0;

	int ori_S_alpha_size = S_alpha.size, ori_S_alpha1_size = S_alpha1.size, ori_S_alpha_bounty = S_alpha_bounty, ori_S_alpha1_bounty = S_alpha1_bounty;
	for (auto eid : edges_to_be_updated) {
		count_updated_edges++;

		//if (eid == 48288) // DEBUG
		//	eid = eid;

		//printf("eid = %d\n", eid); // DEBUG
		//if (eid == 53954) { //DEBUG
		//	printf("Before delete\n");
		//	for (int i = 0; i < S_alpha1.size; i++) {
		//		int u = S_alpha1.nodes[i];
		//		printf("u = %d, bounty[u] = %d\n", u, bounty[u]);
		//	}
		//}

		delete_timer.start();
		edge_delete(eid);
		delete_timer.end(); delete_total_time += delete_timer.time();
		allquery_ave_delete_time+=delete_timer.time();
		//if (eid == 53954) { //DEBUG
		//	printf("After delete\n");
		//	for (int i = 0; i < S_alpha1.size; i++) {
		//		int u = S_alpha1.nodes[i];
		//		printf("u = %d, bounty[u] = %d\n", u, bounty[u]);
		//	}
		//}
		//if (eid == 5838) { // DEBUG
		//	printf("After delete\n");
		//	for (int j = 0; j < adj_length[353]; j++) {
		//		Edge& ne = e[adj[353][j]];
		//		int t = ne.t1 == 353 ? ne.t2 : ne.t1;
		//		printf("eid = %d, t = %d, ne.to = %d, S_alpha.in[t] = %d, S_alpha1.in[t] = %d, bounty[t] = %d\n", adj[353][j], t, ne.to, S_alpha.in[t], S_alpha1.in[t], bounty[t]);
		//	}
		//}
		// check_correctness();

		//check_correctness();

		insert_timer.start();
		edge_insert(eid);
		insert_timer.end(); insert_total_time += insert_timer.time();
		allquery_ave_insert_time+=insert_timer.time();
		//if (eid == 53954) { //DEBUG
		//	printf("After Insert\n");
		//	for (int i = 0; i < S_alpha1.size; i++) {
		//		int u = S_alpha1.nodes[i];
		//		printf("u = %d, bounty[u] = %d\n", u, bounty[u]);
		//	}
		//}
		//if (eid == 5838) { // DEBUG
		//	printf("After insert\n");
		//	for (int j = 0; j < adj_length[353]; j++) {
		//		Edge& ne = e[adj[353][j]];
		//		int t = ne.t1 == 353 ? ne.t2 : ne.t1;
		//		printf("eid = %d, t = %d, ne.to = %d, S_alpha.in[t] = %d, S_alpha1.in[t] = %d, bounty[t] = %d\n", adj[353][j], t, ne.to, S_alpha.in[t], S_alpha1.in[t], bounty[t]);
		//	}
		//}

		//check_correctness();

		check(ori_S_alpha_size == S_alpha.size && ori_S_alpha1_size == S_alpha1.size && ori_S_alpha_bounty == S_alpha_bounty && ori_S_alpha1_bounty == S_alpha1_bounty, "original != processed");
	}

	printf("- %-30s: %d\n", "Number of updated edges", count_updated_edges);
	printf("- %-30s: %lf\n", "Dynamic Total insert runtime", insert_total_time);
	printf("- %-30s: %lf\n", "Dynamic Total delete runtime", delete_total_time);
	// printf("- %-30s: %lf\n", "Dynamic Average insert runtime", insert_total_time / count_updated_edges);
	// printf("- %-30s: %lf\n", "Dynamic Average delete runtime", delete_total_time / count_updated_edges);
}
void Graph::edge_delete(int eid) {
	check(e[eid].exists && e[eid + 1].exists, "The edge to be deleted does not exist");
	int u = e[eid].t1, v = e[eid].t2;
	undeg[u]--, undeg[v]--;
	if (!R.in[u]) IncBounty(u); if (!R.in[v]) IncBounty(v);
	DecBounty(e[eid].to, eid);
	DecBounty(e[eid + 1].to, eid + 1);

	if (bounty[u] > 0 && !V_A.in[u]) V_A.insert(u);
	if (bounty[v] > 0 && !V_A.in[v]) V_A.insert(v);

	if (S_alpha1_bounty > alpha * S_alpha1.size) {
		S_alpha = S_alpha1, S_alpha_bounty = S_alpha1_bounty;
		ReTest(alpha + 2); S_alpha1 = ReTestSubgraph, S_alpha1_bounty = ReTestSubgraphBounty;
		alpha++;
	}
	else if (S_alpha_bounty <= (alpha - 1) * S_alpha.size) {
		S_alpha1 = S_alpha, S_alpha1_bounty = S_alpha_bounty;
		ReTest(alpha - 1); S_alpha = ReTestSubgraph, S_alpha_bounty = ReTestSubgraphBounty;
		alpha--;
	}
}
void Graph::edge_insert(int eid) {
	check(!e[eid].exists && !e[eid + 1].exists, "The edge to be inserted exists");
	int u = e[eid].t1, v = e[eid].t2;
	undeg[u]++, undeg[v]++;
	if (!R.in[u]) DecBounty(u, -1); if (!R.in[v]) DecBounty(v, -1);

	int x, y;
	if ((S_alpha1.in[u] && !S_alpha1.in[v]) || (S_alpha.in[u] && !S_alpha.in[v])) x = u, y = v;
	else x = v, y = u;
	e[eid].exists = true, e[eid].to = y;
	IncBounty(y);

	if ((S_alpha1.in[u] && !S_alpha1.in[v]) || (S_alpha.in[u] && !S_alpha.in[v])) x = u, y = v;
	else x = v, y = u;
	e[eid + 1].exists = true, e[eid + 1].to = y;
	IncBounty(y);

	if (bounty[u] > 0 && !V_A.in[u]) V_A.insert(u);
	if (bounty[v] > 0 && !V_A.in[v]) V_A.insert(v);

	if (S_alpha1_bounty > alpha * S_alpha1.size) {
		S_alpha = S_alpha1, S_alpha_bounty = S_alpha1_bounty;
		ReTest(alpha + 2); S_alpha1 = ReTestSubgraph, S_alpha1_bounty = ReTestSubgraphBounty;
		alpha++;
	}
	else if (S_alpha_bounty <= (alpha - 1) * S_alpha.size) {
		S_alpha1 = S_alpha, S_alpha1_bounty = S_alpha_bounty;
		ReTest(alpha - 1); S_alpha = ReTestSubgraph, S_alpha_bounty = ReTestSubgraphBounty;
		alpha--;
	}
}
void Graph::IncBounty(int x) {
	bounty[x]++; if (S_alpha.in[x]) S_alpha_bounty++; if (S_alpha1.in[x]) S_alpha1_bounty++;
	if (bounty[x] == alpha && !S_alpha.in[x]) {
		Q.clear(), parent.clear(); Q.push(x), parent[x] = -1; int y = -1; bool break_loop = false;
		while (!Q.empty()) {
			int u = Q.pop();
			for (int j = 0; j < adj_length[u]; j++) {
				Edge& ne = e[adj[u][j]]; if (ne.to != u || !ne.exists) continue; int from = ne.to == ne.t1 ? ne.t2 : ne.t1;
				if (S_alpha.in[from] || parent.in[from]) continue; parent[from] = adj[u][j];
				if (bounty[from] <= alpha - 2) { y = from, break_loop = true; break; }
				Q.push(from);
			}
			if (break_loop) break;
		}
		if (y != -1) { int now = y; while (now != x) swap(e[parent[now]].to, now); bounty[x]--, bounty[y]++; }
		else {
			for (int i = 0; i < parent.size; i++) check(!S_alpha.in[parent.nodes[i]], "S_alpha error"), S_alpha.insert(parent.nodes[i]);
			S_alpha_bounty = count_bounty_in_set(S_alpha);
		}
	}
	else if (bounty[x] == alpha + 1 && S_alpha.in[x] && !S_alpha1.in[x]) {
		Q.clear(), parent.clear(); Q.push(x), parent[x] = -1; int y = -1; bool break_loop = false;
		while (!Q.empty()) {
			int u = Q.pop();
			for (int j = 0; j < adj_length[u]; j++) {
				Edge& ne = e[adj[u][j]]; if (ne.to != u || !ne.exists) continue; int from = ne.to == ne.t1 ? ne.t2 : ne.t1;
				if (S_alpha1.in[from] || parent.in[from]) continue; parent[from] = adj[u][j];
				check(S_alpha.in[from], "S_alpha error V \\ S_alpha");
				if (bounty[from] <= alpha - 1) { y = from, break_loop = true; break; }
				Q.push(from);
			}
			if (break_loop) break;
		}
		if (y != -1) { int now = y; while (now != x) swap(e[parent[now]].to, now); bounty[x]--, bounty[y]++; }
		else {
			for (int i = 0; i < parent.size; i++) check(!S_alpha1.in[parent.nodes[i]], "S_alpha error"), S_alpha1.insert(parent.nodes[i]);
			S_alpha1_bounty = count_bounty_in_set(S_alpha1);
		}
	}
}
void Graph::DecBounty(int x, int eid) {
	bounty[x]--; if (S_alpha.in[x]) S_alpha_bounty--; if (S_alpha1.in[x]) S_alpha1_bounty--;
	if (S_alpha1.in[x]) {
		if (bounty[x] == alpha - 1) {
			Q.clear(), parent.clear(); Q.push(x), parent[x] = -1; int y = -1; bool break_loop = false;
			while (!Q.empty()) {
				int u = Q.pop();
				for (int j = 0; j < adj_length[u]; j++) {
					Edge& ne = e[adj[u][j]]; if (ne.to == u || !ne.exists) continue;
					if (!S_alpha1.in[ne.to] || parent.in[ne.to]) continue; parent[ne.to] = adj[u][j];
					if (bounty[ne.to] >= alpha + 1) { y = ne.to, break_loop = true; break; }
					Q.push(ne.to);
				}
				if (break_loop) break;
			}
			check(y != -1, "can not find y");
			int now = y;
			while (now != x) { Edge& ne = e[parent[now]]; int from = ne.to == ne.t1 ? ne.t2 : ne.t1; ne.to = now = from; }
			bounty[x]++, bounty[y]--;
		}
		if (eid != -1) e[eid].exists = false;
		Q.clear(), vis.clear(); for (int i = 0; i < S_alpha1.size; i++) if (bounty[S_alpha1.nodes[i]] >= alpha + 1) Q.push(S_alpha1.nodes[i]), vis.insert(S_alpha1.nodes[i]);
		while (!Q.empty()) {
			int u = Q.pop();
			for (int j = 0; j < adj_length[u]; j++) {
				Edge& ne = e[adj[u][j]]; if (ne.to != u || !ne.exists) continue; int from = ne.to == ne.t1 ? ne.t2 : ne.t1;
				if (vis.in[from]) continue; check(S_alpha1.in[from], "S_alpha1 error");
				Q.push(from), vis.insert(from);
			}
		}
		S_alpha1 = vis; S_alpha1_bounty = count_bounty_in_set(S_alpha1);
	}
	else if (S_alpha.in[x] && !S_alpha1.in[x]) {
		if (bounty[x] == alpha - 2) {
			Q.clear(), parent.clear(); Q.push(x), parent[x] = -1; int y = -1; bool break_loop = false;
			while (!Q.empty()) {
				int u = Q.pop();
				for (int j = 0; j < adj_length[u]; j++) {
					Edge& ne = e[adj[u][j]]; if (ne.to == u || !ne.exists) continue;
					if (!S_alpha.in[ne.to] || parent.in[ne.to]) continue; parent[ne.to] = adj[u][j];
					check(!S_alpha1.in[ne.to], "E_\\times(S_alpha, S_alpha1) error");
					if (bounty[ne.to] >= alpha) { y = ne.to, break_loop = true; break; }
					Q.push(ne.to);
				}
				if (break_loop) break;
			}
			check(y != -1, "can not find y");
			int now = y;
			while (now != x) { Edge& ne = e[parent[now]]; int from = ne.to == ne.t1 ? ne.t2 : ne.t1; ne.to = now = from; }
			bounty[x]++, bounty[y]--;
		}
		if (eid != -1) e[eid].exists = false;
		Q.clear(), vis.clear(); for (int i = 0; i < S_alpha.size; i++) if (!S_alpha1.in[S_alpha.nodes[i]] && bounty[S_alpha.nodes[i]] >= alpha) Q.push(S_alpha.nodes[i]), vis.insert(S_alpha.nodes[i]);
		while (!Q.empty()) {
			int u = Q.pop();
			for (int j = 0; j < adj_length[u]; j++) {
				Edge& ne = e[adj[u][j]]; if (ne.to != u || !ne.exists) continue; int from = ne.to == ne.t1 ? ne.t2 : ne.t1;
				if (vis.in[from] || S_alpha1.in[from]) continue; check(S_alpha.in[from], "S_alpha error");
				Q.push(from), vis.insert(from);
			}
		}
		S_alpha = S_alpha1;
		for (int i = 0; i < vis.size; i++) check(!S_alpha.in[vis.nodes[i]], "S_alpha and vis error"), S_alpha.insert(vis.nodes[i]);
		S_alpha_bounty = count_bounty_in_set(S_alpha);
	}
	else {
		if (eid != -1) e[eid].exists = false;
	}
}
int Graph::count_bounty_in_set(Set<int>& S) {
	int count = 0;
	for (int i = 0; i < S.size; i++) {
		int u = S.nodes[i];
		if (A.in[u]) count += bounty[u] - INF;
		else count += bounty[u];
	}
	return count;
}
void Graph::check_correctness() {
	int numerator = 0;
	for (int i = 0; i < S_alpha.size; i++) {
		int u = S_alpha.nodes[i];
		if (!R.in[u])
			numerator -= undeg[u];
		for (int j = 0; j < adj_length[u]; j++) {
			Edge& ne = e[adj[u][j]];
			if (!ne.exists) continue;
			int v = ne.t1 == u ? ne.t2 : ne.t1;
			if (S_alpha.in[u] && !S_alpha.in[v])
				check(ne.to == v, "E_\\times(S_alpha, V \\ S_alpha) are not all point to V \\ S_alpha");
			if (S_alpha1.in[u] && !S_alpha1.in[v])
				check(ne.to == v, "E_\\times(S_alpha1, V \\ S_alpha1) are not all point to V \\ S_alpha1");
			if (S_alpha1.in[v] && !S_alpha1.in[u])
				check(ne.to == u, "E_\\times(S_alpha1, V \\ S_alpha1) are not all point to V \\ S_alpha1");
			if (S_alpha.in[v] && ne.t1 == u)
				numerator++;
		}
	}
	check(numerator > (alpha - 1) * S_alpha.size && numerator <= alpha * S_alpha.size, "round-up density of S_alpha does not equal to alpha");
	for (int i = 0; i < S_alpha1.size; i++) check(S_alpha.in[S_alpha1.nodes[i]], "S_alpha1 is not contained in S_alpha");

	Q.clear(), vis.clear();
	for (int i = 0; i < V_A.size; i++) if (bounty[V_A.nodes[i]] >= alpha + 1) Q.push(V_A.nodes[i]), vis.insert(V_A.nodes[i]);
	while (!Q.empty()) {
		int u = Q.pop();
		for (int j = 0; j < adj_length[u]; j++) {
			Edge& ne = e[adj[u][j]]; if (!ne.exists) continue;
			int from = ne.to == ne.t1 ? ne.t2 : ne.t1; if (vis.in[from]) continue;
			vis.insert(from), Q.push(from);
		}
	}
	check(vis.size == S_alpha1.size, "S_alpha1 does comply with its definition");
	for (int i = 0; i < vis.size; i++) check(S_alpha1.in[vis.nodes[i]], "S_alpha1 does comply with its definition");

	Q.clear(), vis.clear();
	for (int i = 0; i < V_A.size; i++) if (bounty[V_A.nodes[i]] >= alpha) Q.push(V_A.nodes[i]), vis.insert(V_A.nodes[i]);
	while (!Q.empty()) {
		int u = Q.pop();
		for (int j = 0; j < adj_length[u]; j++) {
			Edge& ne = e[adj[u][j]]; if (!ne.exists) continue;
			int from = ne.to == ne.t1 ? ne.t2 : ne.t1; if (vis.in[from]) continue;
			vis.insert(from), Q.push(from);
		}
	}
	check(vis.size == S_alpha.size, "S_alpha does comply with its definition");
	for (int i = 0; i < S_alpha.size; i++) check(vis.in[S_alpha.nodes[i]], "S_alpha does comply with its definition");
}
void Graph::output_AADS() {
	if (true) { // output density of S_alpha
		int edge_number = 0, penalty = 0;
		for (int i = 0; i < S_alpha.size; i++) {
			int u = S_alpha.nodes[i];
			if (!R.in[u])
				penalty += undeg[u];
			for (int j = 0; j < adj_length[u]; j++) {
				Edge& ne = e[adj[u][j]];
				if (!ne.exists) continue;
				int v = ne.t1 == u ? ne.t2 : ne.t1;
				if (S_alpha.in[v] && ne.t1 == u)
					edge_number++;
			}
		}
		printf("- %-30s: (%d - %d) / %d = %lf\n", "R-density", edge_number, penalty, S_alpha.size, double(edge_number - penalty) / S_alpha.size);
	}

	if (true) { // output nodes in S_alpha
		sort(S_alpha.nodes, S_alpha.nodes + S_alpha.size);
		printf("- %-30s: ", "Nodes in S_alpha");
		for (int i = 0; i < S_alpha.size; i++) {
			printf("%d ", S_alpha.nodes[i]);
		}
		printf("\n");
	}
}

int main(int argc, char** argv) {
	srand(16);
	check(argc == 5, "The number of arguments of main are incorrect. ");
	// strcpy(argv[1], "../AnchoredGraphs/amazon/dataset.txt");
	// strcpy(argv[2], "../AnchoredGraphs/amazon/R.txt");
	// strcpy(argv[3], "../AnchoredGraphs/amazon/A.txt");
	char dataset_address[1000]; strcpy(dataset_address, argv[1]);

	printf("----------Now processing %s----------\n", dataset_address);
	printf("- %-30s: %s\n", "Algorithm used", "BasicIns and BasicDel");

	FILE* R_file = fopen(argv[2], "r"), * A_file = fopen(argv[3], "r");
	check(R_file != NULL && A_file != NULL, "Can not open R or A sets file.");
	int number_of_cases = read_number(R_file); read_number(A_file);

	FILE* dataset_file = fopen(dataset_address, "r");
	check(dataset_file != NULL, "Can not open file dataset_address\n");


	FILE* updatededge_file = fopen(argv[4], "r");
	check(updatededge_file != NULL, "Can not open file updatededge_file\n");

	Graph G; Timer timer;
	G.read_graph_from_dataset(dataset_file);
	G.read_update_edge(updatededge_file);

	printf("- %-30s: %d, %d\n", "n, m", G.n, G.m);
	int nnow_case=0;
	for (int now_case = 0; now_case < number_of_cases; now_case++) {
		if (now_case == 10) // DEBUG
			now_case = now_case;

		printf("-----Now case %d-----\n", now_case);
		G.reset();
		if(G.read_R_and_A_set(R_file, A_file)==false){
			continue;
		}

		// G.read_R_and_A_set(R_file, A_file);

		timer.start();
		G.initialize();
		G.get_AADS(); if (G.alpha == 1) { printf("This case should be skipped.\n"); continue; }
		timer.end(); printf("- %-30s: %lf\n", "Static algorithm runtime", timer.time());

		// G.output_AADS();

		// G.get_update_edge();
		G.maintain();
		nnow_case++;

		// if(nnow_case==10){
		// 	break;
		// }

	}
	fclose(R_file), fclose(A_file);

	printf("- %-30s: %lf\n", "Dynamic Average insert runtime", G.allquery_ave_delete_time / nnow_case);
	printf("- %-30s: %lf\n", "Dynamic Average delete runtime", G.allquery_ave_insert_time / nnow_case);

	return 0;
}
