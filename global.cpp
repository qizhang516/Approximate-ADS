#include <bits/stdc++.h>
#include <fstream>
using namespace std;
typedef int type;

const type INF = 2000000000;
double timer;
inline void timer_start() { timer = clock(); }
inline double timer_end() { double t = (clock() - timer) / CLOCKS_PER_SEC; return t; }


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
};
// Edge structure, t1 and t2 are two endpoints, "to" is the endpoint the edge points to. 
struct Edge { type t1, t2, to; };


vector<Set> RSet;
vector<Set> ASet;
vector<bool> RSetflag;
type Rsettruecnt=0;


// Most computation is in the Graph structure. 
struct Graph {
	Graph() {}
	void read_edge(char* address); // Read the graph file, input the address of the file. 
	void read_R(char* address); // Read the R nodes set. 
	void read_A(char* address); // Read the A nodes set. 
	void get_densest(); // Get the nearly densest anchored subgraph using re-orientation network. 

	Edge* e; // Edges
	type n, m; // Number of nodes and edges. 
	type* d; // Indegree plus weight of nodes. 
	type** adj; // The adjacent edges of nodes. 
	type* adj_length; // The length of the adjacent list of every single node, i.e., doubled degree of nodes. 

	Set R, A; // R and A nodes set. 
	Set NearDensest; // The nearly densest anchored subgraph we obtain. 

	Queue Q; // A queue used in BFS. 

	type test_value, pivot; // The test value used in the re-orientation network, pivot = test_value - 1. 
	void ReTest(type test_value);
	Set ReTestSubgraph; // The dense subgraph that ReTest has found. 
	type ReTestSubgraphEdge;
	type* dist; // dist array denote "distance", used in BFS to store the distance between source and nodes. Besides, dist array also additionally replace the "visit" array, where dist[u] = -1 if u has not been visited. 
	type* cur; // used in DinicDFS, functioning as a well-known optimizing method of Dinic algorithm. 
	type* p; // p array denote "parent", used in DFS to reverse paths, storing the edge index of the BFS search tree. Besides, p array also additionally replace the "visit" array, where p[u] = -1 if u has not been visited. 
	Set s_node; // Used in re-orientation, contain all nodes directly connect to source (i.e. d[u] < pivot). 
	bool DinicBFS(); // Return whether source s can reach sink t. 
	bool DinicDFS(type); // Return whether this recursion has found an augment path. 

	void output(char* output_file_address);
	void output();
	void display_all_edges();

	double mybd=0.0;
	double getlbfromSetR();

	type get_max_d(); // A simple function that returns the max indegree. 
	bool check(); // Test whether the result graph is right, if so, output True. 

	~Graph() {
		free(e); free(d); free(adj); free(adj_length); free(cur); free(dist); free(p);
	}
};


Graph G1;

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
		e[m].to = e[m].t2;
		m++;
		e[m].t2 = e[m - 1].t1, e[m].t1 = e[m - 1].t2;
		e[m].to = e[m].t2;
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
	d = (type*)malloc(n * sizeof(type));
	dist = (type*)malloc(n * sizeof(type));
	cur = (type*)malloc(n * sizeof(type));
	p = (type*)malloc(n * sizeof(type));
	adj = (type**)malloc(n * sizeof(type*));
	adj_length = (type*)malloc(n * sizeof(type));

	//mem of d, dist, cur, p

	// Allocate memory for our array-implemented sets and queue. 
	R.alloc(n); A.alloc(n); NearDensest.alloc(n); ReTestSubgraph.alloc(n); s_node.alloc(n); Q.alloc(n);

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

void read_RSet(char* address, type num) {
	type querytimes=0;
	fstream finA(address, std::ios::in);
	char buff[4194304];
	finA >> querytimes;
	// cout<< querytimes<<endl;
	RSet.resize(querytimes);
	RSetflag.resize(querytimes, 1);
	for(int i=0; i<querytimes; i++){
		RSet[i].alloc(num); 
	}
	// cout<<"emmmm"<<endl;
	for(int i=0; i<querytimes; i++){
		finA >> buff;
		// std::cout<<buff<<"\n";
		char delims[]=",";
		char *res=strtok(buff, delims);
		if(res){
			// std::cout<<"res="<<res<<std::endl;
			// std::cout<<"res="<<stoi(res)<<std::endl;
			RSet[i].insert(std::stoi(res));
		}
		while(res=strtok(NULL, delims)){
			// std::cout<<"res="<<res<<std::endl;
			RSet[i].insert(std::stoi(res));
		}
				// std::cout<<"check R"<<i<<std::endl;
		// at least one edge check
		bool flag=false;
		for(int j=0; j<RSet[i].size; j++)
		{
			for(int k=0; k<RSet[i].size; k++)
			{
				if(j<k)
				{
					type u=RSet[i].nodes[j];
					type v=RSet[i].nodes[k];
					for(int x=0; x<G1.m; x++)
					{
						Edge& ne = G1.e[x];
						if( (ne.t1==u && ne.t2==v) || (ne.t1==v && ne.t2==u)){
							flag=true;
							break;
						}
					}
				}
				if(flag){
					break;
				}
			}
			if(flag){
				break;
			}
		}
		if(!flag){
			RSetflag[i]=0;
			Rsettruecnt++;
		}
		// std::cout<<"R"<<i<<"\t"<<flag<<std::endl;
		// break;
	}
	finA.close();
}



void read_ASet(char* address, type num) {

	type querytimes=0;
	fstream finA(address, std::ios::in);
	char buff[4194304];
	finA >> querytimes;
	// cout<< querytimes<<endl;
	ASet.resize(querytimes);
	for(int i=0; i<querytimes; i++){
		ASet[i].alloc(num); 
	}
	for(int i=0; i<querytimes; i++){
		finA >> buff;
		// std::cout<<buff<<"\n";
		char delims[]=",";
		char *res=strtok(buff, delims);
		if(res){
			// std::cout<<"res="<<res<<std::endl;
			ASet[i].insert(std::stoi(res));
		}
		while(res=strtok(NULL, delims)){
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
	// Initiate the d array. 
	for (type i = 0; i < n; i++) {
		if (A.in[i]) d[i] = adj_length[i] / 2 + INF;
		else if (R.in[i]) d[i] = adj_length[i] / 2;
		else d[i] = 0;
	}

	type du, dl;
	// First get the binary search range du and dl, the answer is in (dl, du]. 
	dl = du = 0;
	for (type i = 0; i < R.size; i++) du = max(du, adj_length[R.nodes[i]] / 2);


	//dl
	dl = ceil(getlbfromSetR());
	// cout<<dl<<endl;
	// Start binary search. 
	while (dl < du) {
		type dm = (dl + du + 1) / 2;
		ReTest(dm);
		// if (ReTestSubgraph.size != A.size) dl = dm;
		if (ReTestSubgraphEdge > (dm - 1) * ReTestSubgraph.size) dl = dm;
		else du = dm - 1;
	}
	type ans = dl;

	// Now we have the ans, which equals to the round-up value of the nearly densest anchored subgraph. 
	// Next invoke ReTest(ans) once and the ReTestSubgraph is the NearDensest we want. 
	ReTest(ans);
	NearDensest.insert(ReTestSubgraph);
	// cout<<ans<<endl;

	return;
}

void Graph::ReTest(type test_value_) {
	test_value = test_value_; pivot = test_value - 1;
	// cout<<pivot<<endl;

	while (DinicBFS()) {
		for (type i = 0; i < s_node.size; i++) {
			DinicDFS(s_node.nodes[i]);
		}
	}

	// Now, the ReTestSubgraph equals to all nodes reachable to a node u, where d[u] >= pivot + 1. Next we compute it. 
	ReTestSubgraph.clear(); Q.clear();
	memset(dist, 0, n * sizeof(type));
	for (type i = 0; i < n; i++)
		if (d[i] >= pivot + 1)
			dist[i] = 1, Q.push(i), ReTestSubgraph.insert(i);

	while (!Q.empty()) { // Start from nodes with d[] > pivot, traversing reversely. 
		type u = Q.pop();
		for (type i = 0; i < adj_length[u]; i++) {
			Edge& ne = e[adj[u][i]];
			if (ne.to != u) continue;
			type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
			if (dist[from]) continue;
			dist[from] = 1;
			Q.push(from); ReTestSubgraph.insert(from);
		}
	}

	// Compute the number of edges in ReTestSubgraph, for computing its density.
	type Edgenumber = 0, penalty = 0;
	for (type i = 0; i < ReTestSubgraph.size; i++) {
		type u = ReTestSubgraph.nodes[i];
		if (!R.in[u]) penalty += adj_length[u] / 2;
		for (type j = 0; j < adj_length[u]; j++) {
			Edge& ne = e[adj[u][j]];
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
	for (type i = 0; i < n; i++) if (d[i] < pivot) s_node.insert(i);

	type dist_t = INF; // The distance from source s to sink t. 

	// Initialize queue and some arrays. 
	Q.clear();
	for (type i = 0; i < n; i++) p[i] = dist[i] = -1, cur[i] = 0;
	for (type i = 0; i < s_node.size; i++) dist[s_node.nodes[i]] = 1, Q.push(s_node.nodes[i]);

	// Start traversal. 
	while (!Q.empty()) {
		type u = Q.pop();
		// if (dist[u] >= dist_t) break;
		for (type i = 0; i < adj_length[u]; i++) {
			Edge& ne = e[adj[u][i]];
			if (ne.to == u) continue;
			if (dist[ne.to] != -1) continue;
			dist[ne.to] = dist[u] + 1;
			if (d[ne.to] > pivot && dist_t == INF) dist_t = dist[ne.to] + 1;
			Q.push(ne.to);
		}
	}
	return dist_t != INF;
}

bool Graph::DinicDFS(type u)
{
	if (d[u] > pivot) {
		type from = e[p[u]].t1 == e[p[u]].to ? e[p[u]].t2 : e[p[u]].t1;
		d[from]++; d[u]--; e[p[u]].to = from;
		return true;
	}
	for (type& i = cur[u]; i < adj_length[u]; i++) {
		Edge& ne = e[adj[u][i]];
		if (ne.to == u) continue;
		if (dist[ne.to] != dist[u] + 1) continue;
		p[ne.to] = adj[u][i];
		if (DinicDFS(ne.to)) {
			if (s_node.in[u]) {
				if (d[u] == pivot) return true; // This means arc (s, u) has been saturated, so we can return. 
				continue; // This means arc (s, u) has not been saturated, so we continue to search for augment path from u. 
			}
			type from = e[p[u]].t1 == e[p[u]].to ? e[p[u]].t2 : e[p[u]].t1;
			d[from]++; d[u]--; e[p[u]].to = from;
			return true;
		}
	}
	return false;
}


double Graph::getlbfromSetR(){

	int sn=R.size;
	int sm=0;
	for(int i=0; i<m; i++){
		Edge& ne = e[i];
		if(R.in[ne.t1] && R.in[ne.t2]){
			sm++;
		}
	}

	// printf("sn:%ld, sm:%ld\n", sn, sm);
	double mybd=(float) sm/sn;
 	// printf("mybd:%lf\n", mybd);
	return mybd;
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


void Graph::output() {
	// Customize what you want to output. 
	const bool OUTPUT_ANCHORED_DENSITY = true;
	const bool OUTPUT_NODES = true;
	
	if (OUTPUT_ANCHORED_DENSITY) {
		printf("Anchored Subgraph: \n");
		printf("Number of Nodes: %d. \n", NearDensest.size);
		type numerator = 0;
		for (type i = 0; i < NearDensest.size; i++)
			numerator += A.in[NearDensest.nodes[i]] ? d[NearDensest.nodes[i]] - INF : d[NearDensest.nodes[i]];
		type edge_cnt = numerator;
		for (type i = 0; i < NearDensest.size; i++)
			if (!R.in[NearDensest.nodes[i]])
				edge_cnt += adj_length[NearDensest.nodes[i]] / 2;

		printf("Number of Edges: %d * 2 = %d. \n", edge_cnt / 2, edge_cnt);
		printf("%lf\t", 1.0 * numerator / NearDensest.size);

		printf("R-density:%lf\n", 1.0 * numerator / NearDensest.size);
	}

	if (OUTPUT_NODES) {
		printf("Near Densest Anchored Subgraph: ");
		sort(NearDensest.nodes, NearDensest.nodes + NearDensest.size);
		for (type i = 0; i < NearDensest.size; i++)	printf("%d ", NearDensest.nodes[i]);
		printf("\n");
	}

	// long long memory = m * sizeof(Edge) + n * sizeof(type) + m * sizeof(type) + n * sizeof(type) + 7 * n * (sizeof(type) + sizeof(bool)) + n * sizeof(type) + 3 * n * sizeof(type);
	//                 e                  d                  adj                adj_length         7 * Set                                 Q                  dist+cur+p

	// printf("\t%lf\n", 1.0 * memory / 1024 / 1024);

	// printf("Memory: %lf MB = %lf GB\n", 1.0 * memory / 1024 / 1024, 1.0 * memory / 1024 / 1024 / 1024);

}

void Graph::display_all_edges() {
	printf("All edges:------------------------\n");
	for (type u = 0; u < n; u++)
		for (type j = 0; j < adj_length[u]; j++) {
			Edge& ne = e[adj[u][j]];
			if (ne.to != u) continue;
			type from = ne.to == ne.t1 ? ne.t2 : ne.t1;
			printf("%d %d\n", from, ne.to);
		}
	printf("----------------------------------\n");
}
bool Graph::check() {
	// This function check the correctness of NearDensest, denoted as Ans for simplicity. 
	// 1. Check all edge between Ans and (V - Ans) should point to (V - Ans);
	// 2. Check the sum of d[] in Ans should equal to (2|E(Ans)| - \sum degree) (the definition of anchored density);
	// 3. Note that this function does not check whether Ans is the densest subgraph, if you want to check, please check if the result is the same as using the global algorithm.

	// 1. Check all edge between Ans and (V - Ans) should point to (V - Ans);
	for (type i = 0; i < NearDensest.size; i++) {
		type u = NearDensest.nodes[i];
		for (type j = 0; j < adj_length[u]; j++) {
			Edge& ne = e[adj[u][j]];
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
	for (type i = 0; i < NearDensest.size; i++){
		if (!R.in[NearDensest.nodes[i]]){
			numerator -= adj_length[NearDensest.nodes[i]] / 2; // numerator -= \sum degree
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
	char adjacent_list_address[100], R_set_address[100], A_set_address[100];
	if (argc == 4) {
		strcpy(adjacent_list_address, argv[1]); strcpy(R_set_address, argv[2]); strcpy(A_set_address, argv[3]);
	}
	else {
		printf("Usage:\n");
		printf("./main adjacent_list_address R_set_address A_set_address\n");
		return 0;
	}
	double runtime;

	timer_start();
	G1.read_edge(adjacent_list_address);
	runtime = timer_end();
	printf("Reading edges runtime:\t\t%.2lf seconds\n", runtime);

	// G.read_R(R_set_address);
	// G.read_A(A_set_address);

	read_RSet(R_set_address, G1.n);
	read_ASet(A_set_address, G1.n);
    Graph *G;
	double total_runtime=0.0; 
	type realcnt=0;

	for(type i=0; i<RSet.size(); i++){
		if(RSetflag[i]==0){
			continue;
		}
		realcnt++;
		G=new Graph();
		G->read_edge(adjacent_list_address);
		for(type j=0; j<RSet[i].size; j++){
			G->R.insert(RSet[i].nodes[j]);
		}
		for(type j=0; j<ASet[i].size; j++){
			G->A.insert(ASet[i].nodes[j]);
		}
		timer_start();
		G->get_densest();
		runtime = timer_end();
		total_runtime+=runtime;
		G->output();
		// printf("query %d \t time %lf\n", i, 1.0 * runtime);
		// printf("query %d \t avetime %lf\n", i, 1.0 * total_runtime/realcnt);
		if(G->check()==false){
			printf("the %d query time check error\n", i);
			delete G;
			continue;
		}
        delete G;
	}
	cout<<argv[1]; 
	printf(" average time:\t%.6f seconds\n", total_runtime / (RSet.size()-Rsettruecnt));
	return 0;
}
