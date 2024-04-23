#ifndef GRAPHLIBRARY_H
#define GRAPHLIBRARY_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <algorithm>
#include <fstream>

class Graph
{
private:
    int V;
    bool isDirected;
    bool isWeighted;
    std::vector<std::vector<int>> adj;
    std::vector<std::vector<int>> weight;
    std::unordered_map<int, std::string> vertexNames;
    std::unordered_map<std::string, int> vertexIndices;

public:
    bool isCyclicUtil(int v, std::vector<bool> &visited, std::vector<bool> &recStack)
    {
        if (!visited[v])
        {
            visited[v] = true;
            recStack[v] = true;

            // Check for all vertices adjacent to this vertex
            for (auto i : adj[v])
            {
                if (!visited[i] && isCyclicUtil(i, visited, recStack))
                {
                    return true; // Found cycle
                }
                else if (recStack[i])
                {
                    return true; // Found back edge
                }
            }
        }
        recStack[v] = false; // Remove vertex from recursion stack
        return false;
    }
    bool isCyclic()
    {
        std::vector<bool> visited(V, false);
        std::vector<bool> recStack(V, false);

        for (int i = 0; i < V; i++)
        {
            if (!visited[i] && isCyclicUtil(i, visited, recStack))
            {
                return true;
            }
        }
        return false;
    }

public:
    Graph(int vertices, bool directed = false, bool weighted = false) : V(vertices), isDirected(directed), isWeighted(weighted)
    {
        adj.resize(V);
        if (weighted)
        {
            weight.resize(V, std::vector<int>(V, std::numeric_limits<int>::max()));
        }
    }

    // Add directed edge from u to v with Default Weight 1
    void addDirectedEdge(int u, int v)
    {
        adj[u].push_back(v);
    }

    // Add directed edge from u to v
    void addDirectedEdge(int u, int v, int w = 1)
    {
        adj[u].push_back(v);
        if (isWeighted)
        {
            weight[u][v] = w;
        }
    }

    // Add undirected edge from u to v with Default Weight 1
    void addUndirectedEdge(int u, int v)
    {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    // Add undirected edge between u and v
    void addUndirectedEdge(int u, int v, int w = 1)
    {
        adj[u].push_back(v);
        adj[v].push_back(u);
        if (isWeighted)
        {
            weight[u][v] = w;
            weight[v][u] = w;
        }
    }

    // Add vertex name for better visualization
    void addVertexName(int v, const std::string &name)
    {
        vertexNames[v] = name;
        vertexIndices[name] = v;
    }

    void printGraph()
    {
        // Print column headers
        std::cout << "Vertex";
        if (isWeighted)
        {
            std::cout << "\tConnected To (Weight)";
        }
        else
        {
            std::cout << "\tConnected To";
        }
        std::cout << std::endl;

        // Print each vertex and its connections
        for (int i = 0; i < V; ++i)
        {
            std::cout << "(" << i;
            if (vertexNames.find(i) != vertexNames.end())
            {
                std::cout << ") " << vertexNames[i] << ": ";
            }
            else
            {
                std::cout << "): ";
            }

            for (int v : adj[i])
            {
                std::cout << v;
                if (isWeighted)
                {
                    std::cout << "(" << weight[i][v] << ")";
                }
                std::cout << " ";
            }
            std::cout << std::endl;
        }
    }

    void printVertexConnections()
    {
        // Print column headers
        std::cout << "Vertex";
        if (isWeighted)
        {
            std::cout << "\tConnected To (Weight)";
        }
        else
        {
            std::cout << "\tConnected To";
        }
        std::cout << std::endl;

        // Print each vertex and its connections
        for (int i = 0; i < V; ++i)
        {
            std::cout << "(" << i;
            if (vertexNames.find(i) != vertexNames.end())
            {
                std::cout << ") " << vertexNames[i] << ":";
            }
            else
            {
                std::cout << "):";
            }

            for (int v : adj[i])
            {
                std::cout << " (" << v;
                if (isWeighted)
                {
                    std::cout << ", " << weight[i][v];
                }
                std::cout << ")";
            }
            std::cout << std::endl;
        }
    }
    std::string getVertexName(int v) {
    if (vertexNames.find(v) != vertexNames.end()) {
        return vertexNames[v];
    } else {
        // Return a default value or a placeholder to indicate the name is not available
        return "Vertex " + std::to_string(v);
    }
}


    std::vector<int> getAdjacentVertices(int v)
    {
        return adj[v];
    }

    int getWeight(int u, int v)
    {
        return weight[u][v];
    }

    // Topological Sort for Directed Acyclic Graphs (DAGs)
    std::vector<int> topologicalSort()
    {
        std::vector<int> inDegree(V, 0);
        for (int u = 0; u < V; ++u)
        {
            for (int v : adj[u])
            {
                inDegree[v]++;
            }
        }

        std::queue<int> q;
        for (int u = 0; u < V; ++u)
        {
            if (inDegree[u] == 0)
            {
                q.push(u);
            }
        }

        std::vector<int> topoOrder;
        while (!q.empty())
        {
            int u = q.front();
            q.pop();
            topoOrder.push_back(u);

            for (int v : adj[u])
            {
                if (--inDegree[v] == 0)
                {
                    q.push(v);
                }
            }
        }

        if (topoOrder.size() != V)
        {
            // Graph has a cycle
            return std::vector<int>();
        }

        return topoOrder;
    }

    // Kosaraju's Algorithm for Strongly Connected Components (SCCs)
    std::vector<std::vector<int>> getStronglyConnectedComponents()
    {
        std::vector<std::vector<int>> SCCs;
        std::vector<int> order = topologicalSort();
        std::vector<bool> visited(V, false);

        Graph reverseGraph = getReverseGraph();

        for (int i = V - 1; i >= 0; --i)
        {
            int v = order[i];
            if (!visited[v])
            {
                std::vector<int> SCC;
                reverseGraph.DFSUtil(v, visited, SCC);
                SCCs.push_back(SCC);
            }
        }

        return SCCs;
    }

    // Minimum Spanning Tree (MST) using Kruskal's Algorithm
    std::vector<std::pair<int, std::pair<int, int>>> findMST()
    {
        std::vector<std::pair<int, std::pair<int, int>>> edges;
        for (int u = 0; u < V; ++u)
        {
            for (int v : adj[u])
            {
                if (!isDirected || (isDirected && u < v))
                {
                    edges.push_back({weight[u][v], {u, v}});
                }
            }
        }

        std::sort(edges.begin(), edges.end());

        std::vector<std::pair<int, std::pair<int, int>>> MST;
        std::vector<int> parent(V, -1);
        for (auto edge : edges)
        {
            int w = edge.first;
            int u = edge.second.first;
            int v = edge.second.second;
            int setU = find(parent, u);
            int setV = find(parent, v);
            if (setU != setV)
            {
                MST.push_back({w, {u, v}});
                parent[setU] = setV;
            }
        }

        return MST;
    }

    // Eulerian Path and Cycle
    bool hasEulerianPath()
    {
        if (isDirected)
        {
            int outDegreeDiff = 0;
            for (int u = 0; u < V; ++u)
            {
                int outDegree = adj[u].size();
                if (outDegree == 0)
                {
                    return false;
                }
                if (outDegreeDiff == 0)
                {
                    outDegreeDiff = outDegree - getInDegree(u);
                }
                else if (outDegreeDiff != outDegree - getInDegree(u))
                {
                    return false;
                }
            }
            return outDegreeDiff == 1 || outDegreeDiff == -1;
        }
        else
        {
            int oddDegreeCount = 0;
            for (int u = 0; u < V; ++u)
            {
                int degree = adj[u].size();
                if (degree % 2 != 0)
                {
                    oddDegreeCount++;
                }
            }
            return oddDegreeCount == 0 || oddDegreeCount == 2;
        }
    }

    // Method to find maximum flow from source to sink
    int findMaxFlow(int source, int sink)
    {
        int maxFlow = 0;
        std::vector<std::vector<int>> residualGraph(V, std::vector<int>(V, 0));

        // Initialize residual graph with given capacities in the original graph
        for (int u = 0; u < V; ++u)
        {
            for (int v = 0; v < adj[u].size(); ++v)
            {
                int adjV = adj[u][v];
                residualGraph[u][adjV] = weight[u][adjV]; // assuming weight[u][v] holds capacity
            }
        }

        std::vector<int> parent(V); // Parent array to store augmenting path

        // Augment the flow while there is path from source to sink
        while (int pathFlow = findAugmentedPath(residualGraph, parent, source, sink))
        {
            // Update residual capacities of the edges and reverse edges along the path
            updateResidualGraph(residualGraph, parent, source, sink, pathFlow);

            // Add path flow to overall flow
            maxFlow += pathFlow;
        }

        return maxFlow;
    }

    // Method to find an augmenting path and return path flow
    int findAugmentedPath(std::vector<std::vector<int>> &residualGraph, std::vector<int> &parent, int source, int sink)
    {
        fill(parent.begin(), parent.end(), -1);
        parent[source] = -2; // Source is its own parent
        std::queue<std::pair<int, int>> q;
        q.push({source, INT_MAX});

        while (!q.empty())
        {
            int u = q.front().first;
            int flow = q.front().second;
            q.pop();

            for (int v = 0; v < V; ++v)
            {
                if (u != v && parent[v] == -1 && residualGraph[u][v])
                {
                    parent[v] = u;
                    int new_flow = std::min(flow, residualGraph[u][v]);
                    if (v == sink)
                        return new_flow;
                    q.push({v, new_flow});
                }
            }
        }

        return 0;
    }

    // Method to update the residual graph
    void updateResidualGraph(std::vector<std::vector<int>> &residualGraph, std::vector<int> &parent, int source, int sink, int pathFlow)
    {
        for (int v = sink; v != source; v = parent[v])
        {
            int u = parent[v];
            residualGraph[u][v] -= pathFlow;
            residualGraph[v][u] += pathFlow; // Increase reverse flow
        }
    }

    // Serialize the graph to a file
    void saveToFile(const std::string &filename)
    {
        std::ofstream outfile(filename);
        if (!outfile)
        {
            std::cerr << "Error opening file for writing: " << filename << std::endl;
            return;
        }

        outfile << V << " " << isDirected << " " << isWeighted << std::endl;
        for (int u = 0; u < V; ++u)
        {
            for (int v : adj[u])
            {
                outfile << u << " " << v;
                if (isWeighted)
                {
                    outfile << " " << weight[u][v];
                }
                outfile << std::endl;
            }
        }

        outfile.close();
    }

    // Deserialize the graph from a file
    // Deserialize the graph from a file
void loadFromFile(const std::string &filename)
{
    std::ifstream infile(filename);
    if (!infile)
    {
        std::cerr << "Error opening file for reading: " << filename << std::endl;
        return;
    }

    int numVertices, directed, weighted;
    infile >> numVertices >> directed >> weighted;
    V = numVertices;
    isDirected = (directed == 1);
    isWeighted = (weighted == 1);

    adj.clear();
    adj.resize(V);
    if (isWeighted)
    {
        weight.clear();
        weight.resize(V, std::vector<int>(V, std::numeric_limits<int>::max()));
    }

    int u, v;
    int w;
    while (infile >> u >> v)
    {
        if (isWeighted)
        {
            infile >> w;
        }
        else
        {
            w = 1; // Default weight if not specified
        }
        // Here we check if the graph is directed or not and add the edge accordingly
        if (isDirected)
        {
            addDirectedEdge(u, v, w);
        }
        else
        {
            addUndirectedEdge(u, v, w);
        }
    }

    infile.close();
}


    // Remove a vertex from the graph
    void removeVertex(int v)
    {
        for (int u = 0; u < V; ++u)
        {
            if (u == v)
            {
                adj[u].clear();
            }
            else
            {
                adj[u].erase(std::remove(adj[u].begin(), adj[u].end(), v), adj[u].end());
            }
        }
        adj.erase(adj.begin() + v);

        if (isWeighted)
        {
            weight.erase(weight.begin() + v);
            for (auto &row : weight)
            {
                row.erase(row.begin() + v);
            }
        }

        if (vertexNames.find(v) != vertexNames.end())
        {
            vertexNames.erase(v);
        }
        updateVertexIndices();
        V--;
    }

    // Remove an edge from the graph
    void removeEdge(int u, int v)
    {
        auto it = std::find(adj[u].begin(), adj[u].end(), v);
        if (it != adj[u].end())
        {
            adj[u].erase(it);
        }

        if (isWeighted)
        {
            weight[u][v] = std::numeric_limits<int>::max();
        }

        if (!isDirected)
        {
            it = std::find(adj[v].begin(), adj[v].end(), u);
            if (it != adj[v].end())
            {
                adj[v].erase(it);
            }

            if (isWeighted)
            {
                weight[v][u] = std::numeric_limits<int>::max();
            }
        }
    }


    std::vector<int> dijkstra(int src) {
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;
    std::vector<int> dist(V, std::numeric_limits<int>::max());

    pq.push(std::make_pair(0, src));
    dist[src] = 0;

    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();

        // Traverse all adjacent vertices of u
        for (const auto &i : adj[u]) {
            int v = i; // The vertex label
            int weight = getWeight(u, v);

            // If there's a shorter path to v through u.
            if (dist[v] > dist[u] + weight) {
                // Updating distance of v
                dist[v] = dist[u] + weight;
                pq.push(std::make_pair(dist[v], v));
            }
        }
    }

    return dist;
}


private:
    // Find parent of a set (for Union-Find)
    int find(std::vector<int> &parent, int u)
    {
        if (parent[u] == -1)
        {
            return u;
        }
        return find(parent, parent[u]);
    }

    // Find in-degree of a vertex
    int getInDegree(int v)
    {
        int inDegree = 0;
        for (int u = 0; u < V; ++u)
        {
            if (std::find(adj[u].begin(), adj[u].end(), v) != adj[u].end())
            {
                inDegree++;
            }
        }
        return inDegree;
    }

    // Utility function for DFS
    void DFSUtil(int v, std::vector<bool> &visited, std::vector<int> &result)
    {
        visited[v] = true;
        result.push_back(v);
        for (int u : adj[v])
        {
            if (!visited[u])
            {
                DFSUtil(u, visited, result);
            }
        }
    }

    // Get the reverse of the graph
    Graph getReverseGraph()
    {
        Graph reverseGraph(V, isDirected, isWeighted);
        for (int u = 0; u < V; ++u)
        {
            for (int v : adj[u])
            {
                reverseGraph.addDirectedEdge(v, u, isWeighted ? weight[u][v] : 1);
            }
        }
        return reverseGraph;
    }

    // Update vertex indices after removal
    void updateVertexIndices()
    {
        vertexIndices.clear();
        for (const auto &entry : vertexNames)
        {
            vertexIndices[entry.second] = entry.first;
        }
    }
};

#endif