package graphs_in_bioinformatic;

import java.util.*;

public class DirectedGraph {
    private HashMap<String, Vertex> vertices;

    public DirectedGraph() {
        this.vertices = new HashMap<>();
    }

    public void addEdge(String source, String destination, int weight) {
        Vertex source_v = vertices.get(source);
        Vertex destination_v = vertices.get(destination);

        if (source_v != null && destination_v != null && source_v.hasEdge(destination)) {
            System.out.println("This edge has already added!");
        } else {
            if (vertices.get(source) == null) {
                source_v = new Vertex(source);
                vertices.put(source, source_v);
            }

            if (vertices.get(destination) == null) {
                destination_v = new Vertex(destination);
                vertices.put(destination, destination_v);
            }

            Edge edge = new Edge(source_v, destination_v, weight);
            source_v.addEdge(edge);
        }
    }

    public void addVertex(String id, String preferredName, int size, String annotation) {
        if (!vertices.containsKey(id)) {
            Vertex v = new Vertex(id, preferredName, size, annotation);
            vertices.put(id, v);
        }
    }

    public void print() {
        for (Vertex v : vertices.values()) {
            System.out.print(v.getName() + " -> ");
            Iterator<Vertex> neighbors = v.getNeighborIterator();
            while (neighbors.hasNext()) {
                Vertex n = neighbors.next();
                System.out.print(n.getName() + " ");
            }
            System.out.println();
        }
    }

    public Iterable<Vertex> vertices() {
        return vertices.values();
    }

    public int size() {
        return vertices.size();
    }

    private void resetVertices() {
        for (Vertex v : vertices.values()) {
            v.unvisit();
            v.setCost(0);
            v.setParent(null);
        }
    }

    public Queue<String> getBreadthFirstTraversal(String origin) {
        resetVertices();
        Queue<String> traversalOrder = new LinkedList<>();
        Queue<Vertex> vertexQueue = new LinkedList<>();

        Vertex originVertex = vertices.get(origin);
        originVertex.visit();

        traversalOrder.add(origin);
        vertexQueue.add(originVertex);

        while (!vertexQueue.isEmpty()) {
            Vertex frontVertex = vertexQueue.remove();
            Iterator<Vertex> neighbors = frontVertex.getNeighborIterator();

            while (neighbors.hasNext()) {
                Vertex nextNeighbor = neighbors.next();
                if (!nextNeighbor.isVisited()) {
                    nextNeighbor.visit();
                    traversalOrder.add(nextNeighbor.getName());
                    vertexQueue.add(nextNeighbor);
                }
            }
        }

        return traversalOrder;
    }

    public Queue<String> getDepthFirstTraversal(String origin) {
        resetVertices();
        Queue<String> traversalOrder = new LinkedList<>();
        Stack<Vertex> vertexStack = new Stack<>();

        Vertex originVertex = vertices.get(origin);
        if (originVertex == null) {
            return traversalOrder;
        }

        originVertex.visit();
        vertexStack.push(originVertex);

        while (!vertexStack.isEmpty()) {
            Vertex currentVertex = vertexStack.pop();
            traversalOrder.add(currentVertex.getName());

            Iterator<Vertex> neighbors = currentVertex.getNeighborIterator();
            Stack<Vertex> tempStack = new Stack<>();

            while (neighbors.hasNext()) {
                Vertex neighbor = neighbors.next();
                if (!neighbor.isVisited()) {
                    tempStack.push(neighbor);
                }
            }

            while (!tempStack.isEmpty()) {
                Vertex neighbor = tempStack.pop();
                neighbor.visit();
                vertexStack.push(neighbor);
            }
        }

        return traversalOrder;
    }

    public Vertex getVertex(String id) {
        return vertices.get(id);
    }

    public boolean containsVertex(String id) {
        return vertices.containsKey(id);
    }

    public Collection<Vertex> getAllVertices() {
        return vertices.values();
    }

   

    //  ALL SHORTEST PATHS (Tüm kısa yolları bul)
    public List<List<Vertex>> findAllShortestPaths(String startId, String endId) {
        if (!vertices.containsKey(startId) || !vertices.containsKey(endId)) {
            return null;
        }

        resetVertices();
        Map<Vertex, Integer> distances = new HashMap<>();
        Map<Vertex, List<Vertex>> predecessors = new HashMap<>();
        Queue<Vertex> queue = new LinkedList<>();

        Vertex start = vertices.get(startId);
        Vertex end = vertices.get(endId);

        for (Vertex v : vertices.values()) {
            distances.put(v, Integer.MAX_VALUE);
            predecessors.put(v, new ArrayList<>());
        }

        distances.put(start, 0);
        queue.add(start);

        while (!queue.isEmpty()) {
            Vertex current = queue.remove();
            int currentDist = distances.get(current);

            Iterator<Vertex> neighbors = current.getNeighborIterator();
            while (neighbors.hasNext()) {
                Vertex neighbor = neighbors.next();

                if (distances.get(neighbor) > currentDist + 1) {
                    distances.put(neighbor, currentDist + 1);
                    predecessors.get(neighbor).clear();
                    predecessors.get(neighbor).add(current);
                    queue.add(neighbor);
                } else if (distances.get(neighbor) == currentDist + 1) {
                    predecessors.get(neighbor).add(current);
                }
            }
        }

        List<List<Vertex>> allPaths = new ArrayList<>();
        if (distances.get(end) != Integer.MAX_VALUE) {
            findAllPathsDFS(end, start, predecessors, new ArrayList<>(), allPaths);
        }

        return allPaths;
    }

    private void findAllPathsDFS(Vertex current, Vertex start,
                               Map<Vertex, List<Vertex>> predecessors,
                               List<Vertex> currentPath,
                               List<List<Vertex>> allPaths) {
        currentPath.add(0, current);

        if (current == start) {
            allPaths.add(new ArrayList<>(currentPath));
        } else {
            for (Vertex pred : predecessors.get(current)) {
                findAllPathsDFS(pred, start, predecessors, currentPath, allPaths);
            }
        }

        currentPath.remove(0);
    }

   
    

    //Dijkstra 
  /*  public ArrayList<Vertex> dijkstra(String startId, String endId) {
        Vertex start = getVertex(startId);
        Vertex goal = getVertex(endId);

        if (start == null || goal == null) return null;

        for (Vertex v : vertices.values()) {
            v.setCost(Double.POSITIVE_INFINITY);
            v.setParent(null);
            v.unvisit();
        }

        start.setCost(0);

        PriorityQueue<Vertex> pq = new PriorityQueue<>(
            (a, b) -> Double.compare(a.getCost(), b.getCost())
        );

        pq.add(start);

        while (!pq.isEmpty()) {
            Vertex current = pq.poll();
            current.visit();

            if (current == goal) break;

            for (Edge e : current.getEdges()) {
                Vertex n = e.getDestination();
                int score = e.getWeight();
                
                // Yüksek skor = düşük maliyet için
                double w = 1000 - score; // Ters çevir
                double newCost = current.getCost() + w;

                if (newCost < n.getCost()) {
                    n.setCost(newCost);
                    n.setParent(current);
                    
                    pq.remove(n);
                    pq.add(n);
                }
            }
        }

        ArrayList<Vertex> path = new ArrayList<>();
        Vertex step = goal;

        if (step.getParent() == null) return null;

        while (step != null) {
            path.add(0, step);
            step = step.getParent();
        }

        return path;
    }*/
}