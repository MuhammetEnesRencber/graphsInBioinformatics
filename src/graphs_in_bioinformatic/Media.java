package graphs_in_bioinformatic;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class Media {
    private DirectedGraph graph;
    
    public Media() {
        graph = new DirectedGraph();
    }
    
    public void dataLoad(double threshold) {
        String filePath = "C:\\Users\\rencb\\OneDrive\\Masaüstü\\9606.protein.info.v12.0.txt";
        String edgePath = "C:\\Users\\rencb\\OneDrive\\Masaüstü\\9606.protein.links.v12.0.txt";
        
        graph = new DirectedGraph();
        
        double normalizedThreshold = threshold;
        if (threshold >= 1) {
            normalizedThreshold = threshold / 1000.0;
        }
        
       
        
        Set<String> proteinsWithValidEdges = new HashSet<>();
        
        try (BufferedReader br = new BufferedReader(new FileReader(edgePath))) {
            String line;
            int totalEdges = 0;
            int filteredEdges = 0;
            
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#") || line.startsWith("protein1") || line.trim().isEmpty())
                    continue;
                
                totalEdges++;
                String[] parts = line.split(" ");
                
                String id1 = parts[0];
                String id2 = parts[1];
                int rawScore = Integer.parseInt(parts[2]);
                double score = rawScore / 1000.0;
                
                if (score >= normalizedThreshold) {
                    proteinsWithValidEdges.add(id1);
                    proteinsWithValidEdges.add(id2);
                    filteredEdges++;
                }
            }
            
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }
        
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            int loadedVertices = 0;
            
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#") || line.trim().isEmpty()) continue;
                
                String[] p = line.split("\t");
                String id = p[0];
                
                if (proteinsWithValidEdges.contains(id)) {
                    String name = p[1];
                    int size = Integer.parseInt(p[2]);
                    String annot = p[3];
                    
                    graph.addVertex(id, name, size, annot);
                    loadedVertices++;
                }
            }
            
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }
        
        try (BufferedReader br = new BufferedReader(new FileReader(edgePath))) {
            String line;
            int addedEdges = 0;
            
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#") || line.startsWith("protein1") || line.trim().isEmpty())
                    continue;
                
                String[] parts = line.split(" ");
                
                String id1 = parts[0];
                String id2 = parts[1];
                int rawScore = Integer.parseInt(parts[2]);
                double score = rawScore / 1000.0;
                
                if (score >= normalizedThreshold) {
                    if (graph.getVertex(id1) != null && graph.getVertex(id2) != null) {
                        graph.addEdge(id1, id2, rawScore);
                        addedEdges++;
                    }
                }
            }
            
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void searchProtein(String id) {
        if (graph == null || graph.size() == 0) {
            System.out.println("GRAF BOŞ! Önce dataLoad() metodunu çağırın.");
            return;
        }
        
        Vertex protein = graph.getVertex(id);
        if (protein == null) {
            System.out.println("Protein bulunamadı: " + id);
            return;
        }
        
        System.out.println("ID: " + protein.getName());
        System.out.println("İsim: " + protein.getPreferredName());
        System.out.println("Boyut: " + protein.getProteinSize());
        System.out.println("Açıklama: " + protein.getAnnotation());
    }
    
    public void interaction(String id1, String id2) {
        Vertex protein1 = graph.getVertex(id1);
        Vertex protein2 = graph.getVertex(id2);
        
        if (protein1 == null) {
            System.out.println("Protein1 bulunamadı: " + id1);
            return;
        }
        if (protein2 == null) {
            System.out.println("Protein2 bulunamadı: " + id2);
            return;
        }
        
        if (protein1.hasEdge(id2)) {
            System.out.println(id1 + " → " + id2 + " ETKİLEŞİMİ VAR");
        } else {
            System.out.println(id1 + " → " + id2 + " ETKİLEŞİMİ YOK");
        }
        
        if (protein2.hasEdge(id1)) {
            System.out.println(id2 + " → " + id1 + " ETKİLEŞİMİ VAR");
        }
    }
    
    //Doğru Path Bulma Metodu
    public void findCorrectPaths(String id1, String id2) {
        if (graph == null || graph.size() == 0) {
            System.out.println("Graf boş! Önce dataLoad() çağırın.");
            return;
        }

        Vertex v1 = graph.getVertex(id1);
        Vertex v2 = graph.getVertex(id2);
        
        if (v1 == null || v2 == null) {
            System.out.println("Proteinler bulunamadı!");
            return;
        }
        
        System.out.println("===PATH ANALİZİ===");
        System.out.println("Başlangıç: " + id1 + " (" + v1.getPreferredName() + ")");
        System.out.println("Hedef: " + id2 + " (" + v2.getPreferredName() + ")");
        System.out.println();
        
        // Tüm shortest path'leri al
        List<List<Vertex>> allShortestPaths = graph.findAllShortestPaths(id1, id2);
        
        if (allShortestPaths == null || allShortestPaths.isEmpty()) {
            System.out.println("❌ Yol bulunamadı!");
            return;
        }
        
        // 1. SHORTEST PATH 
        System.out.println("1. SHORTEST PATH ");
        System.out.println("-".repeat(50));
        
        List<Vertex> shortestPath = allShortestPaths.get(0); // İlk shortest path
        System.out.println("Cost = " + (shortestPath.size() - 1) );
        for (int i = 0; i < shortestPath.size(); i++) {
            Vertex v = shortestPath.get(i);
            System.out.println((i + 1) + "- " + v.getName());
        }
        
        System.out.println();
        
        // 2. MOST CONFIDENT PATH 
        System.out.println("2. MOST CONFIDENT PATH ");
        System.out.println("-".repeat(50));
        
        List<Vertex> mostConfidentPath = null;
        double maxScore = -1;
        
        for (List<Vertex> path : allShortestPaths) {
            double totalScore = 0;
            for (int i = 0; i < path.size() - 1; i++) {
                Vertex current = path.get(i);
                Vertex next = path.get(i + 1);
                
                for (Edge edge : current.getEdges()) {
                    if (edge.getDestination() == next) {
                        totalScore += edge.getWeight();
                        break;
                    }
                }
            }
            
            if (totalScore > maxScore) {
                maxScore = totalScore;
                mostConfidentPath = path;
            }
        }
        
        System.out.println("Cost = " + maxScore);
        for (int i = 0; i < mostConfidentPath.size(); i++) {
            Vertex v = mostConfidentPath.get(i);
            System.out.println((i + 1) + "- " + v.getName());
        }
        
        System.out.println();
        
     // 3. CHEAPEST PATH 
        System.out.println("3. CHEAPEST PATH ");
        System.out.println("-".repeat(50));
        
        List<Vertex> cheapestPath = null;
        double minScore = Double.MAX_VALUE;
        
        for (List<Vertex> path : allShortestPaths) {
            double totalScore = 0;
            for (int i = 0; i < path.size() - 1; i++) {
                Vertex current = path.get(i);
                Vertex next = path.get(i + 1);
                
                for (Edge edge : current.getEdges()) {
                    if (edge.getDestination() == next) {
                        totalScore += edge.getWeight();
                        break;
                    }
                }
            }
            
            if (totalScore < minScore) {
                minScore = totalScore;
                cheapestPath = path;
            }
        }
        
        System.out.println("Cost = " + minScore);
        for (int i = 0; i < cheapestPath.size(); i++) {
            Vertex v = cheapestPath.get(i);
            System.out.println((i + 1) + "- " + v.getName());
           }
        }
        
     
    
    
    private List<Vertex> dijkstraForCheapest(String id1, String id2) {
		// TODO Auto-generated method stub
		return null;
	}

	public void calculateAndPrintMetrics() {
        if (graph == null || graph.size() == 0) {
            System.out.println("❌ Graf boş! Önce dataLoad() çağırın.");
            return;
        }
        
        System.out.println("=== GRAF METRIKLERI ===");
        
        int vertexCount = graph.size();
        System.out.println("1. VertexCount: " + vertexCount);
        
        int edgeCount = 0;
        Set<String> countedEdges = new HashSet<>();
        
        for (Vertex vertex : graph.vertices()) {
            for (Edge edge : vertex.getEdges()) {
                String edgeKey = vertex.getName() + "->" + edge.getDestination().getName();
                
                if (!countedEdges.contains(edgeKey)) {
                    countedEdges.add(edgeKey);
                    edgeCount++;
                }
            }
        }
        
        System.out.println("2. EdgeCount: " + edgeCount);
        
        double averageDegree = (vertexCount > 0) ? (double) edgeCount / vertexCount : 0;
        System.out.printf("3. Average Degree: %.2f\n", averageDegree);
        
        // Reciprocity hesapla
        int reciprocalPairs = 0;
        for (Vertex vertex : graph.vertices()) {
            String sourceId = vertex.getName();
            for (Edge edge : vertex.getEdges()) {
                String targetId = edge.getDestination().getName();
                Vertex targetVertex = graph.getVertex(targetId);
                if (targetVertex != null && targetVertex.hasEdge(sourceId)) {
                    reciprocalPairs++;
                }
            }
        }
        double reciprocity = (edgeCount > 0) ? (double) reciprocalPairs / (2 * edgeCount) : 0;
        System.out.printf("4. Reciprocity: %.4f\n", reciprocity);
        
        // Diameter 
        System.out.println("5. Diameter (yaklaşık):");
        int maxReachable = 0;
        String farthestVertex = "";
        int testCount = Math.min(10, vertexCount);
        int tested = 0;

        for (Vertex start : graph.vertices()) {
            Queue<String> bfsResult = graph.getBreadthFirstTraversal(start.getName());
            int reachableCount = bfsResult.size();
            
            if (reachableCount > maxReachable) {
                maxReachable = reachableCount;
                farthestVertex = start.getName();
            }
            
            tested++;
            if (tested >= testCount) break;
        }
        
        System.out.printf("   En fazla ulaşılabilen vertex: %d (başlangıç: %s)\n", 
                          maxReachable, farthestVertex);
    }
    
    public void traverseProteinNetwork(String originProteinId) {
        if (graph == null || graph.size() == 0) {
            System.out.println("Graf boş! Önce dataLoad() çağırın.");
            return;
        }
        
        Vertex origin = graph.getVertex(originProteinId);
        if (origin == null) {
            System.out.println("Protein bulunamadı: " + originProteinId);
            return;
        }
        
        Queue<String> bfsOrder = graph.getBreadthFirstTraversal(originProteinId);
        Queue<String> dfsOrder = graph.getDepthFirstTraversal(originProteinId);
        
        System.out.println("\n=== BREADTH-FIRST TRAVERSAL (BFS) ===");
        System.out.println("Başlangıç: " + originProteinId + " (" + origin.getPreferredName() + ")");
        System.out.println("=".repeat(60));
        System.out.println("BFS Sırası (Toplam " + bfsOrder.size() + " protein):");
        
        int count = 1;
        int maxShow = 30;
        for (String proteinId : bfsOrder) {
            if (count > maxShow) {
                System.out.println("... ve " + (bfsOrder.size() - maxShow) + " protein daha");
                break;
            }
            
            Vertex protein = graph.getVertex(proteinId);
            System.out.printf("%3d. %s", count, proteinId);
            
            if (!protein.getPreferredName().isEmpty()) {
                System.out.printf(" - %s", protein.getPreferredName());
            }
            
            System.out.printf(" [%d etkileşim]", protein.getEdges().size());
            
            if (count == 1) {
                System.out.print(" ⭐ (BAŞLANGIÇ)");
            }
            
            System.out.println();
            count++;
        }
        
        System.out.println("\n" + "=".repeat(60) + "\n");
        
        System.out.println("=== DEPTH-FIRST TRAVERSAL (DFS) ===");
        System.out.println("Başlangıç: " + originProteinId + " (" + origin.getPreferredName() + ")");
        System.out.println("=".repeat(60));
        System.out.println("DFS Sırası (Toplam " + dfsOrder.size() + " protein):");
        
        count = 1;
        for (String proteinId : dfsOrder) {
            if (count > maxShow) {
                System.out.println("... ve " + (dfsOrder.size() - maxShow) + " protein daha");
                break;
            }
            
            Vertex protein = graph.getVertex(proteinId);
            System.out.printf("%3d. %s", count, proteinId);
            
            if (!protein.getPreferredName().isEmpty()) {
                System.out.printf(" - %s", protein.getPreferredName());
            }
            
            System.out.printf(" [%d etkileşim]", protein.getEdges().size());
            
            if (count == 1) {
                System.out.print(" ⭐ (BAŞLANGIÇ)");
            }
            
            System.out.println();
            count++;
        }
        
        System.out.println("\n=== KARŞILAŞTIRMA ===");
        System.out.println("BFS ile ulaşılan: " + bfsOrder.size() + " protein");
        System.out.println("DFS ile ulaşılan: " + dfsOrder.size() + " protein");
        
        if (bfsOrder.size() == dfsOrder.size()) {
            System.out.println("Her iki traversal da AYNI sayıda protein buldu");
        } else {
            System.out.println("Traversallar FARKLI sayıda protein buldu!");
        }
    }
    
    public void menu() {
        Scanner scan = new Scanner(System.in);

        while (true) {
            System.out.println("\n=== PROTEIN ETKİLEŞİM AĞI MENÜSÜ ===");
            System.out.println("1. Proteinleri Yükle (Threshold ile)");
            System.out.println("2. Protein Ara");
            System.out.println("3. İki Protein Arası Etkileşim Sorgula");
            System.out.println("4. Shortest Path,Cheapest Path,Most Path");
            System.out.println("5. Graf Metriklerini Hesapla");
            System.out.println("6. BFS ve DFS Traversal");
            System.out.println("0. Çıkış");
            System.out.print("Seçiminiz: ");

            int choice = scan.nextInt();
            scan.nextLine();

            switch (choice) {
                case 1:
                    System.out.print("Threshold değeri girin (0-1000 arası): ");
                    int rawThreshold = scan.nextInt();
                    
                    if (rawThreshold < 0 || rawThreshold > 1000) {
                        System.out.println("Hata: Threshold 0-1000 arasında olmalı!");
                        break;
                    }
                    
                    double normalizedThreshold = rawThreshold / 1000.0;
                    System.out.println("Veri yükleniyor (Threshold: " + normalizedThreshold + ")...");
                    
                    long startTime = System.currentTimeMillis();
                    dataLoad(normalizedThreshold);
                    long loadTime = System.currentTimeMillis() - startTime;
                    
                    System.out.println("Graf yüklendi!");
                    System.out.println("  Yükleme süresi: " + loadTime + " ms");
                    System.out.println("  Vertex Count: " + graph.size());
                    
                    int edgeCount = 0;
                    Set<String> countedEdges = new HashSet<>();
                    for (Vertex vertex : graph.vertices()) {
                        for (Edge edge : vertex.getEdges()) {
                            String edgeKey = vertex.getName() + "->" + edge.getDestination().getName();
                            if (!countedEdges.contains(edgeKey)) {
                                countedEdges.add(edgeKey);
                                edgeCount++;
                            }
                        }
                    }
                    
                    System.out.println("   Edge Count: " + edgeCount);
                    break;
                    
                case 2:
                    System.out.print("Protein ID girin: ");
                    String id = scan.nextLine().trim();
                    searchProtein(id);
                    break;

                case 3:
                    System.out.print("Birinci Protein ID: ");
                    String p1 = scan.nextLine().trim();
                    System.out.print("İkinci Protein ID: ");
                    String p2 = scan.nextLine().trim();
                    interaction(p1, p2);
                    break;

                case 4:
                	System.out.print("Başlangıç Protein ID: ");
                    String startId = scan.nextLine().trim();
                    System.out.print("Hedef Protein ID: ");
                    String endId = scan.nextLine().trim();
                   
                        findCorrectPaths(startId, endId);
                    
                    break;
                    
                    
                case 5:
                    calculateAndPrintMetrics();
                    break;
                    
                case 6:
                    System.out.print("Başlangıç Protein ID: ");
                    String originId = scan.nextLine().trim();
                    traverseProteinNetwork(originId);
                    break;
                    
                case 0:
                    System.out.println("Programdan çıkılıyor...");
                    scan.close();
                    return;
                    
                default:
                    System.out.println(" Geçersiz seçim! Lütfen 0-6 arası bir sayı girin.");
            }
            
            if (choice != 0) {
                System.out.print("\nDevam etmek için ENTER'a basın...");
                scan.nextLine();
            }
        }
    }
    
    
    
    
    
    
    public void debugCheapestPath(String id1, String id2) {
        System.out.println("=== DEBUG CHEAPEST PATH ===");
        
        Vertex start = graph.getVertex(id1);
        Vertex end = graph.getVertex(id2);
        
        if (start == null || end == null) {
            System.out.println("Proteinler bulunamadı!");
            return;
        }
        
        System.out.println("Başlangıç: " + id1 + " (" + start.getPreferredName() + ")");
        System.out.println("Hedef: " + id2 + " (" + end.getPreferredName() + ")");
        System.out.println();
        
        // 1. BFS ile bağlantı var mı kontrol et
        System.out.println("1. BFS ile bağlantı kontrolü:");
        Queue<String> bfs = graph.getBreadthFirstTraversal(id1);
        boolean reachable = false;
        for (String proteinId : bfs) {
            if (proteinId.equals(id2)) {
                reachable = true;
                break;
            }
        }
        System.out.println("   BFS ile ulaşılabilir: " + reachable);
        System.out.println("   BFS ulaştığı protein sayısı: " + bfs.size());
        
        // 2. Basit DFS ile path bulmayı dene
        System.out.println("\n2. Basit DFS ile path bulma:");
        List<Vertex> simplePath = findPathSimpleDFS(start, end, new HashSet<>(), new ArrayList<>());
        if (simplePath != null) {
            System.out.println("   Basit DFS ile path bulundu!");
            System.out.print("   Path: ");
            for (Vertex v : simplePath) {
                System.out.print(v.getName() + " -> ");
            }
            System.out.println("END");
            
            // Score'u hesapla
            double totalScore = 0;
            for (int i = 0; i < simplePath.size() - 1; i++) {
                Vertex current = simplePath.get(i);
                Vertex next = simplePath.get(i + 1);
                for (Edge edge : current.getEdges()) {
                    if (edge.getDestination() == next) {
                        totalScore += edge.getWeight();
                        break;
                    }
                }
            }
            System.out.println("   Toplam score: " + totalScore);
        } else {
            System.out.println("   Basit DFS ile de yol bulunamadı!");
        }
        
        // 3. Tüm olası path'leri bul (küçük bir DFS)
        System.out.println("\n3. Tüm olası path'ler (max 1000):");
        List<List<Vertex>> allPaths = new ArrayList<>();
        findAllPathsDFS(start, end, new HashSet<>(), new ArrayList<>(), allPaths, 1000);
        System.out.println("   Bulunan path sayısı: " + allPaths.size());
        
        if (!allPaths.isEmpty()) {
            // En düşük score'lu path'i bul
            List<Vertex> cheapest = null;
            double minScore = Double.MAX_VALUE;
            
            for (List<Vertex> path : allPaths) {
                double totalScore = 0;
                for (int i = 0; i < path.size() - 1; i++) {
                    Vertex current = path.get(i);
                    Vertex next = path.get(i + 1);
                    for (Edge edge : current.getEdges()) {
                        if (edge.getDestination() == next) {
                            totalScore += edge.getWeight();
                            break;
                        }
                    }
                }
                
                if (totalScore < minScore) {
                    minScore = totalScore;
                    cheapest = path;
                }
            }
            
            System.out.println("   En düşük score: " + minScore);
            System.out.print("   Cheapest path: ");
            for (Vertex v : cheapest) {
                System.out.print(v.getName() + " -> ");
            }
            System.out.println("END");
        }
    }

    // Basit DFS path bulma
    private List<Vertex> findPathSimpleDFS(Vertex current, Vertex end, 
                                          Set<Vertex> visited, List<Vertex> path) {
        visited.add(current);
        path.add(current);
        
        if (current == end) {
            return new ArrayList<>(path);
        }
        
        for (Edge edge : current.getEdges()) {
            Vertex neighbor = edge.getDestination();
            if (!visited.contains(neighbor)) {
                List<Vertex> result = findPathSimpleDFS(neighbor, end, visited, path);
                if (result != null) {
                    return result;
                }
            }
        }
        
        path.remove(path.size() - 1);
        visited.remove(current);
        return null;
    }

    // Tüm path'leri bulan DFS
    private void findAllPathsDFS(Vertex current, Vertex end, 
                               Set<Vertex> visited, List<Vertex> path,
                               List<List<Vertex>> allPaths, int maxPaths) {
        if (allPaths.size() >= maxPaths) return;
        
        visited.add(current);
        path.add(current);
        
        if (current == end) {
            allPaths.add(new ArrayList<>(path));
        } else {
            for (Edge edge : current.getEdges()) {
                Vertex neighbor = edge.getDestination();
                if (!visited.contains(neighbor)) {
                    findAllPathsDFS(neighbor, end, visited, path, allPaths, maxPaths);
                }
            }
        }
        
        path.remove(path.size() - 1);
        visited.remove(current);
    }
    
    
    
    
    
    
    
    
}