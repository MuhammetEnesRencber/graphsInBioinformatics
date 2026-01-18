package graphs_in_bioinformatic;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class Vertex {
	
	private String preferredName;
	private int proteinSize;
	private String annotation;
	
	private String name;
	private ArrayList<Edge> edges;
	private Vertex parent;
	private boolean visited;  
	private double cost;  

	public Vertex(String name,String prefferedName, int proteinSize, String annotation) {
		this.name = name;
		this.preferredName = prefferedName;
		this.proteinSize = proteinSize;
		this.annotation = annotation;
		edges = new ArrayList<Edge>();
		parent = null;
		visited = false;
	}
	
	public Vertex(String id) {
	    this(id, "", 0, "");
	}

	public void addEdge(Edge e) {
		edges.add(e);
	}

	public ArrayList<Edge> getEdges() {
		return this.edges;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public Vertex getParent() {
		return parent;
	}

	public void setParent(Vertex parent) {
		this.parent = parent;
	}

	public double getCost() {
		return cost;
	}

	public void setCost(double cost) {
		this.cost = cost;
	}

	public void visit() {
		this.visited = true;
	}

	public void unvisit() {
		this.visited = false;
	}

	public boolean isVisited() {
		return this.visited;
	}

	public Vertex getUnvisitedNeighbor() {
		Vertex result = null;

		Iterator<Vertex> neighbors = getNeighborIterator();
		while (neighbors.hasNext() && (result == null))
		{
			Vertex nextNeighbor = neighbors.next();
			if (!nextNeighbor.isVisited())
				result = nextNeighbor;
		} // end while

			return result;
	}

	public boolean hasEdge(String neighbor) {
		boolean found = false;
		Iterator<Vertex> neighbors = getNeighborIterator();
		while (neighbors.hasNext())
		{
			Vertex nextNeighbor = neighbors.next();
			if (nextNeighbor.getName().equalsIgnoreCase(neighbor))
			{
				found = true;
				break;
			}
		} 

		return found;
	}

	public Iterator<Vertex> getNeighborIterator()
	{
		return new NeighborIterator();
	} // end getNeighborIterator

	private class NeighborIterator implements Iterator<Vertex>
	{
		int edgeIndex = 0;  
		private NeighborIterator()
		{
			edgeIndex = 0; 
		} // end default constructor

		public boolean hasNext()
		{
			return edgeIndex < edges.size();
		} // end hasNext

		public Vertex next()
		{
			Vertex nextNeighbor = null;

			if (hasNext())
			{
				nextNeighbor = edges.get(edgeIndex).getDestination();
				edgeIndex++;
			}
			else
				throw new NoSuchElementException();

			return nextNeighbor;
		} // end next

		public void remove()
		{
			throw new UnsupportedOperationException();
		} // end remove
	} // end NeighborIterator

	public int getProteinSize() {
		// TODO Auto-generated method stub
		return proteinSize;
	}

	public String getPreferredName() {
		return preferredName;
	}

	public void setPreferredName(String preferredName) {
		this.preferredName = preferredName;
	}

	public String getAnnotation() {
		return annotation;
	}

	public void setAnnotation(String annotation) {
		this.annotation = annotation;
	}

	public void setProteinSize(int proteinSize) {
		this.proteinSize = proteinSize;
	}
}

