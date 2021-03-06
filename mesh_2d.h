#pragma once
#ifndef MESH_2D
#define MESH_2D
#include <vector>
#include "configure.h"
const int EDGE_VERTEX_SIZE = 2;

class Edge
{
public:
	Edge(const int e1, const int e2)
	{
		indices[0] = e1;
		indices[1] = e2;
	}
	int indices[EDGE_VERTEX_SIZE];
	
private:
};


class Indices
{
public:
	std::vector<int> indices;

	Indices(){}

	Indices(const int i0, const int i1, const int i2)
	{
		indices.emplace_back(i0);
		indices.emplace_back(i1);
		indices.emplace_back(i2);
	}

	Indices(const int i0, const int i1, const int i2,const int i3)
	{
		indices.emplace_back(i0);
		indices.emplace_back(i1);
		indices.emplace_back(i2);
		indices.emplace_back(i3);
	}

private:
};

class InterpolateVertex
{
public:
	int polygon;
	std::vector<double> weights;

	InterpolateVertex()
	{
		polygon = -1;
	}

	InterpolateVertex(const InterpolateVertex& v)
	{
		polygon = v.polygon;
		weights = v.weights;
	}

	InterpolateVertex(const int _polygon, const int w0, const int w1, const int w2)
	{
		polygon = _polygon;
		weights.emplace_back(w0);
		weights.emplace_back(w1);
		weights.emplace_back(w2);
	}

	InterpolateVertex(const int _polygon, std::vector<double>& _weight)
	{
		polygon = _polygon;
		weights = _weight;
	}
private:
};

class Mesh2D
{
public:
	int nw, nh;
	double lw, lh;
	Mesh2D(const int _cols, const int _rows);
	virtual ~Mesh2D();
	virtual const std::vector<Point2>& getVertices()const = 0;
	virtual const std::vector<Edge>& getEdges()const = 0;
	virtual const std::vector<Indices>& getPolygonsIndices()const = 0;
	virtual const std::vector<Indices>& getPolygonsNeighbors()const = 0;
	virtual const std::vector<Indices>& getPolygonsEdges()const = 0;
	virtual const std::vector<Indices>& getVertexStructures()const = 0;
	virtual const std::vector<Indices>& getEdgeStructures()const = 0;
	virtual const std::vector<Indices>& getTriangulationIndices()const = 0;
	virtual const int getPolygonVerticesCounts()const = 0;
	virtual const std::vector<int>& getBoundaryVertexIndices()const = 0;/* clockwise order */
	virtual const std::vector<int>& getBoundaryEdgeIndices()const = 0;

	virtual InterpolateVertex getInterpolateVertex(const cv::Point_<float>& p)const = 0;
	virtual InterpolateVertex getInterpolateVertex(const cv::Point_<double>& p)const = 0;

	virtual const std::vector<Point2>& getPolygonsCenter()const;

	template<typename T>
	int getGridIndexOfPoint(const cv::Point_<T>& p)const;

protected:
	mutable std::vector<Point2> m_vertices;
	mutable std::vector<Point2> m_polygon_center;
	mutable std::vector<Edge> m_edges;
	mutable std::vector<Indices> m_polygons_indices;
	mutable std::vector<Indices> m_polygons_neighbors;
	mutable std::vector<Indices> m_polygons_edges;
	mutable std::vector<Indices> m_vertex_structures;
	mutable std::vector<Indices> m_edge_structures;
	mutable std::vector<Indices> m_triangulation_indices;
	mutable std::vector<int> m_boundary_vertex_indices;
	mutable std::vector<int> m_boundary_edge_indices;
};


#endif // !MESH_2D
