#pragma once
#ifndef MESH_GRID
#define MESH_GRID
#include "mesh_2d.h"
#include "configure.h"
class MeshGrid :public Mesh2D
{
public:
	MeshGrid(const int _cols, const int _rows);

	const std::vector<Point2>& getVertices()const;
	const std::vector<Edge>& getEdges()const;
	const std::vector<Indices>& getPolygonsIndices()const;
	const std::vector<Indices>& getPolygonsNeighbors()const;
	const std::vector<Indices>& getPolygonsEdges()const;
	const std::vector<Indices>& getVertexStructures()const;
	const std::vector<Indices>& getEdgeStructures()const;
	const std::vector<Indices>& getTriangulationIndices()const;
	const int getPolygonVerticesCounts()const;
	const std::vector<int>& getBoundaryVertexIndices()const;/* clockwise order */
	const std::vector<int>& getBoundaryEdgeIndices()const;

	InterpolateVertex getInterpolateVertex(const cv::Point_<float>& p)const;
	InterpolateVertex getInterpolateVertex(const cv::Point_<double>& p)const;

	template<typename T>
	InterpolateVertex getInterpolateVertexTemplate(const cv::Point_<T>& p)const;
private:
};


#endif // !MESH_GRID
