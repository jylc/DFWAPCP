#include "mesh_grid.h"

MeshGrid::MeshGrid(const int _cols, const int _rows) :Mesh2D(_cols, _rows) {}

const std::vector<Point2>& MeshGrid::getVertices()const
{
	if (m_vertices.empty())
	{
		const int memory = (nw + 1) * (nh + 1);
		m_vertices.reserve(memory);
		for (int h = 0; h <= nh; ++h)
		{
			for (int w = 0; w <= nw; ++w)
			{
				m_vertices.emplace_back(w * lw, h * lh);
			}
		}
		assert(memory == m_vertices.size());
	}
	return m_vertices;
}

const std::vector<Edge>& MeshGrid::getEdges()const
{
	if (m_edges.empty())
	{
		const std::vector<cv::Point2i> nexts = { cv::Point2i(1,0),cv::Point2i(0,1) };
		const int memory = DIMENSION_2D * nw * nh + nw + nh;
		m_edges.reserve(memory);
		for (int h = 0; h <= nh; ++h)
		{
			for (int w = 0; w <= nw; ++w)
			{
				const cv::Point2i p1(w, h);
				for (int i = 0; i < nexts.size(); ++i)
				{
					const cv::Point2i p2 = p1 + nexts[i];
					if (p2.x >= 0 && p2.y >= 0 && p2.x <= nw && p2.y <= nh)
					{
						m_edges.emplace_back(p1.x + p1.y * (nw + 1), p2.x + p2.y * (nw + 1));
					}
				}
			}
		}
		assert(memory == m_edges.size());
	}
	return m_edges;
}

const std::vector<Indices>& MeshGrid::getPolygonsIndices()const
{
	if (m_polygons_indices.empty())
	{
		const cv::Point2i nexts[GRID_VERTEX_SIZE] = {
			cv::Point2i(0,0),cv::Point2i(1,0),cv::Point2i(1,1),cv::Point2i(0,1)
		};
		const int memory = nw * nh;
		m_polygons_indices.resize(memory);
		int index = 0;
		for (int h = 0; h < nh; ++h)
		{
			for (int w = 0; w < nw; ++w)
			{
				const cv::Point2i p1(w, h);
				m_polygons_indices[index].indices.reserve(GRID_VERTEX_SIZE);
				for (int n = 0; n < GRID_VERTEX_SIZE; ++n)
				{
					const cv::Point2i p2 = p1 + nexts[n];
					m_polygons_indices[index].indices.emplace_back(p2.x + p2.y * (nw + 1));
				}
				++index;
			}
		}

		assert(memory == m_polygons_indices.size());
	}
	return m_polygons_indices;
}


const std::vector<Indices>& MeshGrid::getPolygonsNeighbors()const
{
	if (m_polygons_neighbors.empty())
	{
		const std::vector<cv::Point2i> nexts = {
			cv::Point2i(1,0),cv::Point2i(0,1),cv::Point2i(-1,0),cv::Point2i(0,-1) };
		const int memory = nw * nh;
		int index = 0;
		m_polygons_neighbors.resize(memory);
		for (int h = 0; h < nh; ++h)
		{
			for (int w = 0; w < nw; ++w)
			{
				const cv::Point2i p1(w, h);
				for (int n = 0; n < nexts.size(); ++n)
				{
					const cv::Point2i p2 = p1 + nexts[n];
					if (p2.x >= 0 && p2.y >= 0 && p2.x < nw && p2.y < nh)
					{
						m_polygons_neighbors[index].indices.emplace_back(p2.x + p2.y * nw);
					}
				}
				++index;
			}
		}
		assert(memory == m_polygons_neighbors.size());
	}
	return m_polygons_neighbors;
}

const std::vector<Indices>& MeshGrid::getPolygonsEdges()const
{
	if (m_polygons_edges.empty())
	{
		const std::vector<int> nexts = { 0,1,3,2 * nw + 1 };
		const int memory = nw * nh;
		m_polygons_edges.resize(memory);
		int index = 0, e_index = 0;
		for (int h = 0; h < nh; ++h)
		{
			for (int w = 0; w < nw; ++w)
			{
				for (int n = 0; n < nexts.size(); ++n)
				{
					m_polygons_edges[index].indices.emplace_back(e_index + nexts[n]);
				}
				m_polygons_edges[index].indices.back() = m_polygons_edges[index].indices.back() - (h == nh - 1) * w;
				++index;
				e_index += 2;
			}
			m_polygons_edges[index - 1].indices[2] = m_polygons_edges[index - 1].indices[2] - 1;
			e_index += 1;
		}
		assert(memory == m_polygons_edges.size());
	}
	return m_polygons_edges;
}


const std::vector<Indices>& MeshGrid::getVertexStructures()const
{
	if (m_vertex_structures.empty())
	{
		const std::vector<cv::Point2i> nexts = {
			cv::Point2i(1,0),cv::Point2i(0,1),cv::Point2i(-1,0),cv::Point2i(0,-1)
		};
		const int memory = (nw + 1) * (nh + 1);
		m_vertex_structures.resize(memory);
		int index = 0;
		for (int h = 0; h <= nh; ++h)
		{
			for (int w = 0; w <= nw; ++w)
			{
				cv::Point2i p1(w, h);
				for (int n = 0; n < nexts.size(); ++n)
				{
					cv::Point2i p2 = p1 + nexts[n];
					if (p2.x >= 0 && p2.y >= 0 && p2.x <= nw && p2.y <= nh)
					{
						m_vertex_structures[index].indices.emplace_back(p2.x + p2.y * (nw + 1));
					}
				}
				++index;
			}
		}
		assert(memory == m_vertex_structures.size());
	}
	return m_vertex_structures;
}


const std::vector<Indices>& MeshGrid::getEdgeStructures()const
{
	if (m_edge_structures.empty())
	{
		const std::vector<cv::Point2i> nexts = { cv::Point2i(1,0),cv::Point2i(0,1) };
		const std::vector<cv::Point2i> grid_neighbor = { cv::Point2i(0,-1),cv::Point2i(-1,0) };
		const int memory = DIMENSION_2D * nh * nw + nh + nw;
		m_edge_structures.resize(memory);
		int index = 0;
		for (int h = 0; h <= nh; ++h)
		{
			for (int w = 0; w <= nw; ++w)
			{
				cv::Point2i p1(w, h);
				for (int n = 0; n < nexts.size(); ++n)
				{
					cv::Point2i p2 = p1 + nexts[n];
					if (p2.x >= 0 && p2.y >= 0 && p2.x <= nw && p2.y <= nh)
					{
						for (int j = 0; j < grid_neighbor.size(); ++j)
						{
							cv::Point2i p3 = p1 + grid_neighbor[n] * j;
							if (p3.x >= 0 && p3.y >= 0 && p3.x < nw && p3.y < nh)
							{
								m_edge_structures[index].indices.emplace_back(p3.x + p3.y * nw);
							}
						}
						++index;
					}
				}
			}
		}

		assert(memory == m_edge_structures.size());
	}
	return m_edge_structures;
}

const std::vector<Indices>& MeshGrid::getTriangulationIndices()const
{
	if (m_triangulation_indices.empty())
	{
		m_triangulation_indices.emplace_back(0, 1, 2);
		m_triangulation_indices.emplace_back(0, 2, 3);
	}
	return m_triangulation_indices;
}

const int MeshGrid::getPolygonVerticesCounts()const
{
	return GRID_VERTEX_SIZE;
}

const std::vector<int>& MeshGrid::getBoundaryVertexIndices()const
{
	if (m_boundary_vertex_indices.empty())
	{
		const int memory = DIMENSION_2D * (nw + nh) + 1;
		m_boundary_vertex_indices.reserve(memory);
		for (int tw = 0; tw < nw; ++tw)
		{
			m_boundary_vertex_indices.emplace_back(tw);
		}
		const int right_bottom = nw * nh + nw + nh;
		for (int rh = nw; rh < right_bottom; rh += (nw + 1))
		{
			m_boundary_vertex_indices.emplace_back(rh);
		}

		const int left_bottm = nh * (nw + 1);
		for (int bw = right_bottom; bw > left_bottm; --bw)
		{
			m_boundary_vertex_indices.emplace_back(bw);
		}

		for (int lh = left_bottm; lh >= 0; lh -= (nw + 1))
		{
			m_boundary_vertex_indices.emplace_back(lh);
		}
		assert(memory == m_boundary_vertex_indices.size());
	}
	return m_boundary_vertex_indices;
}

const std::vector<int>& MeshGrid::getBoundaryEdgeIndices()const
{
	if (m_boundary_edge_indices.empty())
	{
		const int memory = DIMENSION_2D * (nh + nw);
		m_boundary_edge_indices.reserve(memory);
		const int bottom_shift = DIMENSION_2D * nw * nh + nh;
		for (int w = 0; w < nw; ++w)
		{
			m_boundary_edge_indices.emplace_back(2 * w);
			m_boundary_edge_indices.emplace_back(bottom_shift + w);
		}

		const int dh = 2 * nw + 1;
		for (int h = 0; h < nh; ++h)
		{
			int tmp = h * dh;
			m_boundary_edge_indices.emplace_back(tmp + 1);
			m_boundary_edge_indices.emplace_back(tmp + dh - 1);
		}
		assert(memory == m_boundary_edge_indices.size());
	}
	return m_boundary_edge_indices;
}

template<typename T>
InterpolateVertex MeshGrid::getInterpolateVertexTemplate(const cv::Point_<T>& _p)const
{
	const std::vector<Point2>& vertices = getVertices();
	const std::vector<Indices>& grids = getPolygonsIndices();

	const int grid_index = getGridIndexOfPoint(_p);
	const Indices& g = grids[grid_index];

	const std::vector<int> diagonal_indices = { 2,3,0,1 };

	assert(g.indices.size() == GRID_VERTEX_SIZE);

	//每个网格的权重
	std::vector<double> weight(GRID_VERTEX_SIZE);
	double sum_inv = 0;
	for (int i = 0; i < diagonal_indices.size(); ++i)
	{
		Point2 tmp(_p.x - vertices[g.indices[diagonal_indices[i]]].x,
			_p.y - vertices[g.indices[diagonal_indices[i]]].y);
		weight[i] = fabs(tmp.x * tmp.y);
		sum_inv += weight[i];
	}
	sum_inv = 1. / sum_inv;
	for (int i = 0; i < GRID_VERTEX_SIZE; ++i)
	{
		weight[i] = weight[i] * sum_inv;
	}

	return	InterpolateVertex(grid_index, weight);
}

InterpolateVertex MeshGrid::getInterpolateVertex(const cv::Point_<float>& p)const
{
	return getInterpolateVertexTemplate(p);
}

InterpolateVertex MeshGrid::getInterpolateVertex(const cv::Point_<double>& p)const
{
	return getInterpolateVertexTemplate(p);
}

template InterpolateVertex MeshGrid::getInterpolateVertexTemplate<float>(const cv::Point_<float>& p)const;
template InterpolateVertex MeshGrid::getInterpolateVertexTemplate<double>(const cv::Point_<double>& p)const;