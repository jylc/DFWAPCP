#include "mesh_2d.h"
Mesh2D::Mesh2D(const int _cols, const int _rows)
{
	nw = _cols / GRID_SIZE + (_cols % GRID_SIZE != 0);//宽度的网格数
	nh = _rows / GRID_SIZE + (_rows % GRID_SIZE != 0);//高度的网格数
	lw = _cols / (double)nw;//网格尺寸
	lh = _rows / (double)nh;//网格尺寸

}

Mesh2D::~Mesh2D() {}

const std::vector<Point2>& Mesh2D::getPolygonsCenter()const
{
	if (m_polygon_center.empty())
	{
		const std::vector<Point2>& vertices = getVertices();
		const std::vector<Indices>& polygons_indices = getPolygonsIndices();
		m_polygon_center.reserve(polygons_indices.size());
		for (int i = 0; i < polygons_indices.size(); ++i)
		{
			Point2 center(0, 0);
			for (int j = 0; j < polygons_indices[i].indices.size(); ++j)
			{
				center += vertices[polygons_indices[i].indices[j]];
			}
			m_polygon_center.emplace_back(center);
		}
	}
	return m_polygon_center;
}

template<typename T>
int Mesh2D::getGridIndexOfPoint(const cv::Point_<T>& p)const
{
	cv::Point2i grid_p(p.x / lw, p.y / lh);
	grid_p.x = grid_p.x - (grid_p.x == nw);
	grid_p.y = grid_p.y - (grid_p.y == nh);
	return grid_p.x + grid_p.y * nw;
}

template int Mesh2D::getGridIndexOfPoint<float>(const cv::Point_<float>& p)const;
template int Mesh2D::getGridIndexOfPoint<double>(const cv::Point_<double>& p)const;


