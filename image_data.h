#pragma once
#ifndef IMAGE_DATA
#define IMAGE_DATA
#include "configure.h"
#include <vector>
#include <string>
#include <opencv2/opencv.hpp>
#include  "mesh_grid.h"
class LineData
{
public:
	Point2 m_data[2];
	double m_width, m_length;
	LineData(const Point2& _a, const Point2& _b,
		const double _width, const double  _length);
private:
};


class ImageData
{
public:
	std::unique_ptr<Mesh2D> m_mesh_2d;
	cv::Mat m_img;
	const std::string m_file_dir,m_file_name,m_mask_file_dir,m_mask_file_name;
	ImageData(
		const std::string& _file_dir,
		const std::string& _file_name,
		const std::string& _mask_file_dir,
		const std::string& _mask_file_name);

	ImageData() = default;

	ImageData& operator=(const ImageData& data)
	{
		ImageData tmp(data.m_file_dir, data.m_file_name,
			data.m_mask_file_dir, data.m_mask_file_name);
		return tmp;
	}

	const std::vector<Point2>& getSrcImageVertices()const;
	const cv::Mat& getMaskImg()const;
	const cv::Mat getIntersectedImg(const std::vector<Point2>& vertices,int flag=0)const;//面部检测框与模板相交
	const std::vector<Point2> getStereoImgVertices()const;//stereo projection后的图像点坐标
	const std::vector<Point2> getOptimizedStereoImgVertices()const;
	const std::vector<Edge>& getEdges()const;
	const std::vector<Indices>& getVertexStructures()const;
	const cv::Mat meshTransform(const std::vector<Point2>& new_vertices)const;
	void clear();
	const std::vector<cv::Rect2i> faceDetected()const;
	const std::vector<bool> faceMaskWeight()const;
	const std::vector<int> getCountOfWAndH()const;
	const std::vector<double> getLittleMeshSize()const;
	const void drawVerticesOnImg(cv::Mat& srcImg,const std::vector<Point2>& oldVertices,const std::vector<Point2>& newVertices,const std::vector<bool>& weights)const;

private:
	const std::vector<Point2> mixedMesh(const std::vector<Point2>& old_vertices, const std::vector<Point2>& new_vertices,const std::vector<bool>& weight)const;

private:
	mutable std::vector<LineData> m_img_lines;
	mutable cv::Mat m_stereo_img;
	mutable cv::Mat m_mask_img;
	mutable std::vector<cv::Rect2i> m_face_region;
	int cols, rows;
};

#endif // !IMAGE_DATA
