#pragma once
#ifndef PROJECTION_H
#define PROJECTION_H
#include "configure.h"
#include "mesh_2d.h"
#include <vector>

class Projection
{
public:
	Projection(const cv::Mat& _src_img)
		:src_img(_src_img) 
	{ 
		height = src_img.rows;
		width = src_img.cols;
	}
	virtual ~Projection() {};

protected:
	const cv::Mat& src_img;
	int height;
	int width;
};

class StereoProjection :protected Projection
{
public:
	StereoProjection(const cv::Mat& _src_img,float _focal_length=600.f) :Projection(_src_img)
	{
		focal_length = _focal_length;
		n = width;
		m = height;
		scale = 1.f;
		fov_rads = 1.f;
	}

	//Stereographic projection
	const cv::Mat stereoTransformation()const;
	const std::vector<Point2> stereoTramsformation(const std::vector<Point2>& vertices)const;
private:
	const float calculateExtent()const;
private:
	float focal_length;
	float fov_rads;
	float scale;
	int n;//列
	int m;//行
};

class MeshOptimization:protected Projection
{
public:
	MeshOptimization(const cv::Mat& _src_img,
		const std::vector<cv::Rect2i>& _face_region,
		const std::vector<Point2>& _stereo_mesh,
		const std::vector<Point2>& _vertices,
		const std::vector<Edge>& _edges,
		const std::vector<bool>& _weights,
		float _focal_length = 600.0f)
		:Projection(_src_img),
		m_face_region(_face_region),
		m_stereo_mesh(_stereo_mesh),
		m_vertices(_vertices),
		m_edges(_edges),
		m_weights(_weights)
	{}
	//默认使用各种参数进行网格优化
	const void getFaceObjectiveTerm(std::vector<Triplet<double>>& triplets,
		std::vector<std::pair<int, double>>& b_vector, bool flag = true)const;
	const void getLineBlendingTerm(std::vector<Triplet<double>>& triplets,
		std::vector<std::pair<int, double>>& b_vector)const;
	std::vector<Point2> getImageVerticesBySolving(std::vector<Triplet<double>>& triplets,
		const std::vector<std::pair<int, double>>& b_vector);

private:
	const void verticeSetsOfFaces()const;
private:
	const std::vector<cv::Rect2i>& m_face_region;//面部区域
	const std::vector<Point2>& m_stereo_mesh;//形变后所有点
	const std::vector<Point2>& m_vertices;//所有点
	const std::vector<Edge>& m_edges;
	const std::vector<bool>& m_weights;//权重
};

#endif // !PROJECTION_H
