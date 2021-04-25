#pragma once
#ifndef PROJECTION_H
#define PROJECTION_H
#include "configure.h"
#include "mesh_2d.h"
#include <vector>
#include "ceres/ceres.h"
#include "glog/logging.h"
using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::CauchyLoss;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;

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
		const std::vector<Indices>& _v_neighbors,
		float _focal_length = 600.0f)
		:Projection(_src_img),
		m_face_region(_face_region),
		m_stereo_mesh(_stereo_mesh),
		m_vertices(_vertices),
		m_edges(_edges),
		m_weights(_weights),
		m_v_neighbors(_v_neighbors)
	{}
	//默认使用各种参数进行网格优化
	/*===================START===================*/
	//解析解求解，未完成
	const void getFaceObjectiveTerm(std::vector<Triplet<double>>& triplets,
		std::vector<std::pair<int, double>>& b_vector, bool flag = true)const;
	const void getLineBlendingTerm(std::vector<Triplet<double>>& triplets,
		std::vector<std::pair<int, double>>& b_vector)const;
	std::vector<Point2> getImageVerticesBySolving(std::vector<Triplet<double>>& triplets,
		const std::vector<std::pair<int, double>>& b_vector);
	/*==================END====================*/

	/*===================START===================*/
	//数值解求解
	const void getImageVerticesBySolving(std::vector<cv::Point2f>& optimized_vertices);
	/*==================END====================*/
private:
	const void verticeSetsOfFaces()const;
	
private:
	const std::vector<cv::Rect2i>& m_face_region;//面部区域
	const std::vector<Point2>& m_stereo_mesh;//形变后网格点
	const std::vector<Point2>& m_vertices;//原网格点
	const std::vector<Edge>& m_edges;
	const std::vector<bool>& m_weights;//权重
	const std::vector<Indices>& m_v_neighbors;
};

struct FaceObjectTermEnergy
{
	FaceObjectTermEnergy(float _s_x, float _s_y, float _weight, float _m, std::vector<double> _S) :
		weight(_weight), m(_m), s_x(_s_x), s_y(_s_y), S(_S)
	{
	};

	template <typename T>
	bool operator()(const T* const x, const T* const y, T* residual) const {

		residual[0] = (T(weight) * T(m)) * (x[0] - S[0] * T(s_x)) * (x[0] - S[0] * T(s_x));
		residual[1] = (T(weight) * T(m)) * (y[0] - S[3] * T(s_y)) * (y[0] - S[3] * T(s_y));
		residual[0] *= T(4);
		residual[1] *= T(4);
		return true;
	}
private:
	float weight;
	float m;
	float s_x;
	float s_y;
	std::vector<double> S;
};

struct LineBlendingTermEnergy
{
	LineBlendingTermEnergy(float _e_x,float _e_y,float _e_z)
	:e_x(_e_x),e_y(_e_y),e_z(_e_z){}

	template<typename T>
	bool operator()(const T* const x1, const T* const y1, const T* const x2, const T* const y2, T* residual)const
	{
		T x = (x1[0] - x2[0]), y = (y1[0] - y2[0]), z = T(0);
		residual[0] = (y * T(e_z) - z * T(e_y))* (y * T(e_z) - z * T(e_y));
		residual[1] = (z * T(e_x) - x * T(e_z))* (z * T(e_x) - x * T(e_z));
		residual[2] = (x * T(e_y) - y * T(e_x))* (x * T(e_y) - y * T(e_x));
		residual[0] *= T(4);
		residual[1] *= T(4);
		residual[2] *= T(4);
		return true;
	}

private:
	float e_x;
	float e_y;
	float e_z;
};


struct RegularizationTermEnergy
{
	RegularizationTermEnergy(std::vector<cv::Point2f> _neighbors):neighbors(_neighbors){}

	template<typename T>
	bool operator()(const T* const x1, const T* const y1, T* residual)const
	{
		T a = T(0.), b = T(0.);
		for (size_t i = 0; i < neighbors.size(); i++)
		{
			a += (x1[0] - T(neighbors[i].x))* (x1[0] - T(neighbors[i].x));
			b += (y1[0] - T(neighbors[i].y))* (y1[0] - T(neighbors[i].y));
		}
		residual[0] = a * T(0.5);
		residual[1] = b * T(0.5);
		return true;
	}
private:
	std::vector<cv::Point2f> neighbors;
};
#endif // !PROJECTION_H
