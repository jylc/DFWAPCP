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
	StereoProjection(const cv::Mat& _src_img,float _fov_rads = 97.f) :Projection(_src_img)
	{
		cols = width;
		rows = height;
		scale = 1.f;
		fov_rads = _fov_rads;
	}

	//Stereographic projection
	const cv::Mat stereoTransformation()const;
	const std::vector<Point2> stereoTramsformation(const std::vector<Point2>& vertices)const;
private:
	const float calculateExtent()const;
	const float fovToFocal(float fov,float d)const;
private:
	float fov_rads;
	float scale;
	int cols;//列
	int rows;//行
};

class MeshOptimization:protected Projection
{
public:
	MeshOptimization(const cv::Mat& _src_img,
		const std::vector<Point2>& _stereo_mesh,
		const std::vector<Point2>& _vertices,
		const std::vector<Edge>& _edges,
		const std::vector<bool>& _weights,
		const std::vector<Indices>& _v_neighbors,
		const std::vector<int>& _w_and_h,
		const std::vector<double> _little_mesh_size)
		:Projection(_src_img),
		m_stereo_mesh(_stereo_mesh),
		m_vertices(_vertices),
		m_edges(_edges),
		m_weights(_weights),
		m_v_neighbors(_v_neighbors),
		little_mesh_size(_little_mesh_size),
		w_and_h(_w_and_h)
	{}

	/*===================START===================*/
	//数值解求解
	const void getImageVerticesBySolving(std::vector<cv::Point2f>& optimized_vertices);
	/*==================END====================*/

	~MeshOptimization() {
		std::cout << "release mesh optimization" << std::endl;
	}
private:
	const std::vector<Point2>& m_stereo_mesh;//形变后网格点
	const std::vector<Point2>& m_vertices;//原网格点
	const std::vector<Edge>& m_edges;
	const std::vector<bool>& m_weights;//权重
	const std::vector<Indices>& m_v_neighbors;
	const std::vector<int> w_and_h;
	const std::vector<double> little_mesh_size;
};

struct FaceObjectTermEnergy
{
	FaceObjectTermEnergy(float _s_x, float _s_y, float _weight, float _m, std::vector<double> _S) :
		weight(_weight), m(_m), s_x(_s_x), s_y(_s_y), S(_S)
	{
	};

	template <typename T>
	bool operator()(const T* const x, const T* const y, T* residual) const {
		T a = T(0.0), b = T(0.0);
		a = (T(weight) * T(m)) * (x[0] - S[0] * T(s_x));
		b = (T(weight) * T(m)) * (y[0] - S[3] * T(s_y));
		residual[0] = a * T(4);
		residual[1] = b * T(4);
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
		//vi=(x1,y1,z) , vj=(x2,y2,z) ,eij is the unit vector along the direction pi-pj
		//(vi-vj) x eij
		T x = (x1[0] - x2[0]), y = (y1[0] - y2[0]), z = T(0), a = T(0), b = T(0), c = T(0);
		a = (y * T(e_z) - z * T(e_y));
		b = (z * T(e_x) - x * T(e_z));
		c = (x * T(e_y) - y * T(e_x));
		residual[0] = a * T(1);
		residual[1] = b * T(1);
		residual[2] = c * T(1);
		return true;
	}

private:
	float e_x;
	float e_y;
	float e_z;
};


struct RegularizationTermEnergy
{
	template<typename T>
	bool operator()(const T* const x1, const T* const y1, const T* const x2, const T* const y2, T* residual)const
	{
		//vi=(x1,y1) ,vj=(x2,y2)
		//vi-vj
		T a = T(0.0), b = T(0.0);
		a += (x1[0] - x2[0]);
		b += (y1[0] - y2[0]);
		residual[0] = a * T(.2);
		residual[1] = b * T(.2);
		return true;
	}
private:
};

struct MeshBoundaryTermEnergy
{
	MeshBoundaryTermEnergy(double _left,double _right,double _top,double _bottom,bool _xp1,bool _yp1, bool _xp2, bool _yp2):
	left(_left),right(_right),top(_top),bottom(_bottom),x_l(_xp1),y_t(_yp1), x_r(_xp2), y_b(_yp2) {}
	template<typename T>
	bool operator()(const T* const x1, const T* const y1, T* residual)const
	{
		T l = T(0.0), r = T(0.0), t = T(0.0), b = T(0.0);

		if (x1[0] > left && x_l == true)
			l = x1[0] * x1[0];
		if (x1[0] < right && x_r == true)
			r = (x1[0] - T(right)) * (x1[0] - T(right));

		if (y1[0] > top && y_t == true)
			t = y1[0] * y1[0];
		if (y1[0] < bottom && y_b == true)
			b = (y1[0] - T(bottom)) * (y1[0] - T(bottom));
		residual[0] = T(2.0) * (l + r + t + b);
		return true;
	}
private:
	double left, right, top, bottom;
	bool x_l, y_t, x_r, y_b;
};

#endif // !PROJECTION_H
