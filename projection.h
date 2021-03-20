#pragma once
#ifndef PROJECTION_H
#define PROJECTION_H
#include "configure.h"
#include <vector>
#include <opencv2/opencv.hpp>
class Projection
{
public:
	Projection(cv::Mat& _src_img)
		:src_img(_src_img) 
	{ 
		height = src_img.rows;
		width = src_img.cols;
	}
	virtual ~Projection() {};

protected:
	cv::Mat& src_img;
	int height;
	int width;
};

class StereoProjection :protected Projection
{
public:
	StereoProjection(cv::Mat& _src_img,float _focal_length=600.f) :Projection(_src_img)
	{
		focal_length = _focal_length;
		n = width;
		m = height;
		scale = 1.f;
		fov_rads = 1.f;
	}


	//Stereographic projection
	cv::Mat stereoTransformation();
	std::vector<Point2> stereoTramsformation(std::vector<Point2>& vertices)const;
private:
	const float calculateExtent()const;
private:
	float focal_length;
	float fov_rads;
	float scale;
	int n;//ÁÐ
	int m;//ÐÐ
};

#endif // !PROJECTION_H
