#pragma once
#ifndef BLENDING_H
#define BLENDING_H
#include <opencv2/opencv.hpp>
#include <vector>
#include "configure.h"

cv::Mat getMatOfLinearBlendWeight(const cv::Mat& image);

std::vector<cv::Mat> getMatsLinearBlendWeight(const std::vector<cv::Mat>& images);
cv::Mat Blending(const std::vector<cv::Mat>& images,
	const std::vector<Point2>& origins,
	const Size2& target_size,
	const std::vector<cv::Mat>& weight_mask,
	const bool ignore_weight_mask = true);


#endif // !BLENDING_H
