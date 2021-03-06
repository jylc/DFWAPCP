#include "blending.h"
#include <cmath>

cv::Mat Blending(const std::vector<cv::Mat>& images,
	const std::vector<Point2>& origins,
	const Size2& target_size,
	const std::vector<cv::Mat>& weight_masks,
	const bool ignore_weight_mask)
{
	cv::Mat result = cv::Mat::zeros(std::round(target_size.width), std::round(target_size.height), CV_8UC4);
	cv::Mat image,weight_mask;
	if (1 == images.size() && 1 == weight_masks.size())
	{
		image = images.front();
		weight_mask = weight_masks.front();
	}
	else
	{
		fprintf(stderr, "there have too much images!\n");
		exit(-1);
	}
	for (int y = 0; y < result.rows; ++y)//行
	{
		for (int x = 0; x < result.cols; ++x)//列
		{
			cv::Point2i p(x, y);
			cv::Vec3b v = image.at<cv::Vec3b>(p);
			////TODO: 待添加权重

			result.at<cv::Vec4b>(p) = cv::Vec4b(std::round(v[0]), std::round(v[1]), std::round(v[2]),255);
		}
	}
	return result;
}