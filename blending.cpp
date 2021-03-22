#include "blending.h"
#include <cmath>

//线性融合
cv::Mat getMatOfLinearBlendWeight(const cv::Mat& image) {
	cv::Mat result(image.size(), CV_32FC1, cv::Scalar::all(0));
	for (int y = 0; y < result.rows; ++y)
	{
		int w_y = std::min(y + 1, result.rows - y);
		for (int x = 0; x < result.cols; ++x)
		{
			result.at<float>(y, x) = std::min(x + 1, result.cols - x) * w_y;
		}
	}
	return result;
}

std::vector<cv::Mat> getMatsLinearBlendWeight(const std::vector<cv::Mat>& images) {
	std::vector<cv::Mat> result;
	result.reserve(images.size());
	for (int i = 0; i < images.size(); ++i)
	{
		result.push_back(getMatOfLinearBlendWeight(images[i]));
	}
	return result;
}

//images：图像集（现在固定2张）
//origins：每张图的原点（左上角）
//target_size：结果图像的大小
//weight_masks：每张图的权重
//ignore_weight_mask：是否使用权重
cv::Mat Blending(const std::vector<cv::Mat>& images,
	const std::vector<Point2>& origins,
	const Size2& target_size,
	const std::vector<cv::Mat>& weight_masks,
	const bool ignore_weight_mask)
{
	if (images.size() != 2)
	{
		fprintf(stderr, "[Blending] Should have 2 images\n");
		exit(-1);
	}
	cv::Mat result = cv::Mat::zeros(round(target_size.height), round(target_size.width),CV_8UC4);

	std::vector<Rect2> rects;
	rects.reserve(images.size());
	for (int i = 0; i < origins.size(); ++i)
		rects.emplace_back(origins[i], images[i].size());

	for (int y = 0; y < result.rows; ++y)
	{
		for (int x = 0; x < result.cols; ++x)
		{
			cv::Point2i p(x, y);
			cv::Vec3f pixel_sum(0, 0, 0);
			float weight_sum = 0.f;
			for (int i = 0; i < rects.size(); ++i)
			{
				cv::Point2i pv(round(x - origins[i].x), round(y - origins[i].y));
				if (pv.x >= 0 && pv.x <= images[i].cols &&
					pv.y >= 0 && pv.y <= images[i].rows)
				{
					cv::Vec4b v = images[i].at<cv::Vec4b>(pv);
					cv::Vec3f value = cv::Vec3f(v[0], v[1], v[2]);
					if (ignore_weight_mask)
					{
						if (v[3] > 127)
						{
							pixel_sum += value;
							weight_sum += 1.f;
						}
					}
					else
					{
						float weight = weight_masks[i].at<float>(pv);
						pixel_sum += value;
						weight_sum += weight;
					}
				}
			}
			if (weight_sum)
			{
				pixel_sum /= weight_sum;
				result.at<cv::Vec4b>(p) = cv::Vec4b(round(pixel_sum[0]), round(pixel_sum[1]), round(pixel_sum[2]), 255);
			}
		}
	}
	return result;
}