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
	cv::Mat m_img, m_rgba_img, alpha_mask, m_grey_img, m_mask_img;
	std::vector<cv::Rect2i> m_face_region;
	
	const std::string m_file_dir,m_file_name,m_mask_file_dir,m_mask_file_name;
	ImageData::ImageData(
		const std::string& _file_dir,
		const std::string& _file_name,
		const std::string& _mask_file_dir,
		const std::string& _mask_file_name,
		const double _focal_length);
	
	const cv::Mat& getGreyImage()const;
	const cv::Mat& getSrcImage()const;
	const cv::Mat& getMaskImg()const;
	const cv::Mat& getAlphaMaskImg()const;
	const cv::Mat getIntersectedImg()const;//面部检测框与模板相交
	cv::Mat getStereoImg();//stereo projection后的图像
	const std::vector<LineData>& getLines()const;
	const std::vector<cv::Point2i> meshTransform()const;
	const void meshInfo()const;
	void clear();
	const std::vector<cv::Rect2i> faceDetected();
	const std::vector<bool> faceMaskWeight();

private:
private:
	mutable std::vector<LineData> m_img_lines;
	mutable double m_focal_length;
};

#endif // !IMAGE_DATA
