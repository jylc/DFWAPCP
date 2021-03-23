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
	cv::Mat m_img, m_rgba_img, alpha_mask, m_grey_img;
	const std::string m_file_dir,m_file_name,m_mask_file_dir,m_mask_file_name;
	ImageData(
		const std::string& _file_dir,
		const std::string& _file_name,
		const std::string& _mask_file_dir,
		const std::string& _mask_file_name,
		const double _focal_length);

	ImageData();

	ImageData& operator=(const ImageData& data)
	{
		ImageData tmp(data.m_file_dir, data.m_file_name,
			data.m_mask_file_dir, data.m_mask_file_name,
			data.m_focal_length);
		return tmp;
	}

	const cv::Mat& getGreyImage()const;
	const cv::Mat& getSrcImage()const;
	const cv::Mat& getMaskImg()const;
	const cv::Mat& getAlphaMaskImg()const;
	const cv::Mat getIntersectedImg(cv::Mat img, int flag=0)const;//面部检测框与模板相交
	const cv::Mat getStereoImg()const;//stereo projection后的图像
	const std::vector<LineData>& getLines()const;
	const cv::Mat meshTransform()const;
	const void meshInfo()const;
	void clear();
	const std::vector<cv::Rect2i> faceDetected()const;
	const std::vector<bool> faceMaskWeight()const;
	const std::vector<int> getCountOfWAndH();
	const void drawVerticesOnImg(cv::Mat& srcImg,std::vector<Point2>& oldVertices, std::vector<Point2>& newVertices, std::vector<bool>& weights)const;
	const void drawLinesOnImg(cv::Mat& srcImg, std::vector<Point2>& oldVertices, std::vector<Point2>& newVertices, std::vector<bool>& weights)const;
	const cv::Mat blendImages()const;
	const std::vector<cv::Mat> getImages()const;
	//TODO: 该方法放在此处不太适合，待修改
	const std::vector<Point2> getVerticesOfStereoImg()const;


private:
	mutable std::vector<LineData> m_img_lines;
	mutable double m_focal_length;
	mutable cv::Mat m_stereo_img;
	mutable cv::Mat m_mask_img;
	mutable std::vector<cv::Rect2i> m_face_region;
	mutable std::vector<cv::Mat> m_images;
};

#endif // !IMAGE_DATA
