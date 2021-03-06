#pragma once
#ifndef PROJECTION_FUNCTION
#define PROJECTION_FUNCTION
#include <opencv2/opencv.hpp>

#define  PI 3.14159

class ProjectionFunction
{

};

//球面投影
template<typename T>
class StereoProjection :public ProjectionFunction
{
public:
	//计算结果
	 static cv::Point_<T> Calculate(int width, int height, double focal_length, cv::Point_<T>& center_point, cv::Point_<T>& p1)
	{
		double d = std::min(width, height);
		double r0 = d / (2.0 * tan(0.5 * atan(d / (2.0 * focal_length))));
		double rp = (p1.x - center_point.x) * (p1.x - center_point.x) + (p1.y - center_point.y) * (p1.y - center_point.y);
		rp = std::sqrt(rp);
		double ru = r0 * tan(0.5 * atan(rp / focal_length));
		//double ru = Calculate1(width, height, focal_length, rp);
		T x, y;
		x = rp / ru * (p1.x - center_point.x) + center_point.x;
		y = rp / ru * (p1.y - center_point.y) + center_point.y;
		return cv::Point_<T>(x, y);
	}

	 //图
	 static cv::Mat cylinder(cv::Mat& src)
	 {
		 cv::Mat img_result = src.clone();
		 for (int i = 0; i < img_result.rows; i++)
		 {
			 for (int j = 0; j < img_result.cols; j++)
			 {
				 img_result.at<cv::Vec3b>(i, j) = 0;
			 }
		 }
		 int W = src.cols;
		 int H = src.rows;
		 float r = W / (2 * tan(PI / 6));
		 float k = 0;
		 float fx = 0;
		 float fy = 0;
		 for (int i = 0; i < src.rows; i++)
		 {
			 for (int j = 0; j < src.cols; j++)
			 {
				 k = sqrt((float)(r * r + (W / 2 - j) * (W / 2 - j)));
				 fx = r * sin(PI / 6) + r * sin(atan((j - W / 2) / r));
				 fy = H / 2 + r * (i - H / 2) / k;
				 int ix = (int)fx;
				 int iy = (int)fy;
				 if (ix < W && ix >= 0 && iy < H && iy >= 0)
				 {
					 img_result.at<cv::Vec3b>(iy, ix) = src.at<cv::Vec3b>(i, j);
				 }
			 }
		 }
		 return img_result;
	 }

	 //点
	 static cv::Point2i cylinder(int W,int H,cv::Point2i& _p)
	 {
		 float r = W / (2 * tan(PI / 6));
		 float k = 0;
		 float fx = 0;
		 float fy = 0;
		 int j = _p.x;
		 int i = _p.y;

		 k = sqrt((float)(r * r + (W / 2 - j) * (W / 2 - j)));
		 fx = r * sin(PI / 6) + r * sin(atan((j - W / 2) / r));
		 fy = H / 2 + r * (i - H / 2) / k;
		 return cv::Point2i(int(fx), int(fy));
	 }
private:
	StereoProjection()=delete;
	~StereoProjection() {}

private:
};
#endif // !PROJECTION_FUNCTION

