#include "transform.h"
#include "configure.h"
#include <cmath>
template<typename T>
T getSubpix(const cv::Mat& img, const cv::Point2f& pt)
{
	cv::Mat patch;
	cv::getRectSubPix(img, cv::Size(1, 1), pt, patch);
	return patch.at<T>(0, 0);
}

template <typename T, size_t n>
cv::Vec<T, n> getSubpix(const cv::Mat& img, const cv::Point2f& pt) {
    cv::Mat patch;
    cv::getRectSubPix(img, cv::Size(1, 1), pt, patch);
    return patch.at<cv::Vec<T, n> >(0, 0);
}

template <typename T>
std::vector<cv::Rect_<T> > getVerticesRects(const std::vector<std::vector<cv::Point_<T> > >& vertices)
{
    vector<cv::Rect_<T> > result;
    result.reserve(vertices.size());
    for (int i = 0; i < vertices.size(); ++i) {
        T min_ix = MAXFLOAT, max_ix = -MAXFLOAT;
        T min_iy = MAXFLOAT, max_iy = -MAXFLOAT;
        for (int j = 0; j < vertices[i].size(); ++j) {
            min_ix = std::min(min_ix, vertices[i][j].x);
            max_ix = std::max(max_ix, vertices[i][j].x);
            min_iy = std::min(min_iy, vertices[i][j].y);
            max_iy = std::max(max_iy, vertices[i][j].y);
        }
        result.emplace_back(min_ix, min_iy,
            max_ix - min_ix, max_iy - min_iy);
    }
    return result;
}

template <typename T>
cv::Point_<T> applyTransform2x3(T x, T y, const cv::Mat& matT) {
    return cv::Point_<T>((matT.at<double>(0, 0) * x + matT.at<double>(0, 1) * y + matT.at<double>(0, 2)),
        (matT.at<double>(1, 0) * x + matT.at<double>(1, 1) * y + matT.at<double>(1, 2)));
}

template cv::Vec< uchar, 3> getSubpix<uchar, 3>(const cv::Mat& img, const cv::Point2f& pt);

template cv::Point_< float> applyTransform2x3< float>(float x, float y, const cv::Mat& matT);
template cv::Point_<double> applyTransform2x3<double>(double x, double y, const cv::Mat& matT);