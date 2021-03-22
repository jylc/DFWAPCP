#pragma once
#include <opencv2/opencv.hpp>
template<typename T>
T getSubpix(const cv::Mat& img, const cv::Point2f& pt);

template <typename T, size_t n>
cv::Vec<T, n> getSubpix(const cv::Mat& img, const cv::Point2f& pt);

template <typename T>
std::vector<cv::Rect_<T> > getVerticesRects(const std::vector<std::vector<cv::Point_<T> > >& vertices);

template <typename T>
cv::Point_<T> applyTransform2x3(T x, T y, const cv::Mat& matT);