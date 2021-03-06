#pragma once
#ifndef CONFIGURE_H
#define CONFIGURE_H
#include <opencv2/core/core.hpp>
#include <opencv2/core/types.hpp>

#define DETECT_BUFFER_SIZE 0x20000

const int GRID_SIZE = 40;
const int GRID_VERTEX_SIZE = 4;
const int DOWN_SAMPLE_IMAGE_SIZE = 800 * 600;

typedef float FLOAT_TYPE;
typedef cv::Size_<FLOAT_TYPE> Size2;
typedef cv::Point_<FLOAT_TYPE> Point2;
typedef cv::Rect_<FLOAT_TYPE> Rect2;

const int DIMENSION_2D = 2;
#endif // !CONFIGURE_H
