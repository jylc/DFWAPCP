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

void bilinearIntertpolatioin(cv::Mat& src, cv::Mat& dst, const int rows, const int cols)
{
    //比例尺
    const double scale_row = static_cast<double>(src.rows) / rows;
    const double scale_col = static_cast<double>(src.rows) / cols;

    //扩展src到dst
    //assert(src.channels() == 1 && dst.channels() == 1);

    for (int i = 0; i < rows; ++i)//dst的行
        for (int j = 0; j < cols; ++j)//dst的列
        {
            if(src.at<cv::Vec3b>(i,j)[0]==0)
            for (int k = 0; k < 3; ++k)
            {
                //求插值的四个点
                double y = (i + 0.5) * scale_row + 1;
                double x = (j + 0.5) * scale_col + 1;
                int x1 = static_cast<int>(x);//col对应x
                if (x1 >= (src.cols - 2)) x1 = src.cols - 2;//防止越界
                int x2 = x1 + 1;
                int y1 = static_cast<int>(y);//row对应y
                if (y1 >= (src.rows - 2))  y1 = src.rows - 2;
                int y2 = y1 + 1;

                assert(0 <= x2 && x2 < src.cols && 0 <= y2 && y2 < src.rows);
                //插值公式,参考维基百科矩阵相乘的公式https://zh.wikipedia.org/wiki/%E5%8F%8C%E7%BA%BF%E6%80%A7%E6%8F%92%E5%80%BC

                cv::Matx12d matx = { x2 - x, x - x1 };
                cv::Matx22d matf = { static_cast<double>(src.at<cv::Vec3b>(y1, x1)[k]), static_cast<double>(src.at<cv::Vec3b>(y2, x1)[k]),
                                        static_cast<double>(src.at<cv::Vec3b>(y1, x2)[k]), static_cast<double>(src.at<cv::Vec3b>(y2, x2)[k]) };
                cv::Matx21d maty = {
                    y2 - y,
                    y - y1
                };

                auto  val = (matx * matf * maty);
                dst.at<cv::Vec3b>(i, j)[k] = val(0, 0);
            }
        }
}


template cv::Vec< uchar, 3> getSubpix<uchar, 3>(const cv::Mat& img, const cv::Point2f& pt);
template cv::Point_< float> applyTransform2x3< float>(float x, float y, const cv::Mat& matT);
template cv::Point_<double> applyTransform2x3<double>(double x, double y, const cv::Mat& matT);