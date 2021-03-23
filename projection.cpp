#include "projection.h"

const cv::Mat StereoProjection::stereoTransformation() const{

    float x, y, z, u, v, r;

    float center_x = n / 2;//列
    float center_y = m / 2;//行
    cv::Mat res = cv::Mat::zeros(src_img.size(), CV_8UC3);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            x = j;
            y = i;
            u = x - center_x;
            v = center_y - y;

            // perspective
            // apply Z-B transformation to (u,v)
            float lambda = .0f;
            float R = MIN(m,n) / 2;
            float focal_length = 600.f;
            float alpha = atan2f(v, u);
            float r = hypotf(u, v);
            float ru = R * (tan(atan2f(r, focal_length) / (2.f - lambda)) / tan(atan2f(R, focal_length) / (2.f - lambda)));
            u = ru * cosf(alpha);
            v = ru * sinf(alpha);
            x = (u + center_x) ;
            y = (center_y - v) ;
            //
            if (x >= 0 && x < n && y >= 0 && y < m)
            {
                res.at<cv::Vec3b>(y, x) = src_img.at<cv::Vec3b>(i, j);
            }
        }
    }

    cv::Mat inpaint_mask = cv::Mat::zeros(res.size(), CV_8U);
    for (size_t y = 0; y < res.rows; y++)
    {
        for (size_t x = 0; x < res.cols; x++)
        {
            int a = res.at<cv::Vec3b>(y, x)[0];
            int b = res.at<cv::Vec3b>(y, x)[1];
            int c = res.at<cv::Vec3b>(y, x)[2];

            if (a == 0 && b == 0 && c == 0)
            {
                inpaint_mask.at<uchar>(y, x) = 255;
            }
        }
    }
    inpaint(res, inpaint_mask, res, 3, cv::INPAINT_TELEA);
    return res;
}

const std::vector<Point2> StereoProjection::stereoTramsformation(const std::vector<Point2>& vertices)const
{
    float phi, x, y, z, u, v, r;

    float center_x = n / 2;//列
    float center_y = m / 2;//行
    std::vector<Point2> newVertices;
    newVertices.reserve(vertices.size());
    
    for (auto& it : vertices)
    {
            x = it.x;
            y = it.y;
            //转为笛卡尔坐标系
            u = x - center_x;
            v = center_y - y;

            // perspective
            // apply Z-B transformation to (u,v)
            float lambda = .0f;
            float R = MIN(m, n) / 2;
            float focal_length = 600.f;
            float alpha = atan2f(v, u);
            float r = hypotf(u, v);
            float ru = R * (tan(atan2f(r, focal_length) / (2.f - lambda)) / tan(atan2f(R, focal_length) / (2.f - lambda)));
            u = ru * cosf(alpha);
            v = ru * sinf(alpha);
            x = u + center_x;
            y = center_y - v;
            if (x >= 0 && x < n && y >= 0 && y < m)
            {
                newVertices.push_back(Point2(x, y));
            }
            else
            {
                newVertices.push_back(Point2(it.x, it.y));
            }
    }
    return newVertices;
}

const float StereoProjection::calculateExtent()const
{
    float lambda, phi, x, y, z, u, v, r, theta;

    //calculating the extent of the projection for the given FOV
    lambda = fov_rads;
    phi = 0.f;
    y = sinf(phi);
    x = -sinf(lambda) * cosf(phi);
    z = -cosf(lambda) * cosf(phi);
    u = 2.f * x / (1.f - z);
    v = 2.f * y / (1.f - z);
    r = hypotf(u, v);
    theta = atan2f(u, v);
    r *= scale;
    u = -r * sinf(theta);
    v = r * cosf(theta);
    x = (4.f * u) / (u * u + v * v + 4.f);
    y = (4.f * v) / (u * u + v * v + 4.f);
    z = (u * u + v * v - 4.f) / (u * u + v * v + 4.f);
    u = x / (-z);
    v = y / (-z);
    return u;
}
