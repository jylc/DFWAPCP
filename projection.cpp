#include "projection.h"

cv::Mat StereoProjection::stereo_transformation() {

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
    float extent = u;

    float center_x = n / 2;//ÁÐ
    float center_y = m / 2;//ÐÐ
    cv::Mat res;
    res.create(m / extent, n / extent, CV_8UC3);
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
            float rhoprime = (lambda * r / R) + (1.f - lambda) * (R * (sqrtf(r * r + 1.f) - 1.f)) / (r * (sqrtf(R * R + 1.f) - 1.f));
            //cout << "rhoprime = " << rhoprime << endl;
            float ru = R * (tan(atan2f(r, focal_length) / (2.f - lambda)) / tan(atan2f(R, focal_length) / (2.f - lambda)));
            u = ru * cosf(alpha);
            v = ru * sinf(alpha);
            x = (u + center_x)/extent;
            y = (center_y - v)/extent;
            //
            if (x >= 0 && x < n/extent && y >= 0 && y < m/extent)
            {
                res.at<cv::Vec3b>(y, x) = src_img.at<cv::Vec3b>(i, j);
            }
        }
    }
    return res;
}