#include "projection.h"
#include <iostream>
#include <set>

const float StereoProjection::fovToFocal(float fov, float d)const
{
    return d / (2.0f * tanf(fov / 2.0f));
}
const cv::Mat StereoProjection::stereoTransformation() const{

    float x, y, z, u, v, r;

    float center_x = cols / 2;//列
    float center_y = rows / 2;//行
    cv::Mat res = cv::Mat::zeros(src_img.size(), CV_8UC3);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            x = j;
            y = i;
            u = x - center_x;
            v = center_y - y;

            // perspective
            // apply Z-B transformation to (u,v)
            float lambda = .0f;
            float R = MIN(rows,cols) / 2;
            float focal_length = 600.f;
            float alpha = atan2f(v, u);
            float r = hypotf(u, v);
            float ru = R * (tan(atan2f(r, focal_length) / (2.f - lambda)) / tan(atan2f(R, focal_length) / (2.f - lambda)));
            u = ru * cosf(alpha);
            v = ru * sinf(alpha);
            x = (u + center_x) ;
            y = (center_y - v) ;
            //
            if (x >= 0 && x < cols && y >= 0 && y < rows)
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

    float center_x = cols / 2;//列
    float center_y = rows / 2;//行
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
            float R = MIN(rows, cols) / 2;
            float focal_length = fovToFocal(fov_rads*M_PI/180.f, R*3);
            float alpha = atan2f(v, u);
            float r = hypotf(u, v);
            float ru = R * (tan(atan2f(r, focal_length) / (2.f - lambda)) / tan(atan2f(R, focal_length) / (2.f - lambda)));
            u = ru * cosf(alpha);
            v = ru * sinf(alpha);
            x = u + center_x;
            y = center_y - v;
            if (x >= 0 && x < cols && y >= 0 && y < rows)
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


const void MeshOptimization::getImageVerticesBySolving(std::vector<cv::Point2f>& optimized_vertices)
{
    google::InitGoogleLogging("DFWAPCP");

    double right, bottom;
    const size_t vertices_length = m_vertices.size();
    std::vector<std::pair<double, double>> vertices;
    vertices.reserve(vertices_length);
    int nw = w_and_h[0];
    int nh = w_and_h[1];
    int lw = little_mesh_size[0];
    int lh = little_mesh_size[1];
    std::vector<double> S = { 1.,0.,0.,1. };
    for (size_t i = 0; i < vertices_length; i++)
    {
        vertices.emplace_back(double(m_vertices[i].x), double(m_vertices[i].y));
    }
    right = vertices.back().first;
    bottom = vertices.back().second;
    Problem problem;
    for (size_t i = 0; i < vertices_length; i++)
    {
        float center_x = width / 2, center_y = height / 2;
        float u = m_vertices[i].x - center_x;
        float v = center_y - m_vertices[i].y;
        float alpha = atan2f(v, u);
        float r = hypotf(u, v);
        float rb = 1.f * MAX(width, height), ra = 1.f * log(99) * rb;
        float m = 1.f / (1.f + exp2f(-(r - ra) / rb));

        float w = m_weights[i] ? 1. : 0.;
        CostFunction* face_energy = new AutoDiffCostFunction<FaceObjectTermEnergy, 2, 1, 1>(
            new FaceObjectTermEnergy(m_stereo_mesh[i].x, m_stereo_mesh[i].y, w, m, S));
        problem.AddResidualBlock(face_energy, nullptr, &vertices[i].first, &vertices[i].second);
    }
    size_t index = 0;
    bool x_l = false, y_t = false, x_r = false, y_b = false;
    for (int h = 0; h <= nh; ++h)
    {
        for (int w = 0; w <= nw; ++w)
        {
            if (h==0||h==nh||w==0||w==nw)
            {
                if (h == 0)
                    y_t = true;
                if (h == nh)
                    y_b = true;
                if (w == 0)
                    x_l = true;
                if (w == nw)
                    x_r = true;
                CostFunction* boundary_energy = new AutoDiffCostFunction<MeshBoundaryTermEnergy, 1, 1, 1>(
                    new MeshBoundaryTermEnergy(0, right, 0, bottom, x_l, y_t, x_r, y_b));
                problem.AddResidualBlock(boundary_energy, nullptr, &vertices[index].first, &vertices[index].second);
                x_l = false, y_t = false, x_r = false, y_b = false;
            }
            index++;
        }
    }

    for (size_t i = 0; i < vertices_length; i++)
    {
        for (int v = 0; v < m_v_neighbors[i].indices.size(); ++v) {
            int v_index = m_v_neighbors[i].indices[v];
            cv::Point2f e_vertice = m_vertices[i] - m_vertices[v_index];
            float distance = hypotf(e_vertice.x, e_vertice.y);
            e_vertice.x /= distance;
            e_vertice.y /= distance;
            CostFunction* line_energy = new AutoDiffCostFunction<LineBlendingTermEnergy, 3, 1, 1, 1, 1>(
                new LineBlendingTermEnergy(e_vertice.x, e_vertice.y, 0));
            problem.AddResidualBlock(line_energy, nullptr, &vertices[i].first, &vertices[i].second, &vertices[v_index].first, &vertices[v_index].second);
            CostFunction* regularization_energy = new AutoDiffCostFunction<RegularizationTermEnergy, 2, 1, 1, 1, 1>(new RegularizationTermEnergy());
            problem.AddResidualBlock(regularization_energy, nullptr, &vertices[i].first, &vertices[i].second, &vertices[v_index].first, &vertices[v_index].second);
        }
    }

    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = 15;
    Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);


    std::cout << summary.BriefReport() << "\n";
    

    for (size_t i = 0; i < vertices_length; i++)
    {
        optimized_vertices.emplace_back(cv::Point2f(float(vertices[i].first), float(vertices[i].second)));
    }

    google::ShutdownGoogleLogging();
}