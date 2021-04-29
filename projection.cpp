#include "projection.h"
#include <iostream>
#include <set>

const float StereoProjection::fovToFocal(float fov, float d)const
{
    return d / (2.0f * tanf(fov / 2.0f));
}
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
            float focal_length = fovToFocal(fov_rads*M_PI/180.f, R*3);
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


const void MeshOptimization::getImageVerticesBySolving(std::vector<cv::Point2f>& optimized_vertices)
{
    google::InitGoogleLogging("DFWAPCP");

    double* X;
    double* Y;
    const size_t vertices_length = m_vertices.size();
    int nw = w_and_h[0];
    int nh = w_and_h[1];
    int lw = little_mesh_size[0];
    int lh = little_mesh_size[1];
    std::vector<double> S = { 1.,0.,0.,1. };
    X = (double*)malloc(vertices_length*sizeof(double));
    Y = (double*)malloc(vertices_length*sizeof(double));
    for (size_t i = 0; i < vertices_length; i++)
    {
        X[i] =double(m_vertices[i].x);
        Y[i] =double(m_vertices[i].y);
    }

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
        if (m_vertices[i].x == 0 || m_vertices[i].x == nw * lw || m_vertices[i].y == 0 || m_vertices[i].y == nh * lh)
            continue;
        CostFunction* face_energy = new AutoDiffCostFunction<FaceObjectTermEnergy, 2, 1, 1>(
            new FaceObjectTermEnergy(m_stereo_mesh[i].x, m_stereo_mesh[i].y, w, m, S));
        problem.AddResidualBlock(face_energy, nullptr, &X[i], &Y[i]);
    }

    size_t edges_length = m_edges.size();
    for (size_t i = 0; i < edges_length; i++)
    {
        Edge edge = m_edges[i];
        const int ind_e1 = edge.indices[0];
        const int ind_e2 = edge.indices[1];
        cv::Point2f e_vertice = m_vertices[ind_e1] - m_vertices[ind_e2];
        float distance = hypotf(e_vertice.x, e_vertice.y);
        e_vertice.x /= distance;
        e_vertice.y /= distance;
        if (m_vertices[ind_e1].x == 0 || m_vertices[ind_e1].x == nw * lw || m_vertices[ind_e1].y == 0 || m_vertices[ind_e1].y == nh * lh)
            continue;

        CostFunction* line_energy = new AutoDiffCostFunction<LineBlendingTermEnergy, 1, 1, 1,1,1>(
            new LineBlendingTermEnergy(e_vertice.x, e_vertice.y, 0));
        problem.AddResidualBlock(line_energy, nullptr, &X[ind_e1], &Y[ind_e1], &X[ind_e2], &Y[ind_e2]);
        CostFunction* regularization_energy = new AutoDiffCostFunction<RegularizationTermEnergy, 1, 1, 1, 1, 1>(
            new RegularizationTermEnergy());
        problem.AddResidualBlock(regularization_energy, nullptr, &X[ind_e1], &Y[ind_e1], &X[ind_e2], &Y[ind_e2]);
    }

    /*for (size_t i = 0; i < vertices_length; i++)
    {
        std::vector<cv::Point2f> n;
        for (int v = 0; v < m_v_neighbors[i].indices.size(); ++v) {
            int v_index = m_v_neighbors[i].indices[v];
            n.emplace_back(cv::Point2f(float(X[v_index]), float(Y[v_index])));
            if (m_vertices[i].x == 0 || m_vertices[i].x == nw * lw || m_vertices[i].y == 0 || m_vertices[i].y == nh * lh)
                continue;
            cv::Point2f e_vertice = m_vertices[i] - m_vertices[v_index];
            float distance = hypotf(e_vertice.x, e_vertice.y);
            e_vertice.x /= distance;
            e_vertice.y /= distance;
            CostFunction* line_energy = new AutoDiffCostFunction<LineBlendingTermEnergy, 1, 1, 1, 1, 1>(
                new LineBlendingTermEnergy(e_vertice.x, e_vertice.y, 0));
            problem.AddResidualBlock(line_energy, nullptr, &X[i], &Y[i], &X[v_index], &Y[v_index]);
            CostFunction* regularization_energy = new AutoDiffCostFunction<RegularizationTermEnergy, 1, 1, 1, 1, 1>(
                new RegularizationTermEnergy());
            problem.AddResidualBlock(regularization_energy, nullptr, &X[i], &Y[i], &X[v_index], &Y[v_index]);
        }
    }*/

    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = 15;
    Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);


    std::cout << summary.BriefReport() << "\n";
    

    /*for (size_t i = 0; i < vertices_length; i++)
    {
        if (m_vertices[i].x == 0 || m_vertices[i].x == nw * lw || m_vertices[i].y == 0 || m_vertices[i].y == nh * lh)
            continue;
        double x=0., y=0.;
        for (int v = 0; v < m_v_neighbors[i].indices.size(); ++v) {
            int v_index = m_v_neighbors[i].indices[v];
            x += X[v_index];
            y += Y[v_index];
        }
        X[i] = x / 4.0;
        Y[i] = y / 4.0;
    }*/

    for (size_t i = 0; i < vertices_length; i++)
    {
        optimized_vertices.emplace_back(cv::Point2f(float(X[i]), float(Y[i])));
    }


    free(X); free(Y);
}