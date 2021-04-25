#include "projection.h"
#include <iostream>

#include <set>
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

const void MeshOptimization::getFaceObjectiveTerm(std::vector<Triplet<double>>& triplets,
    std::vector<std::pair<int, double>>& b_vector,bool flag)const
{
    const size_t face_region_size = m_face_region.size();
    const size_t vertices_size = m_vertices.size();
    std::vector<float> a(face_region_size, 1.f), b(face_region_size, 0.f);
    std::vector<std::pair<float, float>> t(face_region_size, std::pair<float, float>(0.f, 0.f));
    cv::Mat matrix = cv::Mat::zeros(DIMENSION_2D, DIMENSION_2D, CV_32FC1);
    cv::Mat matrix_t;
    
    float ws = 2000.f, st = 1.f;
    //默认使用参数优化
    if (flag)
    {
        for (size_t v_index = 0; v_index < m_vertices.size(); v_index++)
        {
            bool flag = false;
            float stereo_x = m_stereo_mesh[v_index].x;
            float stereo_y = m_stereo_mesh[v_index].y;
            float origin_x = m_vertices[v_index].x;
            float origin_y = m_vertices[v_index].y;
            
            int face_region_index = 0;
            for (auto k = m_face_region.begin(); k != m_face_region.end(); ++k,++face_region_index)
            {
                cv::Rect2i rect = *k;
                int a = rect.x, b = rect.y, c = rect.x + rect.width, d = rect.y + rect.height;
                if (origin_x >= a && origin_y >= b && origin_x <= c && origin_y <= d)
                {
                    flag = true;
                    break;
                }
            }

            if (flag && m_weights[v_index])
            {
                //TODO: wi,mi
                //triplets相当于A

                //求径向距离r
                float center_x = width / 2, center_y = height / 2;
                float u = m_vertices[v_index].x - center_x;
                float v = center_y - m_vertices[v_index].y;
                float alpha = atan2f(v, u);
                float r = hypotf(u, v);

                float rb = 1.f * MAX(width, height), ra = 1.f * log(99) * rb;
                float m = 1.f / (1.f + exp2f(-(r - ra) / rb));
                //float m = 1.f;
                matrix.at<float>(0, 0) = sqrtf(m);
                matrix.at<float>(0, 1) = 0.f;
                matrix.at<float>(1, 0) = 0.f;
                matrix.at<float>(1, 1) = sqrtf(m);
                cv::transpose(matrix, matrix_t);
                cv::Mat res = matrix_t.inv(cv::DECOMP_SVD) * (matrix_t * matrix);
                triplets.emplace_back(DIMENSION_2D * v_index, DIMENSION_2D * v_index, res.at<float>(0, 0));
                triplets.emplace_back(DIMENSION_2D * v_index + 1, DIMENSION_2D * v_index + 1, res.at<float>(1, 1));

                //TODO: 待添加lamda(Sk) Sk[ak bk -bk ak],tk,ws=2000,st=1
                //b_vector相当于B
                float b_x = 0.f, b_y = 0.f;
                b_x = sqrtf(m) * (stereo_x * a[face_region_index] + stereo_y * b[face_region_index] + t[face_region_index].first);
                b_y = sqrtf(m) * (stereo_y * a[face_region_index] - stereo_x * b[face_region_index] + t[face_region_index].second);

                b_vector.emplace_back(DIMENSION_2D * v_index, b_x);
                b_vector.emplace_back(DIMENSION_2D * v_index + 1, b_y);
            }
            else
            {
                triplets.emplace_back(DIMENSION_2D * v_index, DIMENSION_2D * v_index, 1);
                triplets.emplace_back(DIMENSION_2D * v_index + 1, DIMENSION_2D * v_index + 1, 1);
                b_vector.emplace_back(DIMENSION_2D * v_index, origin_x);
                b_vector.emplace_back(DIMENSION_2D * v_index + 1, origin_y);
            }
            float tmp = flag ? sqrtf(ws / (1.f * vertices_size)) * (a[face_region_index] - st) : 0;
            b_vector.emplace_back(DIMENSION_2D * (v_index + vertices_size), 1000);
            b_vector.emplace_back(DIMENSION_2D * (v_index + vertices_size) + 1, 1000);
        }
    }
    else
    {
        for (size_t v_index = 0; v_index < m_vertices.size(); v_index++)
        {//遍历所有点
            bool flag = false;
            float stereo_x = m_stereo_mesh[v_index].x;
            float stereo_y = m_stereo_mesh[v_index].y;
            float origin_x = m_vertices[v_index].x;
            float origin_y = m_vertices[v_index].y;
            //size_t index;
            for (auto k = m_face_region.begin(); k != m_face_region.end(); ++k)
            {
                cv::Rect2i rect = *k;
                int a = rect.x, b = rect.y, c = rect.x + rect.width, d = rect.y + rect.height;
                if (origin_x >= a && origin_y >= b && origin_x <= c && origin_y <= d)
                {
                    flag = true;
                }
            }

            //假设Ax=b:x为[v0......vn]^-1
            if (flag&&m_weights[v_index])
            {//面部矩形框内且在面部区域
                //triplets(i,j)处存放(A^T*A)
                triplets.emplace_back(DIMENSION_2D * v_index, DIMENSION_2D * v_index, 1);
                triplets.emplace_back(DIMENSION_2D * v_index + 1, DIMENSION_2D * v_index + 1, 1);
                b_vector.emplace_back(DIMENSION_2D * v_index, stereo_x);
                b_vector.emplace_back(DIMENSION_2D * v_index + 1,  stereo_y);
            }
            else
            {//面部矩形框外或者在面部矩形框内但不在面部区域
                triplets.emplace_back(DIMENSION_2D * v_index, DIMENSION_2D * v_index, 1);
                triplets.emplace_back(DIMENSION_2D * v_index + 1, DIMENSION_2D * v_index + 1, 1);
                b_vector.emplace_back(DIMENSION_2D * v_index, origin_x);
                b_vector.emplace_back(DIMENSION_2D * v_index + 1, origin_y);
            }
        }
    }
}

const void MeshOptimization::getLineBlendingTerm(std::vector<Triplet<double>>& triplets,
    std::vector<std::pair<int, double>>& b_vector)const
{
    cv::Mat e = cv::Mat::zeros(DIMENSION_2D, 1, CV_32FC1);
    for (size_t v_index=0;v_index<m_edges.size();++v_index)
    {
        b_vector.emplace_back(v_index, 0);
    }
}

//有问题，不能直接将面部点单独抽离出来
//TODO: 待修改
const void MeshOptimization::verticeSetsOfFaces()const
{
    
}

std::vector<Point2> MeshOptimization::getImageVerticesBySolving(std::vector<Triplet<double>>& triplets,
    const std::vector<std::pair<int, double>>& b_vector)
{
    LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
    SparseMatrix<double> A((m_vertices.size() + m_edges.size()) * DIMENSION_2D, m_vertices.size() * DIMENSION_2D);
    VectorXd b = VectorXd::Zero((m_vertices.size() + m_edges.size()) * DIMENSION_2D), x;

    A.setZero();
    A.setFromTriplets(triplets.begin(), triplets.end());
    for (int i = 0; i < b_vector.size(); ++i) {
        b(b_vector[i].first, 0) = b_vector[i].second;
    }

    lscg.compute(A);
    x = lscg.solve(b);
    std::vector<Point2> vertices;
    vertices.reserve(m_vertices.size());
    for (size_t i = 0; i < m_vertices.size(); i++)
    {
        vertices.emplace_back(x[i * DIMENSION_2D], x[i * DIMENSION_2D + 1]);
    }
    return vertices;
}



const void MeshOptimization::getImageVerticesBySolving(std::vector<cv::Point2f>& optimized_vertices)
{
    google::InitGoogleLogging("DFWAPCP");

    double* X;
    double* Y;
    const size_t vertices_length = m_vertices.size();

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
        std::vector<cv::Point2f> n;
        n.emplace_back(cv::Point2f(float(X[ind_e2]), float(Y[ind_e2])));
        CostFunction* line_energy = new AutoDiffCostFunction<LineBlendingTermEnergy, 3, 1, 1,1,1>(
            new LineBlendingTermEnergy(e_vertice.x, e_vertice.y, 0));
        problem.AddResidualBlock(line_energy, nullptr, &X[ind_e1], &Y[ind_e1], &X[ind_e2], &Y[ind_e2]);

        CostFunction* regularization_energy = new AutoDiffCostFunction<RegularizationTermEnergy, 2, 1, 1>(
            new RegularizationTermEnergy(n));
        problem.AddResidualBlock(regularization_energy, nullptr, &X[ind_e1], &Y[ind_e1]);
        n.clear();
    }

    /*for (size_t i = 0; i < vertices_length; i++)
    {
        std::vector<cv::Point2f> n;
        for (int v = 0; v < m_v_neighbors[i].indices.size(); ++v) {
            int v_index = m_v_neighbors[i].indices[v];
            n.emplace_back(cv::Point2f(float(X[v_index]), float(Y[v_index])));
        }
        CostFunction* regularization_energy = new AutoDiffCostFunction<RegularizationTermEnergy, 2, 1,1>(
            new RegularizationTermEnergy(n));
        problem.AddResidualBlock(regularization_energy, nullptr, &X[i], &Y[i]);
    }*/

    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);


    std::cout << summary.BriefReport() << "\n";
    for (size_t i = 0; i < vertices_length; i++)
    {
        optimized_vertices.emplace_back(cv::Point2f(float(X[i]), float(Y[i])));
    }

    free(X); free(Y);
}