#include "test_img.h"
#include <stdio.h>
#include <string>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "image_data.h"
int TestImg::test_img()
{
    std::string file_dir("G:/VSProject/DFWAPCP/data/pics");
    std::string file_name("paper_2.JPG");
    std::string mask_dir("G:/VSProject/DFWAPCP/data/masks");
    std::string mask_name("paper_2.png");

    std::string imgSrcPath(file_dir + file_name);
    std::string imgMaskPath(mask_dir + mask_name);
    cv::Mat imageSrc = cv::imread(imgSrcPath);  //Ô­Ê¼Í¼Ïñ
    cv::Mat imageMask = cv::imread(imgMaskPath);//Í¼ÏñÑÚÄ¤

    if (imageSrc.empty())
    {
        fprintf(stderr, "cannot load image %s", imgSrcPath.c_str());
        return -1;
    }

    if (imageMask.empty())
    {
        fprintf(stderr, "cannot load image %s", imgSrcPath.c_str());
        return -1;
    }

    float focal_length = 0.f;
    ImageData img_data(file_dir, file_name, mask_dir, mask_name, focal_length);
    const std::vector<Point2>& mesh_vertices = img_data.m_mesh_2d->getVertices();
    const  std::vector<Edge> mesh_edges = img_data.m_mesh_2d->getEdges();

    img_data.meshInfo();
    cv::Mat src_img = img_data.getSrcImage();
    cv::Mat intersected_img = img_data.getIntersectedImg();
    cv::Mat stereo_img = img_data.getStereoImg();

    imshow("intersected img", intersected_img);
    cv::waitKey();
    img_data.clear();
    return 0;
}
