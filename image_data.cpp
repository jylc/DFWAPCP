#include "image_data.h"
#include <opencv2/imgproc/types_c.h>
#include <opencv2/photo.hpp>
#include "face_detect.h"
#include "projection.h"
#include "blending.h"
#include "transform.h"
LineData::LineData(const Point2& _a, const Point2& _b,
	const double _width, const double  _length):
	m_width(_width),m_length(_length)
{
	m_data[0] = _a;
	m_data[1] = _b;
}

ImageData::ImageData(const std::string& _file_dir, 
	const std::string& _file_name,
	const std::string& _mask_file_dir,
	const std::string& _mask_file_name,
	const double _focal_length):
m_file_dir(_file_dir),
m_file_name(_file_name),
m_mask_file_dir(_mask_file_dir),
m_mask_file_name(_mask_file_name),
m_focal_length(_focal_length)
{
	m_img = cv::imread(m_file_dir + m_file_name);
	if (m_img.empty())
	{
		fprintf(stderr, "cannot load image %s", (m_file_dir + m_file_name).c_str());
		exit(-1);
	}
	
	m_rgba_img = cv::imread(m_file_dir + m_file_name, cv::IMREAD_UNCHANGED);
	if (m_rgba_img.empty())
	{
		fprintf(stderr, "cannot load image %s", (m_file_dir + m_file_name).c_str());
		exit(-1);
	}

	m_mask_img = cv::imread(m_mask_file_dir + m_mask_file_name);
	if (m_mask_img.empty())
	{
		fprintf(stderr, "cannot load image %s", (m_mask_file_dir + m_mask_file_name).c_str());
		exit(-1);
	}

	float original_img_size = m_img.rows * m_img.cols;
	if (original_img_size > DOWN_SAMPLE_IMAGE_SIZE)
	{
		float scale = sqrt(DOWN_SAMPLE_IMAGE_SIZE / original_img_size);	
		cv::resize(m_img, m_img,cv::Size(), scale, scale);
		cv::resize(m_rgba_img, m_rgba_img, cv::Size(), scale, scale);
		cv::resize(m_mask_img, m_mask_img, cv::Size(), scale, scale);
	}

	assert(m_rgba_img.channels() >= 3);
	if (m_rgba_img.channels() == 3)
	{
		cv::cvtColor(m_rgba_img, m_rgba_img, CV_BGR2BGRA);
	}

	std::vector<cv::Mat> channels;
	cv::split(m_rgba_img, channels);
	alpha_mask = channels[3];
	m_mesh_2d = std::make_unique<MeshGrid>(m_img.cols, m_img.rows);
}

ImageData::ImageData(){}

const cv::Mat& ImageData::getGreyImage()const
{
	if(m_grey_img.empty())
	{
		cv::cvtColor(m_img, m_grey_img, CV_BGR2GRAY);
	}
	return m_grey_img;
}

const cv::Mat& ImageData::getSrcImage()const
{
	return m_img;
}

const cv::Mat& ImageData::getAlphaMaskImg()const
{
	return alpha_mask;
}

const cv::Mat& ImageData::getMaskImg()const
{
	if (m_mask_img.empty())
	{
		fprintf(stderr, "mask img is null\n");
		exit(0);
	}
	return m_mask_img;
}

void ImageData::clear()
{
	m_img.release();
	m_rgba_img.release();
	alpha_mask.release();
	m_grey_img.release();
	m_img_lines.clear();
}

const void ImageData::meshInfo()const
{
	std::cout << "nw = " << m_mesh_2d->nw <<
		" , nh = " << m_mesh_2d->nh <<
		" ,lw = " << m_mesh_2d->lw <<
		" ,lh = " << m_mesh_2d->lh <<
		"\n mesh vertices counts = "<<m_mesh_2d->getVertices().size()<<
		"\n mesh vertices center's count = "<<m_mesh_2d->getPolygonsCenter().size() <<
		"\n mesh edges = "<<m_mesh_2d->getEdges().size()<< std::endl;

}

const cv::Mat ImageData::getIntersectedImg(cv::Mat img,int flag)const
{
	cv::Mat result;
	std::vector<short*> rectangle_infos;
	FaceDetect face(m_img);
	rectangle_infos = face.Detect();    //面部框的信息及其中点信息
	result = face.Processed();

	cv::Mat clone_mask = img.clone();
	for (auto it = rectangle_infos.cbegin(); it != rectangle_infos.cend(); ++it)
	{
		int confidence = (*it)[0];
		int x = (*it)[1];
		int y = (*it)[2];
		int w = (*it)[3];
		int h = (*it)[4];
		if (confidence <= 70)
			continue;

		char s_score[255];
		snprintf(s_score, 255, "%d", confidence);
		cv::putText(clone_mask, s_score, cv::Point(x, y - 3), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 255, 0), 1);
		rectangle(clone_mask, Rect(x - w / 2, y - h, 2 * w, 2 * h), Scalar(0, 255, 0), 2);
	}
	StereoProjection stereo_projection(clone_mask);
	if (flag==0)
	{//原图像网格点
		std::vector<Point2> old_vertices = m_mesh_2d->getVertices();
		std::vector<Point2> new_vertices = old_vertices;
		std::vector<bool> weights = faceMaskWeight();
		drawVerticesOnImg(clone_mask, old_vertices, new_vertices, weights);
	}
	else if(flag==1)
	{//混合点
		std::vector<Point2> old_vertices = m_mesh_2d->getVertices();
		std::vector<Point2> new_vertices = stereo_projection.stereoTramsformation(old_vertices);
		std::vector<bool> weights = faceMaskWeight();
		drawVerticesOnImg(clone_mask, old_vertices, new_vertices, weights);
		//drawVerticesOnImg(clone_mask, new_vertices, new_vertices, weights);
	}
	else if(flag==2)
	{
		std::vector<Point2> old_vertices = m_mesh_2d->getVertices();
		std::vector<Point2> new_vertices = stereo_projection.stereoTramsformation(old_vertices);
		std::vector<bool> weights = faceMaskWeight();
		drawVerticesOnImg(clone_mask, new_vertices, new_vertices, weights);
	}
	
	return clone_mask;
}

const void ImageData::drawLinesOnImg(cv::Mat& srcImg, std::vector<Point2>& oldVertices, std::vector<Point2>& newVertices, std::vector<bool>& weights)const
{
}

const cv::Mat ImageData::getStereoImg()const
{
	//TODO:focal length暂时设置为600.f
	cv::Mat src_img = m_img.clone();
	StereoProjection stereo_projection(src_img);

	m_stereo_img = stereo_projection.stereoTransformation();
	cv::Mat stereo_img = m_stereo_img.clone();
	//绘制点，观察效果
	/*std::vector<Point2> old_vertices = m_mesh_2d->getVertices();
	std::vector<Point2> new_vertices = stereo_projection.stereoTramsformation(old_vertices);
	drawVerticesOnImg(stereo_img, new_vertices, new_vertices, std::vector<bool>(old_vertices.size(), false));*/
	return stereo_img;
}

//stereo image的顶点
const std::vector<Point2> ImageData::getVerticesOfStereoImg()const
{
	cv::Mat src_img = m_img.clone();
	StereoProjection stereo_projection(src_img);
	std::vector<Point2> old_vertices = m_mesh_2d->getVertices();
	std::vector<Point2> new_vertices = stereo_projection.stereoTramsformation(old_vertices);
	return new_vertices;
}

//检测面部区域
const std::vector<cv::Rect2i> ImageData::faceDetected()const
{
	cv::Mat result;
	std::vector<short*> rectangle_infos;
	FaceDetect face(m_img);
	rectangle_infos = face.Detect();    //面部框的信息及其中点信息
	result = face.Processed();

	for (auto it = rectangle_infos.cbegin(); it != rectangle_infos.cend(); ++it)
	{
		//int confidence = int(it[0]);
		int confidence = (*it)[0];
		int x = (*it)[1];
		int y = (*it)[2];
		int w = (*it)[3];
		int h = (*it)[4];

		//char s_score[255];
		//snprintf(s_score, 255, "%d", confidence);
		m_face_region.emplace_back(x - w / 2, y - h, 2 * w, 2 * h);
	}
	return m_face_region;
}

//网格变换
const cv::Mat ImageData::meshTransform()const
{
	std::vector<Point2> old_vertices = m_mesh_2d->getVertices();
	std::vector<Indices> polygons_indices = m_mesh_2d->getPolygonsIndices();
	std::vector<Indices> triangulation_indices = m_mesh_2d->getTriangulationIndices();
	std::vector<cv::Mat> affine_transforms;
	cv::Mat src_img = m_img.clone();
	StereoProjection stereo_projection(src_img);

	std::vector<Point2> new_vertices = stereo_projection.stereoTramsformation(old_vertices);
	affine_transforms.reserve(polygons_indices.size() * (triangulation_indices.size()));
	const Point2 shift(0.5, 0.5);
	const int NO_GRID = -1, TRIANGLE_COUNT = 3, PRECISION = 0;
	cv::Mat polygon_index_mask(src_img.rows + shift.y, src_img.cols + shift.x, CV_32SC1, Scalar::all(NO_GRID));
	int label = 0;

	std::vector<bool> face_mask_weight = faceMaskWeight();
	//TODO:只对面部区域形变！！！
	
	std::vector<Point2> tmp_vertices;
	tmp_vertices.reserve(old_vertices.size());
	for (size_t i = 0; i < old_vertices.size(); ++i)
	{
		if (face_mask_weight[i])
			tmp_vertices.emplace_back(new_vertices[i]);
		else
			tmp_vertices.emplace_back(old_vertices[i]);
	}

	for (int i = 0; i < polygons_indices.size(); ++i)//四边形
	{
		for (int j = 0; j < triangulation_indices.size(); ++j)//三角形顶点
		{
			Indices index = triangulation_indices[j];//每个三角形
			int index_0 = polygons_indices[i].indices[index.indices[0]];
			int index_1 = polygons_indices[i].indices[index.indices[1]];
			int index_2 = polygons_indices[i].indices[index.indices[2]];
			const Point2i contour[] =
			{
				old_vertices[index_0],
				old_vertices[index_1],
				old_vertices[index_2],
			};
			fillConvexPoly(polygon_index_mask, contour, TRIANGLE_COUNT, label, LINE_AA, PRECISION);

			Point2f src[] = {
				old_vertices[index_0],
				old_vertices[index_1],
				old_vertices[index_2]
			};

			Point2f dst[] =
			{
				tmp_vertices[index_0],
				tmp_vertices[index_1],
				tmp_vertices[index_2],
			};
			
			affine_transforms.emplace_back(getAffineTransform(src, dst));
			++label;
		}
	}
	
	float scale = 1;
	cv::Mat image = cv::Mat::zeros((src_img.rows + shift.y)*scale, (src_img.cols + shift.x)*scale, CV_8UC3);

	for (int y = 0; y < image.rows; ++y)
	{
		for (int x = 0; x < image.cols; ++x)
		{
			int polygon_index = polygon_index_mask.at<int>(y, x);
			//在网格范围内
			if (polygon_index != NO_GRID)
			{
				Point2f p_f = applyTransform2x3<FLOAT_TYPE>(x, y, affine_transforms[polygon_index]);
				p_f.x *= scale;
				p_f.y *= scale;
				if (y == 1 && x == 1)
				{
					fprintf(stderr,"test: x1 = %d, y1 = %d, x2 = %f, y2 = %f\n", x, y, p_f.x, p_f.y);
				}
				if (p_f.x >= 0 && p_f.y >= 0 &&
					p_f.x < src_img.cols &&
					p_f.y < src_img.rows) {
					cv::Vec3b c = src_img.at<cv::Vec3b>(y,x);
					image.at<Vec3b>(p_f.y, p_f.x) = cv::Vec3b(c[0], c[1], c[2]);
				}
			}
		}
	}
	cv::Mat clone_image = image.clone();

	//双线性插值
	//bilinearIntertpolatioin(image, clone_image, image.rows, image.cols);
	return clone_image;
}

const std::vector<bool> ImageData::faceMaskWeight()const
{
	const std::vector<Point2> vertices = m_mesh_2d->getVertices();
	const std::vector<cv::Rect2i> face_region = faceDetected();
	std::vector<bool> weight(vertices.size(),false);
	cv::Mat src_img = getSrcImage();
	cv::Mat mask_img = getMaskImg();
	cv::Mat img = cv::Mat::zeros(src_img.size(), CV_8UC3);
	int index = 0;
	for (auto it = vertices.begin(); it != vertices.end(); ++it,++index)
	{
		Point2 p = *it;
		int x = int(p.x);
		int y = int(p.y);
		for (auto k = face_region.begin(); k != face_region.end(); ++k)
		{
			cv::Rect2i rect = *k;
			int a = rect.x, b = rect.y, c = rect.x + rect.width, d = rect.y + rect.height;
			if (x >= a && y >= b && x <= c && y <= d&&!weight[index])
			{
				if (x >= 0 && x < mask_img.cols && 
					y >= 0 && y < mask_img.rows &&
					mask_img.at<cv::Vec3b>(p)[0]!=0)
					weight[index] = true;
			}
		}
	}
	std::cout << "face detected region size " << face_region.size() << std::endl;
	return weight;
}

const std::vector<int> ImageData::getCountOfWAndH()
{
	if (m_mesh_2d == nullptr)
	{
		fprintf(stderr, "[ImageData getCountOfWAndH]: m_mesh_2d is null\n");
		exit(-1);
	}
	std::vector<int> list = { m_mesh_2d->nw,m_mesh_2d->nh };
	return list;
}

const void ImageData::drawVerticesOnImg(cv::Mat& srcImg, std::vector<Point2>& oldVertices, std::vector<Point2>& newVertices, std::vector<bool>& weights)const
{
	for (size_t it = 0; it < oldVertices.size(); ++it)
	{
		if (weights[it])
			cv::circle(srcImg, newVertices[it], 1, cv::Scalar(0, 255, 0));
		else
			cv::circle(srcImg, oldVertices[it], 1, cv::Scalar(0, 0, 255));
	}
}

//网格形变获得形变的顶点之后使用仿射变换应用在每个点上？
const cv::Mat ImageData::blendImages()const
{
	std::vector<cv::Mat> weight_mask, new_weight_mask;
	std::vector<Point2> origins;
	std::vector<cv::Mat> images = getImages();
	weight_mask = getMatsLinearBlendWeight(images);
	std::vector<std::vector<Point2>> all_vertices(images.size());
	for (auto i = 0; i != images.size(); ++i)
	{
		if (0 == i)
			all_vertices.push_back(m_mesh_2d->getVertices());
		else
			all_vertices.push_back(getVerticesOfStereoImg());
	}

	//std::vector<Rect_<float>> rects = getVerticesRects<float>(all_vertices);
	return cv::Mat();
}

const std::vector<cv::Mat> ImageData::getImages()const
{
	getStereoImg();
	if (m_img.empty()&&m_stereo_img.empty())
	{
		fprintf(stderr, "[ImageData_getImages] Cannot load source image or stereo image\n");
		exit(-1);
	}
	if (m_images.empty())
	{
		m_images.reserve(2);
		m_images.emplace_back(m_img);
		m_images.emplace_back(m_stereo_img);
	}
	return m_images;
}