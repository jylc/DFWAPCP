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


const cv::Mat& ImageData::getGreyImage()const
{
	if(m_grey_img.empty())
	{
		cv::cvtColor(m_img, m_grey_img, CV_BGR2GRAY);
	}
	return m_grey_img;
}

const std::vector<Point2>& ImageData::getSrcImage()const
{
	return m_mesh_2d->getVertices();
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

const cv::Mat ImageData::getIntersectedImg(const std::vector<Point2>& vertices, int flag)const
{
	std::vector<short*> rectangle_infos;
	FaceDetect face(m_img);
	rectangle_infos = face.Detect();    //面部框的信息及其中点信息
	
	cv::Mat result_img = meshTransform(vertices);//点化图
	//cv::Mat result_img = imgTransform(vertices);
	/*for (auto it = rectangle_infos.cbegin(); it != rectangle_infos.cend(); ++it)
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
		cv::putText(result_img, s_score, cv::Point(x, y - 3), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 255, 0), 1);
		rectangle(result_img, Rect(x - w / 2, y - h, 2 * w, 2 * h), Scalar(0, 255, 0), 2);
	}*/
	StereoProjection stereo_projection(m_img);
	const std::vector<Point2>& old_vertices = m_mesh_2d->getVertices();
	const std::vector<Point2>& new_vertices = stereo_projection.stereoTramsformation(old_vertices);
	const std::vector<bool>& face_weights = faceMaskWeight();
	//默认显示原图网格
	//if (flag == -1)
	//	return result_img;
	//else if (flag == 0)
	//	drawVerticesOnImg(result_img, old_vertices, old_vertices, face_weights);
	//else if(flag ==1)//混合
	//	drawVerticesOnImg(result_img, old_vertices, new_vertices, face_weights);
	//else if(flag==2)//stereo
	//	drawVerticesOnImg(result_img, new_vertices, new_vertices, face_weights);
	//else if (flag == 3)//optimized
	//drawVerticesOnImg(result_img, new_vertices, new_vertices, face_weights);
	return result_img;
}

const std::vector<Point2> ImageData::getStereoImg()const
{
	//TODO:focal length暂时设置为600.f
	cv::Mat src_img = m_img.clone();
	StereoProjection stereo_projection(src_img);
	return stereo_projection.stereoTramsformation(m_mesh_2d->getVertices());
}

const std::vector<Point2> ImageData::getOptimizedStereoImg()const
{
	const cv::Mat& src_img = m_img.clone();
	const std::vector<Point2>& old_vertices = m_mesh_2d->getVertices();
	StereoProjection stereo_projection(src_img);
	const std::vector<Point2>& new_vertices = stereo_projection.stereoTramsformation(old_vertices);
	const std::vector<bool>& face_weights = faceMaskWeight();
	const std::vector<cv::Rect2i>& face_region = faceDetected();
	const std::vector<Edge>& edges = m_mesh_2d->getEdges();
	const std::vector<Indices>& v_neighbors = m_mesh_2d->getVertexStructures();
	MeshOptimization mesh_Optimization(src_img, face_region, new_vertices, old_vertices, edges, face_weights, v_neighbors);
	//std::vector<Triplet<double>> triplets;
	//std::vector<std::pair<int, double>> b_vector;
	//triplets.reserve((old_vertices.size()+ edges.size()) * DIMENSION_2D * old_vertices.size() * DIMENSION_2D);
	//b_vector.reserve((old_vertices.size()+edges.size()) * DIMENSION_2D);
	//mesh_Optimization.getFaceObjectiveTerm(triplets, b_vector, false);
	//mesh_Optimization.getLineBlendingTerm(triplets, b_vector);
	//mesh_Optimization.getImageVerticesBySolving(triplets, b_vector);
	std::vector<cv::Point2f> optimized_vertices;
	optimized_vertices.reserve(old_vertices.size());
	mesh_Optimization.getImageVerticesBySolving(optimized_vertices);
	return optimized_vertices;
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

		m_face_region.emplace_back(x - w / 2, y - h, 2 * w , 2 * h + h / 4);
	}
	return m_face_region;
}

//网格变换
const cv::Mat ImageData::meshTransform(const std::vector<Point2>& new_vertices)const
{
	const std::vector<Point2> &old_vertices = m_mesh_2d->getVertices();
	const std::vector<Indices> &polygons_indices = m_mesh_2d->getPolygonsIndices();
	const std::vector<Indices> &triangulation_indices = m_mesh_2d->getTriangulationIndices();
	std::vector<cv::Mat> affine_transforms;
	cv::Mat src_img = m_img.clone();
	StereoProjection stereo_projection(src_img);

	affine_transforms.reserve(polygons_indices.size() * (triangulation_indices.size()));
	const Point2 shift(0.5, 0.5);
	const int NO_GRID = -1, TRIANGLE_COUNT = 3, PRECISION = 0;
	cv::Mat polygon_index_mask(src_img.rows + shift.y, src_img.cols + shift.x, CV_32SC1, Scalar::all(NO_GRID));
	int label = 0;
	std::vector<Point2> tmp_vertices = new_vertices;
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
				tmp_vertices[index_2]
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
	

	cv::Mat inpaint_mask = cv::Mat::zeros(image.size(), CV_8U);
	for (size_t y = 0; y < image.rows; y++)
	{
		for (size_t x = 0; x < image.cols; x++)
		{
			int a = image.at<cv::Vec3b>(y, x)[0];
			int b = image.at<cv::Vec3b>(y, x)[1];
			int c = image.at<cv::Vec3b>(y, x)[2];

			if (a==0&&b==0&&c==0)
			{
				inpaint_mask.at<uchar>(y, x) = 255;
			}
		}
	}
	inpaint(image, inpaint_mask, image, 3, INPAINT_TELEA);
	return image;
}

const std::vector<bool> ImageData::faceMaskWeight()const
{
	const std::vector<Point2> vertices = m_mesh_2d->getVertices();
	const std::vector<cv::Rect2i> face_region = faceDetected();
	std::vector<bool> weight(vertices.size(),false);
	cv::Mat src_img = m_img;
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

const std::vector<int> ImageData::getCountOfWAndH()const
{
	if (m_mesh_2d == nullptr)
	{
		fprintf(stderr, "[ImageData getCountOfWAndH]: m_mesh_2d is null\n");
		exit(-1);
	}
	std::vector<int> list = { m_mesh_2d->nw,m_mesh_2d->nh };
	return list;
}

const void ImageData::drawVerticesOnImg(cv::Mat& src_img,const std::vector<Point2>& old_vertices, const std::vector<Point2>& new_vertices,const std::vector<bool>& weights)const
{
	auto& vertices = mixedMesh(old_vertices, new_vertices, weights);
	int i = 0;
	for (auto& vertice : vertices)
	{
		cv::circle(src_img, vertice, 1, weights[i++] ? cv::Scalar(0, 255, 0) : cv::Scalar(0, 0, 255));
	}
}


const std::vector<Point2> ImageData::mixedMesh(
	const std::vector<Point2>& old_vertices, 
	const std::vector<Point2>& new_vertices, 
	const std::vector<bool>& weight)const
{
	std::vector<Point2> mixed_mesh_vertices;
	mixed_mesh_vertices.reserve(old_vertices.size());

	for (size_t i = 0; i < old_vertices.size(); i++)
	{
		mixed_mesh_vertices.emplace_back(weight[i]? new_vertices[i] : old_vertices[i]);
	}
	return mixed_mesh_vertices;
}

