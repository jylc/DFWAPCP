#include "image_data.h"
#include <opencv2/imgproc/types_c.h>
#include "face_detect.h"
#include "projection.h"
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

const cv::Mat ImageData::getIntersectedImg()const
{
	cv::Mat result;
	std::vector<short*> rectangle_infos;
	FaceDetect face(m_img);
	rectangle_infos = face.Detect();    //�沿�����Ϣ�����е���Ϣ
	result = face.Processed();

	cv::Mat clone_mask = getMaskImg();
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
	return clone_mask;
}

cv::Mat ImageData::getStereoImg()
{
	//TODO:focal length��ʱ����Ϊ600.f
	StereoProjection stereo_projection(m_img);
	return stereo_projection.stereo_transformation();
}

//����沿����
const std::vector<cv::Rect2i> ImageData::faceDetected()
{
	cv::Mat result;
	std::vector<short*> rectangle_infos;
	FaceDetect face(m_img);
	rectangle_infos = face.Detect();    //�沿�����Ϣ�����е���Ϣ
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

const std::vector<cv::Point2i> ImageData::meshTransform()const
{
	const std::vector<Point2> vertices = m_mesh_2d->getVertices();
	const std::vector<Indices> polygons_indices = m_mesh_2d->getPolygonsIndices();
	std::vector<cv::Point2i> tram_vertices;
	
	for (int k = 0; k < vertices.size(); ++k)
	{
		//TODO: ����任 
		//tram_vertices.push_back(StereoProjection<int>::cylinder(m_img.cols, m_img.rows,
			//cv::Point2i(round(vertices[k].x), round(vertices[k].y))));
	}

	cv::Mat test_img = cv::Mat::zeros(m_img.size(), CV_8UC3);
	std::vector<cv::Point2i>::iterator it = tram_vertices.begin();
	/*for (; it != tram_vertices.end(); ++it)
	{
		cv::circle(test_img, *it, 1, cv::Scalar(255, 0, 0), 1);
	}*/
	return tram_vertices;
}

const std::vector<bool> ImageData::faceMaskWeight()
{
	const std::vector<Point2> vertices = m_mesh_2d->getVertices();
	const std::vector<cv::Rect2i> face_region = faceDetected();
	std::vector<bool> weight(vertices.size(),false);
	cv::Mat src_img = getSrcImage();
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



void ImageData::setStereoImgFromGL(cv::Mat stereo_img)
{
	m_stereo_img = stereo_img;
}

cv::Mat ImageData::getStereoImgFromGL()
{
	return m_stereo_img;
}