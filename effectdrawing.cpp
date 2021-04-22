#include "effectdrawing.h"
#include <cassert>

effectdrawing::effectdrawing(
	std::shared_ptr<ImageData>& _img_data_ptr,
	QWidget *parent)
	:QWidget(parent),
	img_data_ptr(_img_data_ptr)
{
	ui.setupUi(this);
	cv::Mat src_img = img_data_ptr->getIntersectedImg(img_data_ptr->getSrcImage(),-1);
	cv::Mat middle_img = img_data_ptr->getIntersectedImg(img_data_ptr->getOptimizedStereoImg(), 3);
	cv::Mat stereo_img = img_data_ptr->getIntersectedImg(img_data_ptr->getStereoImg(), 3);
	if (src_img.empty() || stereo_img.empty()||middle_img.empty())
	{
		fprintf(stderr, "[Effectdrawing] cannot get source image or effected image");
	}
	else
	{
		//封装为函数执行会导致界面不能更新
		cv::Mat tmp1, tmp2, tmp3;
		cv::cvtColor(src_img, tmp1, cv::COLOR_BGR2RGB);
		QImage img1 = QImage((const unsigned char*)tmp1.data, tmp1.cols, tmp1.rows, tmp1.cols * tmp1.channels(), QImage::Format_RGB888);

		cv::cvtColor(stereo_img, tmp2, cv::COLOR_BGR2RGB);
		QImage img2 = QImage((const unsigned char*)tmp2.data, tmp2.cols, tmp2.rows, tmp2.cols * tmp2.channels(), QImage::Format_RGB888);

		cv::cvtColor(middle_img, tmp3, cv::COLOR_BGR2RGB);
		QImage img3 = QImage((const unsigned char*)tmp3.data, tmp3.cols, tmp3.rows, tmp3.cols * tmp3.channels(), QImage::Format_RGB888);
		

		ui.srcimg->setPixmap(QPixmap::fromImage(img1));
		ui.srcimg->adjustSize();

		ui.middleimg->setPixmap(QPixmap::fromImage(img3));
		ui.middleimg->adjustSize();

		ui.stereoimg->setPixmap(QPixmap::fromImage(img2));
		ui.stereoimg->adjustSize();
	}
}

effectdrawing::~effectdrawing()
{
}
