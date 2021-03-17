#include "effectdrawing.h"
#include <cassert>

effectdrawing::effectdrawing(
	std::shared_ptr<ImageData>& _img_data_ptr,
	QWidget *parent)
	:QWidget(parent),
	img_data_ptr(_img_data_ptr)
{
	ui.setupUi(this);
	//
	cv::Mat src_img = img_data_ptr->getSrcImage();
	cv::Mat effected_img = img_data_ptr->getStereoImg();
	if (src_img.empty() || effected_img.empty())
	{
		fprintf(stderr, "[Effectdrawing] cannot get source image or effected image");
	}
	else
	{
		//封装为函数执行会导致界面不能更新
		cv::Mat tmp1, tmp2;
		cv::cvtColor(src_img, tmp1, cv::COLOR_BGR2RGB);
		QImage img1 = QImage((const unsigned char*)tmp1.data, tmp1.cols, tmp1.rows, tmp1.cols * tmp1.channels(), QImage::Format_RGB888);

		cv::cvtColor(effected_img, tmp2, cv::COLOR_BGR2RGB);
		QImage img2 = QImage((const unsigned char*)tmp2.data, tmp2.cols, tmp2.rows, tmp2.cols * tmp2.channels(), QImage::Format_RGB888);

		

		ui.srcimg->setPixmap(QPixmap::fromImage(img1));
		ui.srcimg->adjustSize();

		ui.effectedimg->setPixmap(QPixmap::fromImage(img2));
		ui.effectedimg->adjustSize();

		
	}
}

effectdrawing::~effectdrawing()
{
}
