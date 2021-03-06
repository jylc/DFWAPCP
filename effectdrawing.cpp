#include "effectdrawing.h"
#include <cassert>
effectdrawing::effectdrawing(
	std::vector<std::string>& _img_info,
	QWidget *parent)
	:QWidget(parent),
	img_info(_img_info),
	img_dir(_img_info[0]),
	img_name(_img_info[1]),
	mask_dir(_img_info[2]),
	mask_name(_img_info[3]),
	img_data_ptr(img_dir,img_name,mask_dir,mask_name,double(1.0f))
{
	ui.setupUi(this);
}

effectdrawing::~effectdrawing()
{
}
