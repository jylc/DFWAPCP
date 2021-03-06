#pragma once

#include <QWidget>
#include <vector>
#include <string>
#include "ui_effectdrawing.h"
#include "image_data.h"

class effectdrawing : public QWidget
{
	Q_OBJECT

public:
	effectdrawing(std::vector<std::string>& _img_info, QWidget *parent = Q_NULLPTR);
	~effectdrawing();

private:


public slots:

private:
	Ui::effectdrawing ui;

private:
	std::vector<std::string>& img_info;
	std::string& img_dir, & img_name, & mask_dir, & mask_name;
	double focal_length;

	ImageData img_data_ptr;			//ÕºœÒ–≈œ¢
};
