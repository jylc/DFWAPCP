#pragma once

#include <QWidget>
#include <vector>
#include <string>
#include "ui_effectdrawing.h"
#include "image_data.h"
#include <qimage.h>

class effectdrawing : public QWidget
{
	Q_OBJECT

public:
	effectdrawing(std::shared_ptr<ImageData>& _img_data_ptr, QWidget* parent = Q_NULLPTR);
	~effectdrawing();

private:


public slots:

private:
	Ui::effectdrawing ui;

private:
	//ImageData& img_data;			//ÕºœÒ–≈œ¢
	std::shared_ptr<ImageData>& img_data_ptr;
};
