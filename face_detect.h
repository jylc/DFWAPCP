#pragma once
#ifndef FACE_DETECT_H
#define FACE_DETECT_H
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <vector>
#include <string>
//define the buffer size. Do not change the size!
#define DETECT_BUFFER_SIZE 0x20000
using namespace cv;

class FaceDetect
{
public:
	FaceDetect(Mat image);
	~FaceDetect();
	//面部识别，返回面部框的信息
	std::vector<short*> Detect();

	//返回处理后的图片
	Mat Processed();

private:
	std::string m_path;//图片路径
	Mat m_image;//读取图片
	Mat m_resultImage;//处理后的图片
	int* m_pResults;//检测结果
	unsigned char* m_pBuffer;
	std::vector<short*> m_detectResults;//所有检测结果
};

#endif