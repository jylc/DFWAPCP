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
	//�沿ʶ�𣬷����沿�����Ϣ
	std::vector<short*> Detect();

	//���ش�����ͼƬ
	Mat Processed();

private:
	std::string m_path;//ͼƬ·��
	Mat m_image;//��ȡͼƬ
	Mat m_resultImage;//������ͼƬ
	int* m_pResults;//�����
	unsigned char* m_pBuffer;
	std::vector<short*> m_detectResults;//���м����
};

#endif