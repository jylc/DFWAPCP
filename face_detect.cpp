#include "face_detect.h"
#include <iostream>
#include <stdio.h>
#include <facedetection/facedetectcnn.h>
#include <facedetection/facedetection_export.h>
FaceDetect::FaceDetect(Mat image)
	:m_image(image), m_pResults(nullptr), m_pBuffer(nullptr)
{
	if (m_image.empty())
	{
		std::cout << "Can not load the image file " << std::endl;
		return;
	}
	m_pBuffer = (unsigned char*)malloc(DETECT_BUFFER_SIZE);
	if (!m_pBuffer)
	{
		std::cout << "Can not alloc buffer." << std::endl;
		return;
	}
}
FaceDetect::~FaceDetect()
{
	free(m_pBuffer);
}

std::vector<short*> FaceDetect::Detect()
{
	TickMeter cvtm;
	cvtm.start();

	m_pResults = facedetect_cnn(m_pBuffer, (unsigned char*)(m_image.ptr(0)), m_image.cols, m_image.rows, (int)m_image.step);
	cvtm.stop();

	std::cout << "time = " << cvtm.getTimeMilli() << " ms" << std::endl;
	std::cout << (m_pResults ? *m_pResults : 0) << " faces detected." << std::endl;

	m_resultImage = m_image.clone();
	/*
	for (int i = 0; i < (m_pResults ? *m_pResults : 0); i++)
	{
		short* p = ((short*)(m_pResults + 1)) + 142 * i;//每一个检测框的内容
		m_detectResults.push_back(p);
		int confidence = p[0];
		int x = p[1];
		int y = p[2];
		int w = p[3];
		int h = p[4];

		if (confidence <= 70)
			continue;
		//show the score of the face. Its range is [0-100]
		char sScore[256];
		snprintf(sScore, 256, "%d", confidence);
		cv::putText(m_resultImage, sScore, cv::Point(x, y - 3), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 255, 0), 1);
		//draw face rectangle
		rectangle(m_resultImage, Rect(x - w / 2, y - h, 2 * w, 2 * h), Scalar(0, 255, 0), 2);
		//rectangle(m_resultImage, Rect(x, y, w, h), Scalar(0, 255, 0), 2);
		//draw five face landmarks in different colors0
		cv::circle(m_resultImage, cv::Point(p[5], p[5 + 1]), 1, cv::Scalar(255, 0, 0), 2);
		cv::circle(m_resultImage, cv::Point(p[5 + 2], p[5 + 3]), 1, cv::Scalar(0, 0, 255), 2);
		cv::circle(m_resultImage, cv::Point(p[5 + 4], p[5 + 5]), 1, cv::Scalar(0, 255, 0), 2);
		cv::circle(m_resultImage, cv::Point(p[5 + 6], p[5 + 7]), 1, cv::Scalar(255, 0, 255), 2);
		cv::circle(m_resultImage, cv::Point(p[5 + 8], p[5 + 9]), 1, cv::Scalar(0, 255, 255), 2);

		cv::circle(m_resultImage, cv::Point(x, y), 1, cv::Scalar(0, 255, 255), 2);
		//print the result
		printf("face %d: confidence=%d, [%d, %d, %d, %d] (%d,%d) (%d,%d) (%d,%d) (%d,%d) (%d,%d)\n",
			i, confidence, x, y, w, h,
			p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14]);//五个点坐标
	}
	*/

	//imshow("face detection", m_resultImage);
	//waitKey();
	return m_detectResults;
}

Mat FaceDetect::Processed()
{
	return m_resultImage;
}