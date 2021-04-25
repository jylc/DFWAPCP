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


	m_resultImage = m_image.clone();
	
	for (int i = 0; i < (m_pResults ? *m_pResults : 0); i++)
	{
		short* p = ((short*)(m_pResults + 1)) + 142 * i;//每一个检测框的内容
		m_detectResults.push_back(p);
	}

	return m_detectResults;
}

Mat FaceDetect::Processed()
{
	return m_resultImage;
}