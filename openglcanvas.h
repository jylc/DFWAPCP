//--------------------------------------------------
// Panoramic is an interface for the visualization of panoramas capable of
// handling wide fields of view, based on Möbius transformations.
// Copyright (C) 2015 Luis Peñaranda, Luiz Velho and Leonardo Sacht.
// 
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
// 
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------- 

#ifndef OPENGLCANVAS_H
#define OPENGLCANVAS_H

#include <QGLWidget>
#include <QTimer>
#include <QFile>
#include <QFileDialog>
#include <fcntl.h>
#include <opencv2/opencv.hpp>
#include "chronos.h"
#include <QtGui/qopenglfunctions.h>
#include "effectdrawing.h"
#include "image_data.h"
#include <vector>
#include <string>

class OpenGLCanvas:
        public QGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT
public:
    explicit OpenGLCanvas(QWidget *parent = 0);
    static const char* image_types;

protected:
    void read_config_file();
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();
    void define_texture_coordinates(float *texCoord, int m, int n, float min_phi, float max_phi, float min_lambda, float max_lambda);
    void load_sphere_mesh(float *positions, int m, int n,float radio);
    float calculate_extent(float fov_rads);
    void define_triangle_indices(unsigned int * indices, int m, int n);
private:
    int compute_auto_fov_max(int);
    void compute_scale();
    void load_image(std::string new_image);


    char *textFileRead(char *fn);
    void setShaders();
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);

signals:
    void fps(QString newFPS);
    void fov_changed(int new_fov);
    void max_fov_changed(int new_max_fov);

public slots:
    void show_effected_imgs();

    void automaxbutton(bool);
    void shrinkallbutton(bool);
    void change_fov(double f);
    void change_fov(int new_fov);
    void change_fov_max(int new_fov_max);
    void change_p_d(int new_p_d);
    void change_zb_lambda(int new_zb_lambda);
    void change_zb_R(int new_zb_R);
//    void change_scale(double s);
    void change_center_lambda(double lambda);
    void change_center_phi(double phi);
    void re_center();
    void change_fov_scale_relation(QString name);
    void change_visualization(QString name);
    void change_input_image();
    void slotTimer();

private:
    enum FileType get_file_type(const char*);

    int width;
    int height;
    int nChannels;

    double fov;
    double fov_max; // the \phi_{max} on the technote
    double scale;
    double center_lambda;
    double center_phi;
    QString fov_scale_relation;
    QString visualization;
    bool auto_fov_max;
    bool shrink_for_all;
    float pd; // parameter of the Pannini projection
    float zblambda,zbR; // parameters of the Zorin-Barr transformation

    Chronos time_time;
    QTimer time_timer;
    int time_frames;
    double time_start;
    double time_fps;

    unsigned int numberOfIndices;
    unsigned int * triangleIndices;
    float * verticesPositions;
    float * textureCoordinates;
    //float windowWidth;
    //float windowHeight;

    char* shader_dir;
    char* input_image_file;
    char* input_image_dir;
    std::vector<std::string> img_info;
    //ImageData img_data;
    std::shared_ptr<ImageData> img_data_ptr;
    float focal_length;
    float d;

    // input image size
    int image_size_x;
    int image_size_y;

    QPoint lastPos; // mouse click position
};

#endif // OPENGLCANVAS_H
