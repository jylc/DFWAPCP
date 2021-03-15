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
#include <QtGui>
#include "openglcanvas.h"
#include <cmath>
#include <sys/types.h>  // for fstat, getpwuid, getuid
#include <sys/stat.h>   // for fstat
#ifdef _WIN32
#include <stdlib.h>   // for getenv (TODO: use getenv_s)
#else
#include <pwd.h>      // for getpwuid
#endif
#include <cstdio>       // for getc
#include <cstring>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "files.h"

#define VERT_SHADER_FILE "vertex_shader.vert"
#define FRAG_SHADER_FILE "fragment_shader.frag"

#define CONST_PI        3.141592653589793115997963468544185161590576171875
#define CONST_PI_2      1.5707963267948965579989817342720925807952880859375
#define CONST_PI_F      3.1415927410125732421875f
#define CONST_PI_2_F    1.57079637050628662109375f
//#define CONST_PI        (0x1.921fb54442d18p+1)
//#define CONST_PI_2      (0x1.921fb54442d18p+0)
//#define CONST_PI_F      (0x1.921fb6p+1f)
//#define CONST_PI_2_F    (0x1.921fb6p+0f)

// These definitions specify which attributes should be used to store some
// parameters passed to the vertex shader. Nvidia hardware
// only leaves attributes 1 and 7 unreserved; attributes 8 to 15 are
// reserved for textures.
// TODO: we use attributes 14 and 15, which work for Nvidia; we need to
// test with other hardware.

//传递参数
//焦距
#define FL_ATTR 11 
//d=min(W,H)
#define DT_ATTR 12
#define PD_ATTR 13
#define ZBL_ATTR 14
#define ZBR_ATTR 15


const char* OpenGLCanvas::image_types = "All images (*.png *.jpg *.jpeg *.pnm *.pbm *.pgm *.ppm);;PNG images (*.png);;JPG images (*.jpg *.jpeg);;PBM images (*.pnm *.pbm *.pgm *.ppm);;All files (*)";


OpenGLCanvas::OpenGLCanvas(QWidget* parent) :
    QGLWidget(parent)
{
    setFormat(QGL::DoubleBuffer | QGL::DepthBuffer);
    fov = 60.f;
    fov_max = 60.f; // TODO: I only use floats here because Leo did... check
    scale = 1.0f;
    center_lambda = 0.f;
    center_phi = 0.f;
    fov_scale_relation = "Simple";
    visualization = "Zorin-Barr";
    auto_fov_max = false;
    shrink_for_all = false;
    pd = 0.f; // Pannini projection d
    zblambda = .1f; // Zorin-Barr transformation lambda
    zbR = 1.f; // Zorin-Barr transformation R
    //TODO: 待改
    focal_length = 1.f;    //自定焦距

    time_frames = 0;
    time_timer.setInterval(0);
    connect(&time_timer, SIGNAL(timeout()), this, SLOT(slotTimer()));
    time_start = time_time.time();
}

void OpenGLCanvas::slotTimer(void) {
    updateGL();
}

void OpenGLCanvas::change_fov(double f) {

    if (fov != f && f >= 1.f && f <= 360.f)
        fov = f;
    compute_scale();
    //    scale = 0.3f;

    fprintf(stderr, "change fov, fov=%f, fov_max=%f, scale=%f\n", fov, fov_max, scale);
    emit fov_changed((int)fov);

    updateGL();

}

void OpenGLCanvas::automaxbutton(bool autovalue) {
    auto_fov_max = autovalue;
    change_fov_max(compute_auto_fov_max(fov));
}

void OpenGLCanvas::shrinkallbutton(bool shrinkvalue) {
    shrink_for_all = shrinkvalue;
    compute_scale();
    updateGL();
}

void OpenGLCanvas::change_fov(int new_fov) {
    if (new_fov <= 360 && new_fov >= 1) {
        change_fov((double)new_fov);
        if (auto_fov_max)
            change_fov_max(compute_auto_fov_max(new_fov));
    }
}

void OpenGLCanvas::change_fov_max(int new_fov_max) {
    if (new_fov_max <= 360.f && new_fov_max >= 1)
        fov_max = (double)new_fov_max;
    compute_scale();
    fprintf(stderr, "change fov_max, fov=%f, fov_max=%f, new scale=%f\n", fov, fov_max, scale);
    emit max_fov_changed((int)fov_max);
    updateGL();
}

void OpenGLCanvas::change_p_d(int new_p_d) {
    pd = (float)new_p_d / 10;
    fprintf(stderr, "p_d=%f\n", pd);
    glVertexAttrib1f(PD_ATTR, pd);
    updateGL();
}

void OpenGLCanvas::change_zb_lambda(int new_zb_lambda) {
    zblambda = (float)new_zb_lambda / 1000;
    fprintf(stderr, "zb_lambda=%f\n", zblambda);
    glVertexAttrib1f(ZBL_ATTR, zblambda);
    updateGL();
}

void OpenGLCanvas::change_zb_R(int new_zb_R) {
    zbR = (float)new_zb_R / 100;
    fprintf(stderr, "zb_R=%f\n", zbR);
    glVertexAttrib1f(ZBR_ATTR, zbR);
    updateGL();
}

void OpenGLCanvas::change_center_lambda(double lambda) {

    if (center_lambda != lambda && lambda >= -CONST_PI_F && lambda <= CONST_PI_F) {
        center_lambda = lambda;
        updateGL();
    }

}

void OpenGLCanvas::change_center_phi(double phi) {

    if (center_phi != phi && phi >= -CONST_PI_2_F && phi <= CONST_PI_2_F) {
        center_phi = phi;
        updateGL();
    }

}

void OpenGLCanvas::re_center() {
    center_phi = .0f;
    center_lambda = .0f;
    updateGL();
}

void OpenGLCanvas::change_fov_scale_relation(QString name) {

    fov_scale_relation = name;
    compute_scale();
    fprintf(stderr, "changed scale relation, scale=%f, fov_max=%f\n", scale, fov_max);
    updateGL();

}

void OpenGLCanvas::change_visualization(QString name) {
    visualization = name;
    compute_scale();
    updateGL();
}

// This function reads the contents of the ~/.panorc file and stores the
// options in private variables.
void OpenGLCanvas::read_config_file() {
   
    char* filepath = (char*)malloc(512 * sizeof(char));
    strcpy(filepath, "G:\\VSProject\\DFWAPCP");
    strcat(filepath, "\\.panorc");
    shader_dir = (char*)malloc(512 * sizeof(char));
    shader_dir[0] = '\0';
    input_image_file = (char*)malloc(512 * sizeof(char));
    input_image_file[0] = '\0';
    input_image_dir = (char*)malloc(512 * sizeof(char));
    input_image_dir[0] = '\0';
    char* read_line = (char*)malloc(64 * sizeof(char));
    struct stat testbuf;
    if (stat(filepath, &testbuf)) {
        fprintf(stderr, "%s does not exist\n", filepath);
    }
    else {
        FILE* rcfile;
        FOPEN_RO(rcfile, filepath);
        char c;
        char* line = (char*)malloc(512 * sizeof(char));
        while ((c = getc(rcfile)) != EOF) {
            while (c == '\n') // discard empty lines
                c = getc(rcfile);
            if (c == EOF)
                break;
            line[0] = c; // first char on the line was already read
            if (!fgets(line + 1, 511, rcfile)) {
                fprintf(stderr, "error reading rcfile\n");
                exit(-1);
            }
            // check for 'shader_dir' option
            if (!strncmp(line, "shader_dir=", 11)) {
                strcpy(shader_dir, line + 11);
                shader_dir[strlen(line) - 12] = '\0';
                fprintf(stderr, "shader_dir=%s\n", shader_dir);
            }
            // check for 'image_file' option
            if (!strncmp(line, "image_file=", 11)) {
                strcpy(input_image_file, line + 11);
                input_image_file[strlen(line) - 12] = '\0';
                fprintf(stderr, "input_image_file=%s\n", input_image_file);
            }
            // check for 'image_dir' option
            if (!strncmp(line, "image_dir=", 10)) {
                strcpy(input_image_dir, line + 10);
                input_image_dir[strlen(line) - 11] = '\0';
                fprintf(stderr, "input_image_dir=%s\n", input_image_dir);
            }
            // check for 'max_fov' option
            if (!strncmp(line, "max_fov=", 8)) {
                strcpy(read_line, line + 8);
                read_line[strlen(line) - 9] = '\0';
                fov_max = atof(read_line);
                fprintf(stderr, "max_fov=%f\n", fov_max);
            }
            // check for 'fov' option
            if (!strncmp(line, "fov=", 4)) {
                strcpy(read_line, line + 4);
                read_line[strlen(line) - 5] = '\0';
                fov = atof(read_line);
                fprintf(stderr, "fov=%f\n", fov);
            }
            // check for 'auto_max_fov' option
            if (!strncmp(line, "auto_max_fov=", 13)) {
                strcpy(read_line, line + 13);
                read_line[strlen(line) - 14] = '\0';
                auto_fov_max = atof(read_line);
                fprintf(stderr, "auto_max_fov=%d\n", auto_fov_max);
            }
            // check for the Pannini parameter, d
            if (!strncmp(line, "pd=", 3)) {
                strcpy(read_line, line + 3);
                read_line[strlen(line) - 4] = '\0';
                pd = atof(read_line);
                fprintf(stderr, "pd=%f\n", pd);
            }
            // check for the Zorin-Barr transformation parameters, lambda and R
            if (!strncmp(line, "zblambda=", 9)) {
                strcpy(read_line, line + 9);
                read_line[strlen(line) - 10] = '\0';
                zblambda = atof(read_line);
                fprintf(stderr, "zblambda=%f\n", zblambda);
            }
            if (!strncmp(line, "zbR=", 4)) {
                strcpy(read_line, line + 4);
                read_line[strlen(line) - 5] = '\0';
                zbR = atof(read_line);
                fprintf(stderr, "zbR=%f\n", zbR);
            }
        }
        fclose(rcfile);
    }
    free(filepath);
    emit fov_changed((int)fov);
    emit max_fov_changed((int)fov_max);
}

void OpenGLCanvas::load_image(std::string new_image) {
    // 加载并生成纹理

    unsigned char* data = 
        stbi_load(new_image.c_str(), &width, &height, &nChannels,0);
    d = 1.0f * MIN(width, height);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    GLuint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D,
        0,
        GL_RGB,
        width,
        height,
        0,
        GL_RGB,
        GL_UNSIGNED_BYTE,
        data);

    //分割图片路径
    size_t len = new_image.rfind('/');
    strcpy(input_image_file, new_image.c_str() + len + 1);
    img_info.clear();
    img_info.push_back(input_image_dir);
    img_info.push_back(input_image_file);
    std::string tmp(input_image_file);
    size_t pos = tmp.find_first_of('.');
    input_image_file[pos] = '\0';
    char* s = strcat(input_image_file, ".png");
    img_info.push_back("G:/VSProject/DFWAPCP/data/masks/");
    img_info.push_back(s);

    img_data_ptr = std::make_shared<ImageData>(img_info[0], img_info[1], img_info[2], img_info[3], focal_length);
    fprintf(stderr, "input_image_file : %s\n", input_image_file);
    
}

void OpenGLCanvas::change_input_image() {
    std::string fname = QFileDialog::getOpenFileName(this, tr("Choose Panorama File"), input_image_dir, OpenGLCanvas::image_types).toStdString();
    if (fname.length()>0) {
        load_image(fname);
        updateGL();
    }
}

void OpenGLCanvas::initializeGL() {
    initializeOpenGLFunctions();
    glShadeModel(GL_SMOOTH);
    glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
    glClearDepth(1.0f);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

#ifndef _MSC_VER
#ifdef __APPLE__
    const char* const progname = "PROJ_ROOT_DIR";
#else
    // progname is a file name or a path???
    const char* const progname = (char*)(PROGNAME);
#endif
    fprintf(stderr, "progname=%s\n", progname);
#endif // MSC_VER

    read_config_file();
    compute_scale();
    // If the input file does not exist or was not specified.
    struct stat testbuf;
    if (stat(input_image_file, &testbuf) || !strcmp(input_image_file, "")) {
        load_image(QFileDialog::getOpenFileName(this, tr("Choose Panorama File"), input_image_dir, OpenGLCanvas::image_types).toStdString().c_str());
    }
    else {
    load_image(input_image_file);
   }

    // mesh resolution
    int m, n;
    //m = 100;
    //n = 100;
    auto count = img_data_ptr->getCountOfWAndH();
    n = count[0];//宽
    m = count[1];//高

    //defining texture coordinates
    int meshNumTexCoord = m * n;
    float* texCoord = (float*)malloc(2 * meshNumTexCoord * sizeof(float));
    if (texCoord == NULL) {
        printf("problem allocating memory for texture coordinates \n");
    }
    define_texture_coordinates(texCoord, m, n, -CONST_PI_2_F, CONST_PI_2_F, -CONST_PI_F, CONST_PI_F);

    //defining positions of the sphere vertices
    int meshNumVertices = m * n;
    float* positions = (float*)malloc(3 * meshNumVertices * sizeof(float));
    if (positions == NULL) {
        printf("problem allocating memory for positions \n");
    }
    //float fov_rads = (fov / 360.f) * CONST_PI_F;
    //vertex_transformation(positions, m, n, center_lambda, center_phi, fov_rads, scale); //passar pelo vertex shader
    float radio = (1.f * height) / (1.f * width);
    load_sphere_mesh(positions, m, n, radio); //colocar essa e funcoes para textura e triangulos no initializeGL

    //defining triagle indices
    unsigned int meshNumFaces = 2 * (m - 1) * (n - 1);
    unsigned int meshNumIndices = 3 * meshNumFaces;
    unsigned int* indices = (unsigned int*)malloc(meshNumIndices * sizeof(unsigned int));
    define_triangle_indices(indices, m, n);

    // draw setup
    verticesPositions = positions;
    textureCoordinates = texCoord;
    numberOfIndices = meshNumIndices;
    triangleIndices = indices;

    setShaders();
}


void OpenGLCanvas::define_texture_coordinates(float* texCoord,int m, int n, float min_phi, float max_phi, float min_lambda, float max_lambda) {
    float delta_x = 1.0f / (1.0 * (n - 1));
    float delta_y = 1.0f / (1.0 * (m - 1));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            texCoord[2 * (j + i * n)] = delta_x * j;
            texCoord[2 * (j + i * n) + 1] = delta_y * i;
        }
    }
}

// This function makes the same computation GLSL does. It is never called.
void OpenGLCanvas::vertex_transformation(float* positions, int m, int n, float center_lambda, float center_phi, float fov_rads, float scale) {

    float min_lambda = -CONST_PI_F;
    float max_lambda = CONST_PI_F;
    float min_phi = -CONST_PI_2_F;
    float max_phi = CONST_PI_2_F;

    float delta_lambda = (max_lambda - min_lambda) / (1.0 * (n - 1));
    float delta_phi = (max_phi - min_phi) / (1.0 * (m - 1));

    float lambda, phi, x, y, z, u, v, r, theta;

    //calculating the extent of the projection for the given FOV
    lambda = fov_rads;
    phi = 0.f;
    // OpenGL: x is the vertical axes pointg downwards, and y is horizontal axes
    y = sinf(phi);
    x = -sinf(lambda) * cosf(phi);
    z = -cosf(lambda) * cosf(phi);
    u = 2.f * x / (1.f - z);
    v = 2.f * y / (1.f - z);
    r = hypotf(u, v);
    theta = atan2f(u, v);
    r *= scale;
    u = -r * sinf(theta);
    v = r * cosf(theta);
    x = (4.f * u) / (u * u + v * v + 4.f);
    y = (4.f * v) / (u * u + v * v + 4.f);
    z = (u * u + v * v - 4.f) / (u * u + v * v + 4.f);
    u = x / (-z);
    v = y / (-z);
    float extent = u;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {

            lambda = (min_lambda + delta_lambda * j);
            phi = (min_phi + delta_phi * i);

            // OpenGL: x is the vertical axes pointg downwards, and y is horizontal axes
            y = sinf(phi);
            x = -sinf(lambda) * cosf(phi);
            z = -cosf(lambda) * cosf(phi);

            //Rotation 1: (-center_lambda)-rotation on the xz-plane
            float x_copy = x;
            x = cosf(-center_lambda) * x - sinf(-center_lambda) * z;
            //y=1.f*y;
            z = sinf(-center_lambda) * x_copy + cosf(-center_lambda) * z;

            //Rotation 2: (-center_phi)-rotation on the yz-plane
            float y_copy = y;
            //x = 1.f*x;
            y = cosf(-center_phi) * y - sinf(-center_phi) * z;
            z = sinf(-center_phi) * y_copy + cosf(-center_phi) * z;

            u = 2.f * x / (1.f - z);
            v = 2.f * y / (1.f - z);

            r = hypotf(u, v);
            theta = atan2f(u, v);

            // scaling the complex plane according to scale specified in the interface (relate it to FOV)
            r *= scale;

            u = -r * sinf(theta);
            v = r * cosf(theta);

            x = (4.f * u) / (u * u + v * v + 4.f);
            y = (4.f * v) / (u * u + v * v + 4.f);
            z = (u * u + v * v - 4.f) / (u * u + v * v + 4.f);

            lambda = atan2f(x, -z) / CONST_PI_F;
            phi = asinf(y) / CONST_PI_2_F;

            if (visualization == "Moebius" || visualization == "Perspective") {
                u = x / (-z);
                v = y / (-z);
                positions[3 * (j + i * n)] = u / extent;
                positions[3 * (j + i * n) + 1] = v / extent;
                positions[3 * (j + i * n) + 2] = z;
            }
            else if (visualization == "3D Sphere") {
                positions[3 * (j + i * n)] = 0.9f * x;
                positions[3 * (j + i * n) + 1] = 0.9f * y;
                positions[3 * (j + i * n) + 2] = z;
            }
            else if (visualization == "Equi-Rectangular") {
                positions[3 * (j + i * n)] = lambda;
                positions[3 * (j + i * n) + 1] = phi;
                positions[3 * (j + i * n) + 2] = z;
            }
            else if (visualization == "Stereographic") {
                u = 2 * x / (-z + 1);
                v = 2 * y / (-z + 1);
                positions[3 * (j + i * n)] = u / extent;
                positions[3 * (j + i * n) + 1] = v / extent;
                positions[3 * (j + i * n) + 2] = z;
            }
            else if (visualization == "Mercator") {
                u = lambda;
                v = logf((1.0 / cosf(phi)) + tanf(phi));
                positions[3 * (j + i * n)] = u / extent;
                positions[3 * (j + i * n) + 1] = v / extent;
                positions[3 * (j + i * n) + 2] = z;
            }
            else if (visualization == "Orthographic") {
                u = x;
                v = y;
                positions[3 * (j + i * n)] = x / extent;
                positions[3 * (j + i * n) + 1] = y / extent;
                positions[3 * (j + i * n) + 2] = z;
            }
            else if (visualization == "Pannini") {
                u = 2 * x / (-z + pd);
                v = 2 * y / (-z + pd);
                positions[3 * (j + i * n)] = u / extent;
                positions[3 * (j + i * n) + 1] = v / extent;
                positions[3 * (j + i * n) + 2] = z;
            }
            else if (visualization == "Zorin-Barr") {
                // perspective
                u = x / (-z);
                v = y / (-z);
                // apply Z-B transformation to (u,v)
                float lambda = .1f;
                float R = 1.f;
                float alpha = atanf(v / u);
                float r = hypotf(u, v);
                float rhoprime = (lambda * r / R) + (1.f - lambda) * (R * (sqrtf(r * r + 1.f) - 1.f)) / (r * (sqrtf(R * R + 1.f) - 1.f));
                u = rhoprime * cosf(alpha);
                v = rhoprime * sinf(alpha);
                //
                positions[3 * (j + i * n)] = u / extent;
                positions[3 * (j + i * n) + 1] = v / extent;
                positions[3 * (j + i * n) + 2] = z;
            }
        }
    }
}

void OpenGLCanvas::load_sphere_mesh(float* positions, int m, int n,float radio) {

    float delta_x = 2.0f / (1.0 * (n - 1));
    float delta_y = 2.0f / (1.0 * (m - 1)); 

    float lambda, phi, x, y, z;

    fprintf(stderr, "\nradio = %lf\n", radio);
    if (radio < 1) {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                positions[3 * (j + i * n)] = delta_x * j - 1.0f;
                positions[3 * (j + i * n) + 1] = (delta_y * i - 1.0f) * radio;
                positions[3 * (j + i * n) + 2] = -1;
            }
        }
    }
    else {
        radio = 1 / radio;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                positions[3 * (j + i * n)] = (delta_x * j - 1.0f) * radio;
                positions[3 * (j + i * n) + 1] = delta_y * i - 1.0f;
                positions[3 * (j + i * n) + 2] = -1;
            }
        }
    }

}

float OpenGLCanvas::calculate_extent(float fov_rads) {
    float lambda, phi, x, y, z, u, v, r, theta;
    //calculating the extent of the projection for the given FOV
    lambda = fov_rads;
    phi = 0.f;
    // OpenGL: x is the vertical axes pointg downwards, and y is horizontal axes
    y = sinf(phi);
    x = -sinf(lambda) * cosf(phi);
    z = -cosf(lambda) * cosf(phi);
    u = 2.f * x / (1.f - z);
    v = 2.f * y / (1.f - z);
    r = hypotf(u, v);//sqrt(u*u+v*v);
    theta = atan2f(u, v);
    r *= scale;
    u = -r * sinf(theta);
    v = r * cosf(theta);
    x = (4.f * u) / (u * u + v * v + 4.f);
    y = (4.f * v) / (u * u + v * v + 4.f);
    z = (u * u + v * v - 4.f) / (u * u + v * v + 4.f);
    if (visualization == "Moebius" || visualization == "Perspective") {
        u = x / (-z);
        v = y / (-z);
    }
    else if (visualization == "Stereographic") {
        u = 2.f * x / (1.f - z);
        v = 2.f * y / (1.f - z);
    }
    else if (visualization == "Mercator") {
        u = fov_rads;
        // Warning: this extent calculation is wrong.
        // Write now it's olny showing the entire panorama intead of
        // the corresponging FOV.
    }
    else if (visualization == "Orthographic") {
        if (z < 0.f) {
            u = x;
            v = y;
        }
        else {
            u = v = 1.f;
        }
    }
    else if (visualization == "Pannini") {
        u = 2.f * x / (pd - z);
        v = 2.f * y / (pd - z);
    }
    else if (visualization == "Zorin-Barr") {
        // TODO: check whether this is correct
        u = x / (-z);
        v = y / (-z);
    }
    return u;
}

void OpenGLCanvas::define_triangle_indices(unsigned int* indices, int m, int n) {

    for (int i = 0; i < m - 1; i++) {
        for (int j = 0; j < n - 1; j++) {

            unsigned int index = (j + i * n);

            indices[3 * (2 * (j + i * (n - 1)))] = index;
            indices[3 * (2 * (j + i * (n - 1))) + 1] = index + 1;
            indices[3 * (2 * (j + i * (n - 1))) + 2] = index + n;

            indices[3 * (2 * (j + i * (n - 1)) + 1)] = index + 1;
            indices[3 * (2 * (j + i * (n - 1)) + 1) + 1] = index + n + 1;
            indices[3 * (2 * (j + i * (n - 1)) + 1) + 2] = index + n;

        }
    }

}

// This function computes a new value of the maximum fov and must be called
// when the fov is changed and the auto setting is enabled.
int OpenGLCanvas::compute_auto_fov_max(int new_fov) {
    if (new_fov < 60)
        return 60;
    if (new_fov > 179)
        return 1;
    //if(new_fov>60)
    return (90 - new_fov / 2);
}

// This function computes the scale, in function of the projection method,
// the fov/scale ratio and the value of shrink_for_all.
void OpenGLCanvas::compute_scale() {
    if (shrink_for_all || visualization == "Moebius" || visualization == "3D Sphere") {
        if (fov <= fov_max) {
            scale = 1.f;
            //else if (fov>295.f)
            //    scale = 0.02f; // TODO: check this value wrt fov_max
        }
        else {
            if (fov_scale_relation == "Simple")
                scale = fov_max / fov;
            else if (fov_scale_relation == "Square Root")
                scale = sqrtf((360.f - fov_max - fov) / (360.f - 2.f * fov_max));
            else if (fov_scale_relation == "Linear")
                scale = (360.f - fov_max - fov) / (360.f - 2.f * fov_max);
            else if (fov_scale_relation == "Square Power")
                scale = powf((360.f - fov_max - fov) / (360.f - 2.f * fov_max), 2.f);
            else if (fov_scale_relation == "Cubic Power")
                scale = powf((360.f - fov_max - fov) / (360.f - 2.f * fov_max), 3.f);
            else if (fov_scale_relation == "Logarithm")
                scale = logf(expf(1.f) + (1.f - expf(1.f)) * (fov - fov_max) / (360.f - 2.f * fov_max));
        }
    }
    else {
        scale = 1.f;
    }
}

void OpenGLCanvas::resizeGL(int w, int h) {
    if (w > h)
        glViewport(0, (h - w) / 2, w, w);
    else
        glViewport((w - h) / 2, 0, h, h);
}

char* OpenGLCanvas::textFileRead(char* fn) {

    FILE* fp;
    char* content = NULL;
    int f, count;
    f = OPEN_FILE(fn, O_RDONLY);
    count = LSEEK_FD(f, 0, SEEK_END);
    //    close(f);
    if (fn != NULL) {
        FOPEN_RO(fp, fn);
        if (fp != NULL) {
            if (count > 0) {
                content = (char*)malloc(sizeof(char) * (count + 1));
                count = fread(content, sizeof(char), count, fp);
                content[count] = '\0';
            }
            fclose(fp); // maybe this line must be outside the {}
        }
    }
    return content;
}

void OpenGLCanvas::setShaders() {

    char* vs, * fs;
    fprintf(stderr, "1.GLEW initialization.\n");
#ifndef __APPLE__
#ifdef GLEW_VERSION_1_5
    fprintf(stderr, "2.GLEW initialization.\n");
    GLenum err = glewInit();
    if (err != GLEW_OK) {
        fprintf(stderr, "error in GLEW initialization: %s\n", glewGetString(err));
        exit(-1);
    }
#endif
#endif

    GLuint v = glCreateShader(GL_VERTEX_SHADER);
    GLuint f = glCreateShader(GL_FRAGMENT_SHADER);

    // Configure vertex and fragment shader files.
    char* vs_file = (char*)malloc(512 * sizeof(char*));
    char* fs_file = (char*)malloc(512 * sizeof(char*));
    if (!strcmp(shader_dir, "")) { // if shader_dir was not configured
        if (!GET_WORKDIR(vs_file, 512) || !GET_WORKDIR(fs_file, 512)) {
            fprintf(stderr, "error reading shader files\n");
            exit(-1);
        }
        strcat(vs_file, "/shaders/");
        strcat(fs_file, "/shaders/");
    }
    else {
        strcpy(vs_file, shader_dir);
        strcat(vs_file, "/");
        strcpy(fs_file, shader_dir);
        strcat(fs_file, "/");
    }
    strcat(vs_file, VERT_SHADER_FILE);
    strcat(fs_file, FRAG_SHADER_FILE);
    fprintf(stderr, "vs_file=%s\nfs_file=%s\n", vs_file, fs_file);

    struct stat vs_testbuf, fs_testbuf;
    if (stat(vs_file, &vs_testbuf) || stat(fs_file, &fs_testbuf)) {
        fprintf(stderr, "a shader file does not exist!\n");
        free(vs_file);
        free(fs_file);
        exit(-1);
    }

    vs = textFileRead(vs_file);
    fs = textFileRead(fs_file);

    const char* vv = vs;
    const char* ff = fs;

    glShaderSource(v, 1, &vv, NULL);
    glShaderSource(f, 1, &ff, NULL);

    free(vs); free(fs);
    free(vs_file); free(fs_file);

    glCompileShader(v);
    glCompileShader(f);

    GLuint p = glCreateProgram();

    // Bind attributes pd, zblambda and zbR to the vertex shader.
    glVertexAttrib1f(DT_ATTR, d);
    glBindAttribLocation(p, DT_ATTR, "d");
    glVertexAttrib1f(FL_ATTR, focal_length);
    glBindAttribLocation(p, FL_ATTR, "focal_length");

    glVertexAttrib1f(PD_ATTR, pd);
    glBindAttribLocation(p, PD_ATTR, "pd");
    glVertexAttrib1f(ZBL_ATTR, zblambda);
    glBindAttribLocation(p, ZBL_ATTR, "zblambda");
    glVertexAttrib1f(ZBR_ATTR, zbR);
    glBindAttribLocation(p, ZBR_ATTR, "zbR");

    glAttachShader(p, v);
    glAttachShader(p, f);
    
    glLinkProgram(p);
    glUseProgram(p);
}

void OpenGLCanvas::mousePressEvent(QMouseEvent* event) {
    lastPos = event->pos();
    fprintf(stderr,"mouse click\n");
}

void OpenGLCanvas::mouseMoveEvent(QMouseEvent* event) {
    // scroll with the left button
    if (event->buttons() == Qt::LeftButton) {
        fprintf(stderr, "mouse left button click\n");
        // compute the delta and move the image
        center_lambda += (event->x() - lastPos.x()) * CONST_PI_F / width;
        center_phi += (event->y() - lastPos.y()) * CONST_PI_F / height;
        lastPos = event->pos();
        updateGL();
    }
}

void OpenGLCanvas::wheelEvent(QWheelEvent* event) {
    if (event->orientation() == Qt::Vertical) {
        if (event->modifiers() == Qt::ShiftModifier) {
            change_fov_max(fov_max + ((double)event->delta()) / 30);
        }
        else {
            int new_fov = fov + event->delta() / 30;
            change_fov((double)new_fov);
            if (auto_fov_max)
                change_fov_max(compute_auto_fov_max(new_fov));
        }
    }
}

void OpenGLCanvas::paintGL() {

    float fov_rads = (fov / 360.f) * CONST_PI_F;
  
    // defining transformation parameters (that will be passed to the vertex shader)
    float extent = calculate_extent(fov_rads);
    float vis_mode = .0f;
    if (visualization == "Moebius" || visualization == "Perspective") vis_mode = 1.f;
    else if (visualization == "3D Sphere") vis_mode = 2.f;
    else if (visualization == "Equi-Rectangular") vis_mode = 3.f;
    else if (visualization == "Stereographic") vis_mode = 4.f;
    else if (visualization == "Mercator") vis_mode = 5.f;
    else if (visualization == "Orthographic") vis_mode = 4.5f;
    else if (visualization == "Pannini") vis_mode = 4.25f;
    else if (visualization == "Zorin-Barr") vis_mode = 6.f;

    fprintf(stderr, "center_lambda = %f, center_phi = %f , fov = %f ,fov_max = %f , scales = %f , visualization = %f\n", center_lambda, center_phi,fov,fov_max,scale, vis_mode);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, 2.0 / extent, 0.0, 2.0 / scale, 0.0, -2.0 / vis_mode);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glOrtho(0.0, 2.0 / center_lambda, 0.0, 2.0 / center_phi, -1.0, 1.0);
    
    // drawing the mesh
    glClearColor(1.0, 1.0, 1.0, 1.0);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

    glColor3f(1, 0, 0);

    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, verticesPositions);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    glTexCoordPointer(2, GL_FLOAT, 0, textureCoordinates);

    glDrawElements(GL_TRIANGLES, numberOfIndices, GL_UNSIGNED_INT, triangleIndices);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);

    time_frames++;
    //    if (time_frames > 0) {
    double dt = time_time.elapsed();
    //        if (dt > 0.5) {
    time_fps = time_frames / dt;
    time_frames = 0;
    time_time.reset();
    emit fps(QString("%1 fps").arg((int)(time_fps + 0.5)));
}

void OpenGLCanvas::show_effected_imgs()
{
    effectdrawing* drawing= new effectdrawing(img_data_ptr);
    drawing->show();
}

