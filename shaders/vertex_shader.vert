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

#version 120

varying vec2 texcoord0;
vec4 pos;
float u, v, x, y, z;
varying float r, theta, s;
float lambda, phi;
float extent, scale, vis_mode, center_lambda, center_phi,d;
attribute float zblambda,zbR,pd;

void main(void){

    texcoord0 = vec2(gl_MultiTexCoord0);
    pos = gl_Vertex;

    // getting paramters from the interface
    extent = gl_ProjectionMatrix[0].x;
    scale = gl_ProjectionMatrix[1].y;
    vis_mode = gl_ProjectionMatrix[2].z;
    center_lambda = gl_ModelViewMatrix[0].x;
    center_phi = gl_ModelViewMatrix[1].y;

    x = pos.x;
    y = pos.y;
    z = pos.z;

    // perspective projection
    u=x/(-z);
    v=y/(-z);
    // Z-B transformation
    float zbalpha=atan(v,u);
    float zbr=sqrt(u*u+v*v);
    float zbrho=(zblambda*zbr/zbR)+(1.0-zblambda)*(zbR*(sqrt(zbr*zbr+1.0)-1.0))/(zbr*(sqrt(zbR*zbR+1.0)-1.0));
    //
    //float r0= d/(2.0*tan(0.5*atan(d,(2.0f*focal_length))));
    //float rp=zbr;
    //float ru = r0*tan(0.5*atan(rp,focal_lenght));
    u=zbrho*cos(zbalpha);
    v=zbrho*sin(zbalpha);
    //u=zbrho;
    //v=zbrho;
    gl_Position = vec4(u/extent,-v/extent,z,1.0);
}
