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

uniform sampler2D diffuse_texture;
varying vec2 texcoord0;
varying vec4 K_color;
varying float r, theta, s;

void main() {
    gl_FragColor = texture2D(diffuse_texture, texcoord0);
//    gl_FragColor =  0.1*texture2D(diffuse_texture, texcoord0) + s;
//    gl_FragColor = 0.1*texture2D(diffuse_texture, texcoord0) + theta/6.28+0.5;
}
