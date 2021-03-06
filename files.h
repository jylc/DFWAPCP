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

#ifndef FILES_H
#define FILES_H

#ifdef _MSC_VER
  #include <io.h>
  #include <direct.h>
  #define GET_WORKDIR _getcwd
  #define OPEN_FILE _open
  #define LSEEK_FD _lseek
  #define FOPEN_RO(_descriptor,_filename) \
    if(fopen_s(&(_descriptor),_filename,"r")){ \
        fprintf(stderr,"unable to open file '%s'\n",_filename); \
        exit(-1); \
    }
#else
  #include <cstdio>     // for fopen and fclose
  #include <unistd.h>   // for getcwd, getpwuid, getuid
  #define GET_WORKDIR getcwd
  #define OPEN_FILE open
  #define LSEEK_FD lseek
  #define FOPEN_RO(_descriptor,_filename) \
    _descriptor=fopen(_filename,"r"); \
    if((_descriptor)==NULL){ \
        fprintf(stderr,"unable to open file '%s'\n",_filename); \
        exit(-1); \
    }
#endif

#endif // FILES_H
