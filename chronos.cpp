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

#ifdef _WIN32
#include <windows.h>
#else
#include <ctime>
#include <sys/time.h>
#endif

#include "chronos.h"

Chronos::
Chronos() { 
    reset(); 
}

void 
Chronos::
reset(void) { 
    m_reset = time(); 
}

double 
Chronos::
elapsed(void) { 
    return time() - m_reset; 
}

double 
Chronos::time(void) {
#ifdef _WIN32
    LARGE_INTEGER counter, freq;
    QueryPerformanceCounter(&counter);
    QueryPerformanceFrequency(&freq);
    return (1.0*counter.QuadPart)/(1.0*freq.QuadPart);
#else
    struct timeval v;
    gettimeofday(&v, (struct timezone *) NULL);
    return v.tv_sec + v.tv_usec/1.0e6;
#endif
}
