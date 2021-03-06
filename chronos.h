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

#ifndef CHRONOS_H
#define CHRONOS_H

class Chronos {
public:
    Chronos();
    void reset(void);
    double elapsed(void);
    double time(void);
private:
    double m_reset;
};

#endif // CHRONOS_H
