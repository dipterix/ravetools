/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
#ifndef __VCG_EXCEPTION_H
#define __VCG_EXCEPTION_H

#include <stdexcept>
#include <iostream>
// #include <Rcpp.h>

namespace vcg
{
class MissingComponentException : public std::runtime_error
{
public:
  MissingComponentException(const std::string &err):std::runtime_error(err)
  {
    ::Rf_warning("vcglib: Missing Component Exception - %s\n", err.c_str());
  }
    virtual const char *what() const throw ()
    {
      static char buf[128]="Missing Component";
      return buf;
    }
};

class MissingCompactnessException : public std::runtime_error
{
public:
  MissingCompactnessException(const std::string &err):std::runtime_error(err)
  {
    ::Rf_warning("vcglib: Lack of Compactness Exception - %s\n", err.c_str());
  }
    virtual const char *what() const throw ()
    {
      static char buf[128]="Lack of Compactness";
      return buf;
    }
};

class MissingTriangularRequirementException : public std::runtime_error
{
public:
  MissingTriangularRequirementException(const std::string &err):std::runtime_error(err)
  {
    ::Rf_warning("vcglib: Mesh has to be composed by triangle and not polygons - %s\n", err.c_str());
  }

    virtual const char *what() const throw ()
    {
      static char buf[128]="Mesh has to be composed by triangle and not polygons";
      return buf;
    }
};

class MissingPolygonalRequirementException : public std::runtime_error
{
public:
  MissingPolygonalRequirementException(const std::string &err):std::runtime_error(err)
  {
    ::Rf_warning("vcglib: Mesh has to be composed by polygonal faces (not plain triangles) - %s\n", err.c_str());
  }

    virtual const char *what() const throw ()
    {
      static char buf[128]="Mesh has to be composed by polygonal faces (not plain triangles) ";
      return buf;
    }
};

class MissingTetrahedralRequirementException : public std::runtime_error
{
public:
  MissingTetrahedralRequirementException(const std::string &err):std::runtime_error(err)
  {
    ::Rf_warning("vcglib: Mesh has to be composed by tetrahedras - %s\n", err.c_str());
  }

    virtual const char *what() const throw ()
    {
      static char buf[128]="Mesh has to be composed by tetrahedras";
      return buf;
    }
};

class MissingPreconditionException : public std::runtime_error
{
public:
  MissingPreconditionException(const std::string &err):std::runtime_error(err)
  {
    ::Rf_warning("vcglib: Mesh does not satisfy the following precondition: %s\n", err.c_str());
  }

    virtual const char *what() const throw ()
    {
      static char buf[128]="Mesh does not satisfy precondition";
      return buf;
    }
};

} // end namespace vcg
#endif // EXCEPTION_H
