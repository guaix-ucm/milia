/*
 * Copyright 2008-2010 Sergio Pascual
 * 
 * This file is part of Milia
 * 
 * Milia is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Milia is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Milia.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

// $Id$

#include "exception.h"

namespace milia
{

    exception::exception(const std::string& information) :
            m_message(information)
    {}

    exception::~exception() throw()
    {}

    const char* exception::what() const throw()
    {
        return m_message.c_str();
    }


    recollapse::recollapse(const std::string& information) :
            exception(information)
    {}

    recollapse::~recollapse() throw()
    {}


    no_big_bang::no_big_bang(const std::string& information) :
            exception(information)
    {}

    no_big_bang::~no_big_bang() throw()
    {}



}
