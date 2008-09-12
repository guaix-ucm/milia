/*
 * Copyright 2008 Sergio Pascual
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

#ifndef MILIA_EXCEPTION_H
#define MILIA_EXCEPTION_H

#include <exception>
#include <string>

namespace milia
{

    class exception : public std::exception
    {
    public:
        exception(const std::string& information);

        virtual ~exception() throw();

        const char* what() const throw();

    private:
        std::string m_message;
    };

    class recollapse : public milia::exception
    {
    public:
        recollapse(const std::string& information);

        virtual ~recollapse() throw();
            
    };

    class no_big_bang : public milia::exception
    {
    public:
        no_big_bang(const std::string& information);

        virtual ~no_big_bang() throw();
        
    };





} // namespace milia

#endif /* MILIA_EXCEPTION_H */
