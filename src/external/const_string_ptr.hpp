/**
 *  This file is part of pgrep.
 *
 *  pgrep is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  pgrep is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright 2015-2016  Sarah Lutteropp, Bernhard Scheirle and Kevin Zerr
 */

// Source file adapted from https://code.ipd.kit.edu/pmp/pgrep/tree/development/lib/src/util/ConstrStringPtr.h
// see also https://ps.ipd.kit.edu/176_3088.php

#pragma once

#include <string>
#include <stdexcept>
#include <cassert>

namespace seq_correct {

namespace external {

class ConstStringPtr
{
private:
	const char* m_ptr;
	std::size_t m_length;
	const std::string* m_originalString;

public:
	//constructs a ConstStringPtr from other (does not copy)
	//ConstStringPtr(const std::string& other);
	//constructs a ConstStringPtr from other (does not copy)
	ConstStringPtr(const std::string* other)
		: m_ptr(other->data())
		, m_length(other->length())
		, m_originalString(other)
	{}
	//constructs a ConstStringPtr
	ConstStringPtr(const char* _ptr, std::size_t _length, const std::string* _originalString)
		: m_ptr(_ptr)
		, m_length(_length)
		, m_originalString(_originalString)
	{};
	ConstStringPtr(const ConstStringPtr& other)
		: m_ptr(other.m_ptr)
		, m_length(other.m_length)
		, m_originalString(other.m_originalString)
	{};
	ConstStringPtr(ConstStringPtr&& other)
		: m_ptr(other.m_ptr)
		, m_length(other.m_length)
		, m_originalString(other.m_originalString)
	{};

	ConstStringPtr& operator=(const ConstStringPtr& other)
	{
		m_ptr = other.m_ptr;
		m_length = other.m_length;
		m_originalString = other.m_originalString;
		return *this;
	};
	ConstStringPtr& operator=(ConstStringPtr&& other)
	{
		m_ptr = other.m_ptr;
		m_length = other.m_length;
		m_originalString = other.m_originalString;
		return *this;
	};

	//access specified character with bounds checking
	const char& at(std::size_t position) const
	{
		if(position >= m_length)
		{
			std::string error("ConstStringPtr::at: given position (");
			error += std::to_string(position);
			error += ") is out of range (";
			error += std::to_string(m_length);
			error += ").";
			throw std::out_of_range(error);
		}
		return m_ptr[position];
	};
	//access specified character
	const char& operator[](std::size_t position) const
	{
		assert(position <= m_length);	// allow access to one byte after the text (some StringMatcher requires that..)
		return m_ptr[position];
	}

	//access specified character, return value from 0..255
	unsigned char getU(std::size_t index) const {
		return static_cast<unsigned char>(m_ptr[index]);
	}

	//returns the number of characters
	std::size_t size() const
	{
		return m_length;
	};
	//returns the number of characters
	std::size_t length() const
	{
		return m_length;
	};
	// return if size == 0
	bool empty() const
	{
		return m_length == 0;
	};

	//returns a pointer to the first character of a string
	const char* data() const
	{
		return m_ptr;
	};

	/// returns a unsigned char pointer to the first character of the string
	const unsigned char* dataU() const {
		return reinterpret_cast<const unsigned char*>(m_ptr);
	}

	/// returns a unsigned char pointer to the character defined by the offset
	const unsigned char* dataU(std::size_t offset) const {
		assert(offset <= m_length);

		return dataU() + offset;
	}

	// returns a substring from offset until the end (does not copy the string)
	ConstStringPtr substr(std::size_t startOffset) const {
		assert(startOffset <= m_length);

		return ConstStringPtr(this->m_ptr + startOffset, m_length - startOffset, m_originalString);
	}

	//returns a substring (does not copy the string)
	ConstStringPtr substr(std::size_t startOffset, std::size_t newLength) const
	{
		if (startOffset + newLength > this->m_length)
		{
			std::string error("ConstStringPtr::substr given startOffset (");
			error += std::to_string(startOffset);
			error += ") and newLength (";
			error += std::to_string(newLength);
			error += ") are greater equal than the max length (";
			error += std::to_string(this->m_length);
			error += ").";
			throw std::out_of_range(error);
		}
		return ConstStringPtr(this->m_ptr + startOffset, newLength, m_originalString);
	};

	//returns a copy of this string as std::string
	std::string toString() const
	{
		return std::string(m_ptr, m_length);
	};

	//returns the original string
	const std::string* originalString() const
	{
		return m_originalString;
	};

	//returns the offset in the originalString
	std::size_t originalStringOffset() const {
		// assert that the orginal string is in the range
		assert(m_ptr >= m_originalString->data());
		assert(m_ptr <= m_originalString->data() + m_originalString->size());

		return static_cast<std::size_t>(m_ptr - m_originalString->data());
	}

	// returns the offset of the first occurrence of the char, or std::string::npos if there is no match
	std::size_t find(char c, std::size_t pos = 0) const {
		assert(pos <= m_length);

		std::size_t result = m_originalString->find(c, originalStringOffset() + pos);

		if (result == std::string::npos)
			return std::string::npos;

		assert(result >= pos);

		return result - originalStringOffset();
	}
};

} // end of namespace seq_correct::external
} // end of namespace seq_correct
