/**
 * @file palette.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines colors and a palette
 * @version 0.3
 * @date 2023-10-02
 * 
 * @copyright Copyright (c) 2023
 * 
 * MIT License
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef __RACES_PALETTE__
#define __RACES_PALETTE__

#include <vector>
#include <cstdint>

namespace Races 
{

/**
 * @brief The namespace of RACES UI
 */
namespace UI 
{

/**
 * @brief A class to represent colors
 */
struct Color {
	uint8_t red;	//!< The red level
	uint8_t green;	//!< The green level
	uint8_t blue;	//!< The blue level
	uint8_t alpha;	//!< The alpha level

	/**
	 * @brief The empty constructor
	 */
	Color();

	/**
	 * @brief A constructor
	 * 
	 * This method creates a totally opaque color.
	 * 
	 * @param red is the red level of the new color
	 * @param green is the green level of the new color
	 * @param blue is the blue level of the new color
	 */
	Color(const uint8_t red, const uint8_t green, const uint8_t blue);

	/**
	 * @brief A constructor
	 * 
	 * @param red is the red level of the new color
	 * @param green is the green level of the new color
	 * @param blue is the blue level of the new color
	 * @param alpha is the alpha level of the new color
	 */
	Color(const uint8_t red, const uint8_t green, const uint8_t blue, const uint8_t alpha);
};

/**
 * @brief The UI palette
 */
extern std::vector<Color> palette;

} // UI

} // Races
#endif // __RACES_PALETTE__