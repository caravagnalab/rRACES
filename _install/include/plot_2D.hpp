/**
 * @file plot_2D.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a 2D plot window
 * @version 0.8
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

#ifndef __RACES_PLOT_2D__
#define __RACES_PLOT_2D__

#include <string>
#include <cstdint>

#include "palette.hpp"

namespace Races 
{

namespace UI 
{

/**
 * @brief Dummy windows for 2D plotting
 */
class Plot2DWindow {
protected:
	Color font_color;			//!< font color
	Color background_color;		//!< background color

	bool w_closed;				//!< flag reporting whether the window has been closed
public:
	/**
	 * @brief The constructor
	 * 
	 * @param width is the width of the plot window
	 * @param height is the height of the plot window
	 * @param name is the name of the plot window
	 */
	Plot2DWindow(const unsigned int width, const unsigned int height, const std::string& name);

	/**
	 * @brief Clear the window
	 */
	void clear();

	/**
	 * @brief Get the selected color
	 * 
	 * @return the selected color
	 */
	Color get_color() const;

	/**
	 * @brief Use a color
	 * 
	 * @param color is the color to be used
	 */
	void use_color(const Color& color);

	/**
	 * @brief Get the background color
	 * 
	 * @return the background color
	 */
	Color get_font_color() const;

	/**
	 * @brief Set the font color
	 * 
	 * @param color is the color to be set a font color
	 */
	void set_font_color(const Color& color);

	/**
	 * @brief Get the background color
	 * 
	 * @return the background color
	 */
	Color get_background_color() const;

	/**
	 * @brief Set the background color
	 * 
	 * @param color is the color to be set a background color
	 */
	void set_background_color(const Color& color);

	/**
	 * @brief Draw a point
	 * 
	 * @param x is the x-axis position of the point
	 * @param y is the y-axis position of the point 
	 */
	void draw_point(const unsigned int x, const unsigned int y);

	/**
	 * @brief Draw a rectangle
	 * 
	 * @param upper_left_x is the x-axis position of the rectangle upper left corner 
	 * @param upper_left_y is the y-axis position of the rectangle upper left corner 
	 * @param width is the rectangle width
	 * @param height is the rectangle height
	 * @param thickness is the rectangle thickness
	 */
	void draw_rectangle(const unsigned int upper_left_x, const unsigned int upper_left_y,
						const unsigned int width, const unsigned int height, const unsigned int thickness=1);

	/**
	 * @brief Draw a filled rectangle
	 * 
	 * @param upper_left_x is the x-axis position of the rectangle upper left corner 
	 * @param upper_left_y is the y-axis position of the rectangle upper left corner 
	 * @param width is the rectangle width
	 * @param height is the rectangle height
	 */
	void draw_filled_rectangle(const unsigned int upper_left_x, const unsigned int upper_left_y,
							   const unsigned int width, const unsigned int height);

	/**
	 * @brief Get the size of the a string 
	 * 
	 * @param[in] text is the string to be drawn
	 * @param[out] width is the width of the string
	 * @param[out] height is the height of the string
	 */
	void get_text_size(const std::string& text, unsigned int& width, unsigned int& height);

	/**
	 * @brief Draw a string
	 * 
	 * @param text is the string to be drawn
	 * @param upper_left_x is the x-axis position of the string upper left corner
	 * @param upper_left_y is the y-axis position of the string upper left corner
	 */
	void draw_text(const std::string& text, const unsigned int upper_left_x, const unsigned int upper_left_y);

	/**
	 * @brief Delete a point
	 * 
	 * This method delete a point in the plot by drawing over it a point 
	 * whose color is the set background color.
	 * 
	 * @param x is the x-axis position of the point
	 * @param y is the y-axis position of the point
	 */
	void delete_point(const unsigned int x, const unsigned int y);

	/**
	 * @brief Update the plot in the window
	 */
	void update();

	/**
	 * @brief Test whenever the plotting window has been closed
	 * 
	 * @return `true` if and only if the plotting window has been closed
	 */
	const bool& closed() const;

	/**
	 * @brief Test whenever the plotting window is waiting for any event
	 * 
	 * @return `false`
	 */
	bool waiting_end() const;

	/**
	 * @brief Get the number of dimensions plotted
	 * 
	 * @return the number of dimensions supported by 
	 *    this class, i.e., 2
	 */
	static constexpr uint8_t dimensions();
};

inline constexpr uint8_t Plot2DWindow::dimensions()
{
	return 2;
}

inline void Plot2DWindow::clear()
{
}
	
inline void Plot2DWindow::use_color(const Color& color)
{
	(void)color;
}

inline void Plot2DWindow::set_font_color(const Color& color)
{
	font_color = color;
}

inline void Plot2DWindow::set_background_color(const Color& color)
{
	background_color = color;
}
	
inline Color Plot2DWindow::get_color() const
{
	return Color();
}

inline Color Plot2DWindow::get_font_color() const
{
	return font_color;
}

inline Color Plot2DWindow::get_background_color() const
{
	return background_color;
}

inline void Plot2DWindow::draw_point(const unsigned int x, const unsigned int y)
{
	(void)x;
	(void)y;
}

inline void Plot2DWindow::draw_rectangle(const unsigned int upper_left_x, const unsigned int upper_left_y,
						   				 const unsigned int width, const unsigned int height, const unsigned int thickness)
{
	(void)upper_left_x;
	(void)upper_left_y;
	(void)width;
	(void)height;
	(void)thickness;
}

inline void Plot2DWindow::draw_filled_rectangle(const unsigned int upper_left_x, const unsigned int upper_left_y,
								  				const unsigned int width, const unsigned int height)
{
	(void)upper_left_x;
	(void)upper_left_y;
	(void)width;
	(void)height;
}

inline void Plot2DWindow::get_text_size(const std::string& text, unsigned int& width, unsigned int& height)
{
	(void)text;
	(void)width;
	(void)height;
}

inline void Plot2DWindow::draw_text(const std::string& text, const unsigned int upper_left_x, const unsigned int upper_left_y)
{
	(void)upper_left_x;
	(void)upper_left_y;
	(void)text;
}

inline void Plot2DWindow::delete_point(const unsigned int x, const unsigned int y)
{
	(void)x;
	(void)y;
}

inline void Plot2DWindow::update()
{
}

inline const bool& Plot2DWindow::closed() const
{
	return w_closed;
}

inline bool Plot2DWindow::waiting_end() const
{
	return false;
}

}  // UI

}  // Races

#endif // __RACES_PLOT_2D__