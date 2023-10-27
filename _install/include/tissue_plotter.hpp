/**
 * @file tissue_plotter.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a UI window to plot a tissue
 * @version 0.12
 * @date 2023-10-12
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

#ifndef __RACES_TISSUE_PLOTTER__
#define __RACES_TISSUE_PLOTTER__

#include <memory> // smart pointers
#include <sstream>
#include <chrono>

#include "statistics.hpp"

#include "tissue.hpp"
#include "palette.hpp"

namespace Races 
{

namespace UI 
{

template <class Rep, std::intmax_t num, std::intmax_t denom>
auto extract_time_units(std::chrono::duration<Rep, std::ratio<num, denom>> duration)
{
    using namespace std::chrono;

    const auto hrs = duration_cast<hours>(duration);
    duration -= hrs;

    const auto mins = duration_cast<minutes>(duration);
    duration -= mins;

    const auto secs = duration_cast<seconds>(duration);

    return std::make_tuple(hrs, mins, secs);
}

template <class Rep, std::intmax_t num, std::intmax_t denom>
std::string format_duration(const std::chrono::duration<Rep, std::ratio<num, denom>> duration)
{
    const auto [hrs, mins, secs] = extract_time_units(duration);

    std::ostringstream oss;

    oss << std::setfill('0') << std::setw(2) << hrs.count() << ":" 
		<< std::setfill('0') << std::setw(2) << mins.count() << ":" 
		<< std::setfill('0') << std::setw(2) << secs.count();

    return oss.str();
}

/**
 * @brief A tissue plotting class
 * 
 * @tparam PLOT_WINDOW is the type of plotting window
 */
template<typename PLOT_WINDOW>
class TissuePlotter {
	using Tissue = Drivers::Simulation::Tissue;

	const Tissue& tissue;					//!< the tissue to plot

	std::unique_ptr<PLOT_WINDOW> window;	//!< the plotting window

	const unsigned int frame_border;		//!< frame border width
	const unsigned int frame_thickness;		//!< frame line thickness

	const unsigned int legend_rectangle_width;	//!< legend rectangle width
	const unsigned int legend_rectangle_height;	//!< legend rectangle height
	const unsigned int legend_label_height;		//!< legend label height

	std::chrono::milliseconds redraw_interval;	//!< milliseconds between two redraws
	std::chrono::system_clock::time_point last_redraw_time;    		//!< the last redraw time in milliseconds

	/**
	 * @brief Get the total width of the plotting window
	 * 
	 * @return the total width of the plotting window in pixels
	 */
	unsigned int total_width() const
	{
		return static_cast<unsigned int>(tissue.size()[0])
				+3*frame_border+legend_rectangle_width;
	}

	/**
	 * @brief Get the total height of the plotting window
	 * 
	 * @return the total height of the plotting window in pixels
	 */
	unsigned int total_height() const
	{
		unsigned int height =  tissue.num_of_species()*(legend_rectangle_height+frame_border);

		return std::max(height, static_cast<unsigned int>(tissue.size()[1]+2*frame_border+legend_rectangle_height));
	}

	/**
	 * @brief Draw the tissue frame
	 */
	void draw_tissue_frame()
	{
		auto sizes(tissue.size());
		if constexpr(PLOT_WINDOW::dimensions()==2) {
			window->use_color(window->get_font_color());
			window->draw_rectangle(frame_border-frame_thickness, 
								   legend_rectangle_height+frame_border-frame_thickness,
							       sizes[0]+2*frame_thickness, sizes[1]+2*frame_thickness);
		}
	}

	/**
	 * @brief Draw the simulated time
	 * 
	 * @param statistics are the simulation statistics
	 */
	void draw_time(const Drivers::Simulation::TissueStatistics& statistics)
	{
		unsigned int label_width, label_height;

		std::ostringstream oss;
		oss << "Time: " << statistics.get_simulated_time() << "s";

		window->get_text_size(oss.str(), label_width, label_height);
		const unsigned int label_x_pos = frame_border;
		const unsigned int label_y_pos = (frame_border+legend_rectangle_height-label_height)/2;

		window->use_color(window->get_font_color());
		window->draw_text(oss.str(), label_x_pos, label_y_pos);
	}

	/**
	 * @brief Draw the elapsed real time
	 * 
	 * @param statistics are the simulation statistics
	 */
	void draw_elapsed_time(const Drivers::Simulation::TissueStatistics& statistics)
	{
		unsigned int label_width, label_height;

		std::ostringstream oss;
		oss << "Elapsed: " << format_duration(statistics.get_elapsed_time());

		window->get_text_size(oss.str(), label_width, label_height);
		const unsigned int label_x_pos=frame_border+300;
		const unsigned int label_y_pos = (frame_border+legend_rectangle_height-label_height)/2;

		window->use_color(window->get_font_color());
		window->draw_text(oss.str(), label_x_pos, label_y_pos);
	}

	/**
	 * @brief Draw the instantaneous number of events per second
	 * 
	 * @param statistics are the simulation statistics
	 */
	void draw_events_per_second(const Drivers::Simulation::TissueStatistics& statistics)
	{
		using namespace std::chrono;
		unsigned int label_width, label_height;

		std::ostringstream oss;
		oss << "Events/ms: "
			<< std::setprecision(1)
			<< std::fixed
		    << statistics.get_last_recorded_events_over_time<milliseconds>();

		window->get_text_size(oss.str(), label_width, label_height);
		const unsigned int label_x_pos=frame_border+600;
		const unsigned int label_y_pos = (frame_border+legend_rectangle_height-label_height)/2;

		window->use_color(window->get_font_color());
		window->draw_text(oss.str(), label_x_pos, label_y_pos);
	}

	void draw_time_statistics(const Drivers::Simulation::TissueStatistics& statistics)
	{
		draw_time(statistics);
		draw_elapsed_time(statistics);
		draw_events_per_second(statistics);
	}

	/**
	 * @brief Draw the tissue
	 */
	void draw_tissue()
	{	
		draw_tissue_frame();

		for (const auto& species: tissue) {
			const auto& color = palette[species.get_id()];

			window->use_color(color);

			for (const auto& cell: species) {
				if constexpr(PLOT_WINDOW::dimensions()==2) {
					window->draw_point(cell.x+frame_border,
									   cell.y+frame_border+legend_rectangle_height);
				} else {
					window->draw_point(cell.x+frame_border,
									   cell.y+frame_border,
									   cell.z+frame_border);
				}
			}
		}
	}

	void draw_species_legend(const Drivers::Simulation::Species& species,
							 const Drivers::Simulation::SpeciesStatistics& statistics, 
							 const unsigned int& x, const unsigned int& y, const Color& color)
	{
		window->use_color(color);
		window->draw_filled_rectangle(x,y,legend_rectangle_width,legend_rectangle_height);

		window->use_color(window->get_font_color());

		std::ostringstream oss;

		oss << species.get_name() 
			<< ": " << species.num_of_cells() << "/" << statistics.total_cells
			<< "/" << statistics.killed_cells << "/" << statistics.lost_cells;

		unsigned int label_width, label_height;

		window->get_text_size(oss.str(), label_width, label_height);
		const unsigned int label_x_pos = x+frame_border; //+(legend_rectangle_width-label_height)/2;
		const unsigned int label_y_pos = y+(legend_rectangle_height-label_height)/2;
		window->draw_text(oss.str(), label_x_pos, label_y_pos);
	}

	/**
	 * @brief Draw the legend
	 */
	void draw_legend(const Drivers::Simulation::TissueStatistics& statistics)
	{
		auto sizes(tissue.size());
		const unsigned int x = 2*frame_border + static_cast<unsigned int>(sizes[0]);
		unsigned int y = legend_rectangle_height+frame_border-frame_thickness;

		for (const auto& species: tissue) {
			if (statistics.contains_data_for(species.get_id())) {
				draw_species_legend(species, statistics.at(species.get_id()),
									x, y, palette[species.get_id()]);
			} else {
				Drivers::Simulation::SpeciesStatistics s_statistics;
				draw_species_legend(species, s_statistics,
									x, y, palette[species.get_id()]);
			}

			y += legend_rectangle_height+frame_border;
		}
	}

public:
	/**
	 * @brief A constructor
	 * 
	 * @param tissue is the tissue to plot
	 * @param frames_per_second is the number of frames per second
	 */
	explicit TissuePlotter(const Drivers::Simulation::Tissue& tissue, const unsigned int frames_per_second=10):
		TissuePlotter<PLOT_WINDOW>(tissue, "RACES Simulation"+((tissue.get_name()=="")?"":" - "+tissue.get_name()), frames_per_second)
	{
	}

	/**
	 * @brief A constructor
	 * 
	 * @param tissue is the tissue to plot
	 * @param name is the name of the plotting window
	 * @param frames_per_second is the number of frames per second
	 */
	TissuePlotter(const Drivers::Simulation::Tissue& tissue, const std::string& name,
				  const unsigned int frames_per_second=10):
		tissue(tissue), frame_border(20), frame_thickness(5), legend_rectangle_width(400),
		legend_rectangle_height(35), legend_label_height(16), redraw_interval(1000/frames_per_second), 
		last_redraw_time(std::chrono::system_clock::from_time_t(0))
	{
		if (tissue.num_of_species()>palette.size()) {
			throw std::domain_error("The color palette does not support so many species");
		}

		auto sizes(tissue.size());
		if constexpr(PLOT_WINDOW::dimensions()==2) {
			window = std::make_unique<PLOT_WINDOW>(total_width(), total_height(), name);
			
			return;
		}

		throw std::runtime_error("Unsupported plot window type");
	}

	/**
	 * @brief Plot the tissue
	 * 
	 * @param statistics are the simulation statistics
	 */
    void plot(const Drivers::Simulation::TissueStatistics& statistics)
	{
		using namespace std::chrono;

		auto time_delta = system_clock::now()-last_redraw_time;
		if (duration_cast< milliseconds >(time_delta) >=redraw_interval) {
			last_redraw_time = std::chrono::system_clock::now();
			window->clear();
		
			draw_time_statistics(statistics);

			draw_legend(statistics);
			draw_tissue();
			
			window->update();
		}
	}

	/**
	 * @brief Test whether the plotting window has been closed
	 * 
	 * @return `true` if and only if the plotting window has been closed
	 */
	inline const bool& closed() const
	{
		return window->closed();
	}

	/**
	 * @brief Test whenever the plotting window is waiting for any event
	 * 
	 * @return `true` if and only if the plotting window is waiting 
	 *     for any event
	 */
	inline bool waiting_end() const
	{
		return window->waiting_end();
	}

	/**
	 * @brief Set the number of frame per second
	 * 
	 * @param frames_per_second is the number of frame per seconds to 
	 * 		be plot at most
	 */
	inline void set_frames_per_second(const unsigned int frames_per_second)
	{
		redraw_interval = std::chrono::milliseconds(1000/frames_per_second);
	}

	/**
	 * @brief Set the font color
	 * 
	 * @param font_color is the new font color
	 */
	inline void set_font_color(const Color& font_color)
	{
		window->set_font_color(font_color);
	}

	/**
	 * @brief Set the background color
	 * 
	 * @param background_color is the new background color
	 */
	inline void set_background_color(const Color& background_color)
	{
		window->set_background_color(background_color);
	}

	/**
	 * @brief The destroyer
	 */
	~TissuePlotter()
	{}
};

}	// UI

}	// Races

#endif // __RACES_TISSUE_PLOTTER__