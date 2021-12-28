#include "progress_bar.h"

progress_bar::progress_bar(std::string name, int max_steps)
{
	this->max_steps = max_steps;
	progress = 0;
	previous_progress = 0;
	process_name = name;
	std::cout << process_name << std::endl;
	std::cout << "[" << std::string(100,'-') << ']' << "\r";
	previous_progress = progress;
}

progress_bar::~progress_bar()
{
	std::cout << "[" << std::string(100, '#') << ']' << "\r";
	std::cout.flush();
	std::cout << "\n" << "\n";
}

void progress_bar::update(int current_step)
{
	progress = std::round(100 * current_step / max_steps);
	if (progress != previous_progress)
	{
		std::cout << '[' << std::string(progress, '#') << std::string(100 - progress, '-') << ']' << "\r";
		std::cout.flush();
	}
}
