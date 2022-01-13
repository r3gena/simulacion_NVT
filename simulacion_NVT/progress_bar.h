#pragma once
#include <iostream>
#include <cmath>

class progress_bar
{
public:
	int max_steps;
	int progress;
	int previous_progress;
	std::string process_name;

	progress_bar(std::string name, int max_steps);
	~progress_bar();

	void update(int current_step);
};
