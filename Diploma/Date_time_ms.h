#pragma once
class Date_time_ms
{
public:
	Date_time_ms();
	~Date_time_ms();
	void set_date_time(std::string a, int ms);
	std::string get_string();
	void update_date_time();
	Date_time_ms& operator=(const Date_time_ms& right);

	tm Date;
	int ms;
};