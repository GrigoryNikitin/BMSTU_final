#include "pch.h"
#include "Date_time_ms.h"


Date_time_ms::Date_time_ms()
{
	Date = { 0 };
	ms = 0;
}


Date_time_ms::~Date_time_ms()
{
}

void Date_time_ms::set_date_time(std::string a, int ms)
{
	std::istringstream iss(a);
	iss >> std::get_time(&Date, "%Y-%b-%d %H:%M:%S");
	this->ms = ms;
	update_date_time(); // обновляю на всяк случай сразу
}

std::string Date_time_ms::get_string()
{
	char Date_char[32];
	strftime(Date_char, 32, "%Y-%b-%d %H:%M:%S", &Date); // из tm в массив из char
	std::string Date_string = Date_char;
	/*Date_char[20] = '.';
	Date_char[21] = '0' + (int)(ms / 100);
	Date_char[22] = '0' + (int)((ms % 100) / 10);
	Date_char[23] = '0' + (int)(ms % 10);*/
	Date_string += '.';
	Date_string += std::to_string(ms);
	return Date_string;
}

void Date_time_ms::update_date_time()
{
	mktime(&Date);
	if (ms >= 1000) // если миллисекунды больше или равны 1000
	{
		Date.tm_sec += 1;
		mktime(&Date);
		ms -= 1000;
	}
	if (ms < 0) // если миллисекунды меньше 0
	{
		do
		{
			Date.tm_sec -= 1;
			mktime(&Date);
			ms += 1000;
		} while (ms < 0);
	}
}

Date_time_ms& Date_time_ms::operator=(const Date_time_ms& right)
{
	//проверка на самоприсваивание
	if (this == &right) {
		return *this;
	}
	Date = right.Date;
	ms = right.ms;
}