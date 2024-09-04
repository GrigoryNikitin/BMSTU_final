// Course_project.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <iostream>
#include "Spacecraft.h"

#pragma comment(lib, "../libs/cspice/lib/cspice.lib")
extern "C" {
#include "../libs/cspice/include/SpiceCK.h" // must be before SpiceZdf.h and SpiceZpr.h
#include "../libs/cspice/include/SpiceZdf.h"
#include "../libs/cspice/include/SpiceZpr.h"
}

std::vector <Spacecraft> Integration(Spacecraft Sputnik, double integration_time, double h)
{
	std::vector <Spacecraft> SC_Condition;
	SC_Condition.push_back(Sputnik);
	while (Sputnik.time <= integration_time - h)
	{
		Sputnik.Runge_Kutta(h);
		SC_Condition.push_back(Sputnik);
		//if (Sputnik.distance_to_Moon() < 1737) { std::cout << "Fell on the Moon" << '\n'; }
		//if (Sputnik.distance_to_Earth() < 6371) { std::cout << "Fell on the Earth" << '\n'; }
	}
	return SC_Condition;
}

void Integration_void(Spacecraft Sputnik, double integration_time, double h, int b)
{
	std::ofstream fout1("Integ.txt");
	std::ofstream fout2("Kepler.txt");
	fout1 << "Date" << '\t' << "Time" << '\t' << "X" << '\t' << "Y" << '\t' << "Z" << '\t' << "r_KA" << '\t' << "Расстояние до Луны" << '\t' << "Vx" << '\t' << "Vy" << '\t' << "Vz" << '\t' << "V_KA" << '\t' << "Главный центр" << '\n';
	fout2 << "Date" << '\t' << "Time" << '\t' << "a" << '\t' << "ecc" << '\t' << "inc" << '\t' << "RAAN" << '\t' << "omega" << '\t' << "tau" << '\t' << "theta" << '\t' << "Главный центр" << '\n';
	while (Sputnik.time <= integration_time - h)
	{
		int int_time = (int)trunc(Sputnik.time);
		if ((int_time % b == 0) or (Sputnik.time == integration_time - h))
		{
			Sputnik.printParameters(fout1);
			Sputnik.printKepler(fout2);
		}
		Sputnik.Runge_Kutta(h);
	}
}

void set_parameters_Moon_flight(Spacecraft& Sputnik, double flight_time, Date_time_ms& approach_date)
{
	double integer_delta_t_depart, fractional_delta_t_depart;
	Sputnik.Moon_flight_without_perturbations_north(flight_time, approach_date);
	Sputnik.a = Sputnik.a0;
	Sputnik.ecc = Sputnik.ecc0;
	Sputnik.inc = Sputnik.inc0;
	Sputnik.RAAN = Sputnik.RAAN0;
	Sputnik.omega = 0;
	Sputnik.theta = 0;
	Sputnik.tau = 0;
	Sputnik.Period = 2 * Pi * sqrt(Sputnik.a * Sputnik.a * Sputnik.a / mu_Earth);
	Sputnik.from_Kepl_to_Dec();
	Sputnik.from_Dec_to_Kepl();
	// ниже считаю, сколько лететь от u=theta=0 на круговой орбите до точки подачи импульса (u_depart я знаю) и вычитаю это время из tau1_date
	// и устанавливаю полученное как время старта
	Date_time_ms tm_Start_Date; // сюда запишу время старта
	tm_Start_Date = Sputnik.tau1_date;
	fractional_delta_t_depart = modf(Sputnik.delta_t_depart, &integer_delta_t_depart);
	tm_Start_Date.Date.tm_sec -= integer_delta_t_depart;
	tm_Start_Date.update_date_time();
	Sputnik.set_Date(tm_Start_Date.get_string().c_str(), -1000 * fractional_delta_t_depart);
	Sputnik.Date_time.update_date_time();
}

std::vector <Spacecraft> Integration_to_Moon(Spacecraft Sputnik, double flight_time, double h_near_earth, double h_outer_space, double h_near_moon)
{
	// н.у. у Спутника уже должны быть прописаны нужные
	std::vector <Spacecraft> SC_Condition;
	SC_Condition.push_back(Sputnik);

	while (Sputnik.time <= (Sputnik.delta_t_depart - h_near_earth))
	{
		Sputnik.Runge_Kutta(h_near_earth);
		SC_Condition.push_back(Sputnik);
		if (Sputnik.distance_to_Earth() < 6371) { std::cout << "Fell on the Earth" << '\n'; break; }
		// потом сделать, чтобы шаг менялся так, чтобы точно попадать в момент подачи импульса
	}
	Sputnik.Runge_Kutta(Sputnik.delta_t_depart - Sputnik.time);
	SC_Condition.push_back(Sputnik);
	Sputnik.Collinear_impulse(Sputnik.fly_away_impulse);
	Sputnik.from_Dec_to_Kepl();
	SC_Condition.push_back(Sputnik);
	while (Sputnik.time <= (flight_time * 24 * 60 * 60 + Sputnik.delta_t_depart - h_outer_space) and ((Sputnik.distance_to_Moon() > 66000)))
	{
		Sputnik.Runge_Kutta(h_outer_space);
		SC_Condition.push_back(Sputnik);
		if (Sputnik.distance_to_Earth() < 6371) { std::cout << "Fell on the Earth" << '\n'; break; }
	}
	if (Sputnik.time > (flight_time * 24 * 60 * 60 + Sputnik.delta_t_depart - h_outer_space)) // на случай, если не в СД Луны, а время закончилось
	{
		Sputnik.Runge_Kutta(flight_time * 24 * 60 * 60 + Sputnik.delta_t_depart - Sputnik.time);
		SC_Condition.push_back(Sputnik);
		return SC_Condition;
	}
	while (Sputnik.time <= (flight_time * 24 * 60 * 60 + Sputnik.delta_t_depart - h_near_moon))
	{
		if (Sputnik.distance_to_Moon() <= 5000)
		{
			Sputnik.Runge_Kutta(h_near_moon / 10);
			SC_Condition.push_back(Sputnik);
		}
		else
		{
			Sputnik.Runge_Kutta(h_near_moon);
			SC_Condition.push_back(Sputnik);
		}
	}
	Sputnik.Runge_Kutta(flight_time * 24 * 60 * 60 + Sputnik.delta_t_depart - Sputnik.time);
	SC_Condition.push_back(Sputnik);
	return SC_Condition;
}

std::vector <Spacecraft> Integration_to_Moon_periselene(Spacecraft Sputnik, double h_near_earth, double h_outer_space, double h_near_moon)
{
	// н.у. у Спутника уже должны быть прописаны нужные
	std::vector <Spacecraft> SC_Condition;
	SC_Condition.push_back(Sputnik);

	while (Sputnik.time <= (Sputnik.delta_t_depart - h_near_earth))
	{
		Sputnik.Runge_Kutta(h_near_earth);
		SC_Condition.push_back(Sputnik);
		if (Sputnik.distance_to_Earth() < 6371) { std::cout << "Fell on the Earth" << '\n'; break; }
		// потом сделать, чтобы шаг менялся так, чтобы точно попадать в момент подачи импульса
	}
	Sputnik.Runge_Kutta(Sputnik.delta_t_depart - Sputnik.time);
	SC_Condition.push_back(Sputnik);
	Sputnik.Collinear_impulse(Sputnik.fly_away_impulse);
	Sputnik.from_Dec_to_Kepl();
	SC_Condition.push_back(Sputnik);
	double distance_to_Moon_last = Sputnik.distance_to_Moon();
	while ((Sputnik.distance_to_Moon() <= distance_to_Moon_last) and (Sputnik.distance_to_Moon() > 1737.1) /*or (Sputnik.main_center_flag == "Earth")*/)
	{
		if (Sputnik.distance_to_Moon() > 66000)
		{
			distance_to_Moon_last = Sputnik.distance_to_Moon();
			Sputnik.Runge_Kutta(h_outer_space);
			SC_Condition.push_back(Sputnik);
			if (Sputnik.distance_to_Earth() < 6371) { std::cout << "Fell on the Earth" << '\n'; break; }
		}
		if (Sputnik.distance_to_Moon() <= 66000)
		{
			distance_to_Moon_last = Sputnik.distance_to_Moon();
			Sputnik.Runge_Kutta(h_near_moon);
			SC_Condition.push_back(Sputnik);
			//if (Sputnik.distance_to_Moon() < 1737.1) { std::cout << "Fell on the Moon" << '\n'; break; } // лечу только до пов-ти
		}
	}
	if (Sputnik.distance_to_Moon() < 1737.1) // если надо только до пов-ти
	{
		SC_Condition.pop_back(); //удалил последний элемент массива, так как там спутник уже дальше периселения ушёл (либо внутрь Луны)
		double h_new = h_near_moon;
		Sputnik.Runge_Kutta(-2 * h_new);
		double distance_to_Moon_next;
		h_new /= 2;
		do
		{
			do
			{
				distance_to_Moon_last = Sputnik.distance_to_Moon();
				Sputnik.Runge_Kutta(h_new);
				distance_to_Moon_next = Sputnik.distance_to_Moon();
			} while ((distance_to_Moon_next <= distance_to_Moon_last) and (distance_to_Moon_next > 1737.1));
			Sputnik.Runge_Kutta(-2 * h_new);
			h_new /= 4;
		} while (abs(distance_to_Moon_last - distance_to_Moon_next) > pow(10, -5));
		SC_Condition.push_back(Sputnik);
		return SC_Condition;
	}
	else
	{
		SC_Condition.pop_back(); //удалил последний элемент массива, так как там спутник уже дальше периселения ушёл (либо внутрь Луны)
		double h_new;
		if (Sputnik.distance_to_Moon() > 66000)
		{
			h_new = h_outer_space;
		}
		if (Sputnik.distance_to_Moon() <= 66000)
		{
			h_new = h_near_moon;
		}
		Sputnik.Runge_Kutta(-2 * h_new);
		double distance_to_Moon_next;
		h_new /= 2;
		do
		{
			do
			{
				distance_to_Moon_last = Sputnik.distance_to_Moon();
				Sputnik.Runge_Kutta(h_new);
				distance_to_Moon_next = Sputnik.distance_to_Moon();
			} while ((distance_to_Moon_next < distance_to_Moon_last));
			Sputnik.Runge_Kutta(-2 * h_new);
			h_new /= 4;
		} while (abs(distance_to_Moon_last - distance_to_Moon_next) > pow(10, -5) and (h_new > 0.1));
		SC_Condition.push_back(Sputnik);
		return SC_Condition;
	}
}

std::vector <Spacecraft> Integration_to_Moon_aposelene(Spacecraft Sputnik, double h_near_moon)
{
	// н.у. у Спутника уже должны быть прописаны нужные
	std::vector <Spacecraft> SC_Condition;
	SC_Condition.push_back(Sputnik);

	double distance_to_Moon_last = Sputnik.distance_to_Moon();
	while ((Sputnik.distance_to_Moon() >= distance_to_Moon_last))
	{
		distance_to_Moon_last = Sputnik.distance_to_Moon();
		Sputnik.Runge_Kutta(h_near_moon);
		SC_Condition.push_back(Sputnik);
	}
	SC_Condition.pop_back(); //удалил последний элемент массива, так как там спутник уже дальше апоселения ушёл
	double h_new;
	h_new = h_near_moon;
	Sputnik.Runge_Kutta(-2 * h_new);
	double distance_to_Moon_next;
	h_new /= 2;
	do
	{
		do
		{
			distance_to_Moon_last = Sputnik.distance_to_Moon();
			Sputnik.Runge_Kutta(h_new);
			distance_to_Moon_next = Sputnik.distance_to_Moon();
		} while ((distance_to_Moon_next > distance_to_Moon_last));
		Sputnik.Runge_Kutta(-2 * h_new);
		h_new /= 4;
	} while (abs(distance_to_Moon_last - distance_to_Moon_next) > pow(10, -5) and (h_new > 0.1));
	SC_Condition.push_back(Sputnik);
	return SC_Condition;
}

std::vector <Spacecraft> Integration_full_period(Spacecraft Sputnik, double h_near_moon)
{
	// н.у. у Спутника уже должны быть прописаны нужные
	std::vector <Spacecraft> SC_Condition;
	SC_Condition.push_back(Sputnik);

	Sputnik.Runge_Kutta(h_near_moon);
	SC_Condition.push_back(Sputnik);

	while ((Sputnik.theta <= 357 * Pi / 180))
	{
		Sputnik.Runge_Kutta(h_near_moon);
		SC_Condition.push_back(Sputnik);
	}
	double h_new;
	h_new = h_near_moon;
	h_new /= 2;
	do
	{
		do
		{
			Sputnik.Runge_Kutta(h_new);
		} while ((Sputnik.theta > 357 * Pi / 180));
		Sputnik.Runge_Kutta(-2 * h_new);
		h_new /= 4;
	} while (abs(Sputnik.theta) > pow(10, -5) and (h_new > 0.1));
	SC_Condition.push_back(Sputnik);
	return SC_Condition;
}

void Newtons_method_coords(Spacecraft& Sputnik, double flight_time, Date_time_ms& approach_date, double h_near_earth, double h_outer_space, double h_near_moon)
{
	Spacecraft Nominal_Sputnik, RAAN_Sputnik, tau_Sputnik, impulse_Sputnik, Sputnik_last;
	std::vector <Spacecraft> Nominal_trajectory; //номинальная траектория (и другие)
	set_parameters_Moon_flight(Sputnik, flight_time, approach_date); // устанавливаю Спутнику первые приближения

	std::vector <double> Nevyazki(3), Nevyazki_last(3), Variations(3);
	matrix Jacoby(3, 3); // матрица Якоби
	double delta_X_RAAN_plus, delta_Y_RAAN_plus, delta_Z_RAAN_plus,
		delta_X_RAAN_minus, delta_Y_RAAN_minus, delta_Z_RAAN_minus;
	double delta_X_tau_plus, delta_Y_tau_plus, delta_Z_tau_plus,
		delta_X_tau_minus, delta_Y_tau_minus, delta_Z_tau_minus;
	double delta_X_impulse_plus, delta_Y_impulse_plus, delta_Z_impulse_plus,
		delta_X_impulse_minus, delta_Y_impulse_minus, delta_Z_impulse_minus;

	double X_target = 0, Y_target = 0, Z_target = 0; // задаю условия, в которые нужно попасть
	double delta_RAAN = 0.00000001, delta_tau = 0.01 /*1 секунда*/, delta_impulse = 0.0000001; // шаги для производных
	double delta_R, delta_R_last;
	double koeff = 1; // коэффициент для прибавления приращений
	int k = 0; // номер итерации

	std::ofstream fout("Newton_coords.txt");
	fout << "Номер итерации" << '\t' << "Дельта X" << '\t' << "Дельта Y" << '\t' << "Дельта Z" << '\t' << "Дельта R" << '\t' << "Коэффициент K" << '\n';
	Nominal_Sputnik = Sputnik;
	Nominal_trajectory = Integration_to_Moon(Nominal_Sputnik, flight_time, h_near_earth, h_outer_space, h_near_moon);
	if (Nominal_trajectory.back().main_center_flag == "Moon")
	{
		Nevyazki[0] = X_target - Nominal_trajectory.back().get_parameters()[0]; //посчитал невязки
		Nevyazki[1] = Y_target - Nominal_trajectory.back().get_parameters()[1];
		Nevyazki[2] = Z_target - Nominal_trajectory.back().get_parameters()[2];
	}
	if (Nominal_trajectory.back().main_center_flag == "Earth")
	{
		Date_time_ms Date_time = Nominal_trajectory.back().Date_time;
		std::vector <double> Moon_state(6); // вектор параметров Луны
		double X_M, Y_M, Z_M;
		double X_M0, Y_M0, Z_M0;

		Moon_state = Nominal_trajectory.back().get_planet_parameters(Date_time.get_string().c_str(), "MOON", "J2000", "EARTH");
		X_M = Moon_state[0]; // вектор Земля - Луна
		Y_M = Moon_state[1];
		Z_M = Moon_state[2];
		Nevyazki[0] = X_target - (Nominal_trajectory.back().get_parameters()[0] - X_M);
		Nevyazki[1] = Y_target - (Nominal_trajectory.back().get_parameters()[1] - Y_M);
		Nevyazki[2] = Z_target - (Nominal_trajectory.back().get_parameters()[2] - Z_M);
	}
	delta_R = sqrt(Nevyazki[0] * Nevyazki[0] + Nevyazki[1] * Nevyazki[1] + Nevyazki[2] * Nevyazki[2]);

	fout << k << '\t' << Nevyazki[0] << '\t' << Nevyazki[1] << '\t' << Nevyazki[2] << '\t' << delta_R << '\t' << koeff << '\n';

	while (delta_R >= 1731.1)
	{
		delta_R_last = delta_R;
		Sputnik_last = Sputnik;
		Nevyazki_last = Nevyazki;
		// первый столбец матрицы Якоби (пока считаю сразу приращение по РААН, а не по времени старта)
		RAAN_Sputnik = Sputnik;
		RAAN_Sputnik.RAAN += delta_RAAN;
		RAAN_Sputnik.from_Kepl_to_Dec();
		Nominal_trajectory = Integration_to_Moon(RAAN_Sputnik, flight_time, h_near_earth, h_outer_space, h_near_moon);
		delta_X_RAAN_plus = X_target - Nominal_trajectory.back().get_parameters()[0]; //посчитал невязки
		delta_Y_RAAN_plus = Y_target - Nominal_trajectory.back().get_parameters()[1];
		delta_Z_RAAN_plus = Z_target - Nominal_trajectory.back().get_parameters()[2];
		RAAN_Sputnik.RAAN -= 2 * delta_RAAN;
		RAAN_Sputnik.from_Kepl_to_Dec();
		Nominal_trajectory = Integration_to_Moon(RAAN_Sputnik, flight_time, h_near_earth, h_outer_space, h_near_moon);
		delta_X_RAAN_minus = X_target - Nominal_trajectory.back().get_parameters()[0]; //посчитал невязки
		delta_Y_RAAN_minus = Y_target - Nominal_trajectory.back().get_parameters()[1];
		delta_Z_RAAN_minus = Z_target - Nominal_trajectory.back().get_parameters()[2];
		Jacoby.matr[0][0] = (delta_X_RAAN_plus - delta_X_RAAN_minus) / (2 * delta_RAAN);
		Jacoby.matr[1][0] = (delta_Y_RAAN_plus - delta_Y_RAAN_minus) / (2 * delta_RAAN);
		Jacoby.matr[2][0] = (delta_Z_RAAN_plus - delta_Z_RAAN_minus) / (2 * delta_RAAN);

		// второй столбец матрицы Якоби
		tau_Sputnik = Sputnik;
		tau_Sputnik.delta_t_depart += delta_tau;
		Nominal_trajectory = Integration_to_Moon(tau_Sputnik, flight_time, h_near_earth, h_outer_space, h_near_moon);
		delta_X_tau_plus = X_target - Nominal_trajectory.back().get_parameters()[0]; //посчитал невязки
		delta_Y_tau_plus = Y_target - Nominal_trajectory.back().get_parameters()[1];
		delta_Z_tau_plus = Z_target - Nominal_trajectory.back().get_parameters()[2];
		tau_Sputnik.delta_t_depart -= 2 * delta_tau;
		Nominal_trajectory = Integration_to_Moon(tau_Sputnik, flight_time, h_near_earth, h_outer_space, h_near_moon);
		delta_X_tau_minus = X_target - Nominal_trajectory.back().get_parameters()[0]; //посчитал невязки
		delta_Y_tau_minus = Y_target - Nominal_trajectory.back().get_parameters()[1];
		delta_Z_tau_minus = Z_target - Nominal_trajectory.back().get_parameters()[2];
		Jacoby.matr[0][1] = (delta_X_tau_plus - delta_X_tau_minus) / (2 * delta_tau);
		Jacoby.matr[1][1] = (delta_Y_tau_plus - delta_Y_tau_minus) / (2 * delta_tau);
		Jacoby.matr[2][1] = (delta_Z_tau_plus - delta_Z_tau_minus) / (2 * delta_tau);

		// третий столбик матрицы Якоби
		impulse_Sputnik = Sputnik;
		impulse_Sputnik.fly_away_impulse += delta_impulse;
		Nominal_trajectory = Integration_to_Moon(impulse_Sputnik, flight_time, h_near_earth, h_outer_space, h_near_moon);
		delta_X_impulse_plus = X_target - Nominal_trajectory.back().get_parameters()[0]; //посчитал невязки
		delta_Y_impulse_plus = Y_target - Nominal_trajectory.back().get_parameters()[1];
		delta_Z_impulse_plus = Z_target - Nominal_trajectory.back().get_parameters()[2];
		impulse_Sputnik.fly_away_impulse -= 2 * delta_impulse;
		Nominal_trajectory = Integration_to_Moon(impulse_Sputnik, flight_time, h_near_earth, h_outer_space, h_near_moon);
		delta_X_impulse_minus = X_target - Nominal_trajectory.back().get_parameters()[0]; //посчитал невязки
		delta_Y_impulse_minus = Y_target - Nominal_trajectory.back().get_parameters()[1];
		delta_Z_impulse_minus = Z_target - Nominal_trajectory.back().get_parameters()[2];
		Jacoby.matr[0][2] = (delta_X_impulse_plus - delta_X_impulse_minus) / (2 * delta_impulse);
		Jacoby.matr[1][2] = (delta_Y_impulse_plus - delta_Y_impulse_minus) / (2 * delta_impulse);
		Jacoby.matr[2][2] = (delta_Z_impulse_plus - delta_Z_impulse_minus) / (2 * delta_impulse);

		matrix J_obr = Jacoby.Obr_matr();
		Variations = J_obr * Nevyazki;

		Sputnik.RAAN -= koeff * Variations[0];
		Sputnik.delta_t_depart -= koeff * Variations[1];
		Sputnik.fly_away_impulse -= koeff * Variations[2];
		Sputnik.from_Kepl_to_Dec();
		Nominal_trajectory = Integration_to_Moon(Sputnik, flight_time, h_near_earth, h_outer_space, h_near_moon);
		if (Nominal_trajectory.back().main_center_flag == "Moon")
		{
			Nevyazki[0] = X_target - Nominal_trajectory.back().get_parameters()[0]; //посчитал невязки
			Nevyazki[1] = Y_target - Nominal_trajectory.back().get_parameters()[1];
			Nevyazki[2] = Z_target - Nominal_trajectory.back().get_parameters()[2];
		}
		if (Nominal_trajectory.back().main_center_flag == "Earth")
		{
			Date_time_ms Date_time = Nominal_trajectory.back().Date_time;
			std::vector <double> Moon_state(6); // вектор параметров Луны
			double X_M, Y_M, Z_M;
			//double X_M0, Y_M0, Z_M0;

			Moon_state = Nominal_trajectory.back().get_planet_parameters(Date_time.get_string().c_str(), "MOON", "J2000", "EARTH");
			X_M = Moon_state[0]; // вектор Земля - Луна
			Y_M = Moon_state[1];
			Z_M = Moon_state[2];
			Nevyazki[0] = X_target - (Nominal_trajectory.back().get_parameters()[0] - X_M);
			Nevyazki[1] = Y_target - (Nominal_trajectory.back().get_parameters()[1] - Y_M);
			Nevyazki[2] = Z_target - (Nominal_trajectory.back().get_parameters()[2] - Z_M);
		}
		delta_R = sqrt(Nevyazki[0] * Nevyazki[0] + Nevyazki[1] * Nevyazki[1] + Nevyazki[2] * Nevyazki[2]);
		if (koeff >= 0.07)
		{
			if (delta_R < delta_R_last)
			{
				koeff *= 2;
				if (koeff > 1) { koeff = 1; }
				k++;
				fout << k << '\t' << Nevyazki[0] << '\t' << Nevyazki[1] << '\t' << Nevyazki[2] << '\t' << delta_R << '\t' << koeff << '\n';
			}
			else
			{
				koeff /= 2;
				Sputnik = Sputnik_last;
				delta_R = delta_R_last;
				Nevyazki = Nevyazki_last;
				k++;
			}
		}
		else
		{
			if (delta_R < delta_R_last)
			{
				koeff *= 2;
				if (koeff > 1) { koeff = 1; }
				k++;
				fout << k << '\t' << Nevyazki[0] << '\t' << Nevyazki[1] << '\t' << Nevyazki[2] << '\t' << delta_R << '\t' << koeff << '\n';
			}
			else { k++; }
		}
	}
}

void Newtons_method_b_plane(Spacecraft& Sputnik, double flight_time, double h_near_earth, double h_outer_space, double h_near_moon, double Orbit_altitude)
{
	Spacecraft Nominal_Sputnik, RAAN_Sputnik, tau_Sputnik, impulse_Sputnik, Sputnik_last;
	std::vector <Spacecraft> Nominal_trajectory; //номинальная траектория (и другие)
	// у спутника уже должны быть первые приближения

	std::vector <double> Nevyazki(3), Nevyazki_last(3), Variations(3);
	matrix Jacoby(3, 3); // матрица Якоби
	double delta_Be_RAAN_plus, delta_Bn_RAAN_plus, delta_t_RAAN_plus,
		delta_Be_RAAN_minus, delta_Bn_RAAN_minus, delta_t_RAAN_minus;
	double delta_Be_tau_plus, delta_Bn_tau_plus, delta_t_tau_plus,
		delta_Be_tau_minus, delta_Bn_tau_minus, delta_t_tau_minus;
	double delta_Be_impulse_plus, delta_Bn_impulse_plus, delta_t_impulse_plus,
		delta_Be_impulse_minus, delta_Bn_impulse_minus, delta_t_impulse_minus;

	// один раз интегрирую, чтобы получить модуль скорости на  бесконечности
	Nominal_trajectory = Integration_to_Moon_periselene(Sputnik, h_near_earth, h_outer_space, h_near_moon);
	double V_inf = Nominal_trajectory.back().V_inf;
	double R_periselene = Orbit_altitude + 1737.1; // радиус точке периселения
	double abs_B_target = sqrt(R_periselene * R_periselene + 2 * R_periselene * mu_Moon / (V_inf * V_inf));
	double Period_target = 2 * Pi * sqrt(pow(R_periselene, 3) / mu_Moon);

	double Be_target = 0, Bn_target = -abs_B_target, time_target = flight_time * 24 * 60 * 60; // задаю условия, в которые нужно попасть
	double delta_RAAN = 0.00000001, delta_tau = 0.01 /*в секундах*/, delta_impulse = 0.0000001; // шаги для производных
	double delta_F, delta_F_last;
	double koeff = 1; // коэффициент для прибавления приращений
	int k = 0; // номер итерации

	std::ofstream fout("Newton_b_plane.txt");
	fout << "Номер итерации" << '\t' << "Дельта Be" << '\t' << "Дельта Bn" << '\t' << "Дельта time" << '\t' << "Дельта F" << '\t' << "Коэффициент K" << '\n';
	Nominal_Sputnik = Sputnik;
	Nominal_trajectory = Integration_to_Moon_periselene(Nominal_Sputnik, h_near_earth, h_outer_space, h_near_moon);
	if (Nominal_trajectory.back().main_center_flag == "Moon")
	{
		Nevyazki[0] = Be_target - Nominal_trajectory.back().get_b_plane_parameters()[0]; //посчитал невязки
		Nevyazki[1] = Bn_target - Nominal_trajectory.back().get_b_plane_parameters()[1];
		Nevyazki[2] = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
	}
	if (Nominal_trajectory.back().main_center_flag == "Earth")
	{
		Nevyazki[0] = Be_target - Nominal_trajectory.back().get_b_plane_parameters()[0]; //посчитал невязки
		Nevyazki[1] = Bn_target - Nominal_trajectory.back().get_b_plane_parameters()[1];
		Nevyazki[2] = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
	}
	// Считаю модуль вектора невязок для оценивания сходимости, для этого его надо обезразмерить
	// Вектор прицельной дальности поделю на радиус прицельной орбиты ИСЛ, а время на период этой орбиты
	delta_F = sqrt(Nevyazki[0] * Nevyazki[0] / (R_periselene * R_periselene) + Nevyazki[1] * Nevyazki[1] / (R_periselene * R_periselene)
		+ Nevyazki[2] * Nevyazki[2] / (Period_target * Period_target));

	fout << k << '\t' << Nevyazki[0] << '\t' << Nevyazki[1] << '\t' << Nevyazki[2] << '\t' << delta_F << '\t' << koeff << '\n';

	while (delta_F >= 0.05) //0.005
	{
		delta_F_last = delta_F;
		Sputnik_last = Sputnik;
		Nevyazki_last = Nevyazki;
		// первый столбец матрицы Якоби (пока считаю сразу приращение по РААН, а не по времени старта)
		RAAN_Sputnik = Sputnik;
		RAAN_Sputnik.RAAN += delta_RAAN;
		RAAN_Sputnik.from_Kepl_to_Dec();
		Nominal_trajectory = Integration_to_Moon_periselene(RAAN_Sputnik, h_near_earth, h_outer_space, h_near_moon);
		delta_Be_RAAN_plus = Be_target - Nominal_trajectory.back().get_b_plane_parameters()[0]; //посчитал невязки
		delta_Bn_RAAN_plus = Bn_target - Nominal_trajectory.back().get_b_plane_parameters()[1];
		delta_t_RAAN_plus = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
		RAAN_Sputnik.RAAN -= 2 * delta_RAAN;
		RAAN_Sputnik.from_Kepl_to_Dec();
		Nominal_trajectory = Integration_to_Moon_periselene(RAAN_Sputnik, h_near_earth, h_outer_space, h_near_moon);
		delta_Be_RAAN_minus = Be_target - Nominal_trajectory.back().get_b_plane_parameters()[0]; //посчитал невязки
		delta_Bn_RAAN_minus = Bn_target - Nominal_trajectory.back().get_b_plane_parameters()[1];
		delta_t_RAAN_minus = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
		Jacoby.matr[0][0] = (delta_Be_RAAN_plus - delta_Be_RAAN_minus) / (2 * delta_RAAN);
		Jacoby.matr[1][0] = (delta_Bn_RAAN_plus - delta_Bn_RAAN_minus) / (2 * delta_RAAN);
		Jacoby.matr[2][0] = (delta_t_RAAN_plus - delta_t_RAAN_minus) / (2 * delta_RAAN);

		// второй столбец матрицы Якоби
		tau_Sputnik = Sputnik;
		tau_Sputnik.delta_t_depart += delta_tau;
		Nominal_trajectory = Integration_to_Moon_periselene(tau_Sputnik, h_near_earth, h_outer_space, h_near_moon);
		delta_Be_tau_plus = Be_target - Nominal_trajectory.back().get_b_plane_parameters()[0]; //посчитал невязки
		delta_Bn_tau_plus = Bn_target - Nominal_trajectory.back().get_b_plane_parameters()[1];
		delta_t_tau_plus = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
		tau_Sputnik.delta_t_depart -= 2 * delta_tau;
		Nominal_trajectory = Integration_to_Moon_periselene(tau_Sputnik, h_near_earth, h_outer_space, h_near_moon);
		delta_Be_tau_minus = Be_target - Nominal_trajectory.back().get_b_plane_parameters()[0]; //посчитал невязки
		delta_Bn_tau_minus = Bn_target - Nominal_trajectory.back().get_b_plane_parameters()[1];
		delta_t_tau_minus = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
		Jacoby.matr[0][1] = (delta_Be_tau_plus - delta_Be_tau_minus) / (2 * delta_tau);
		Jacoby.matr[1][1] = (delta_Bn_tau_plus - delta_Bn_tau_minus) / (2 * delta_tau);
		Jacoby.matr[2][1] = (delta_t_tau_plus - delta_t_tau_minus) / (2 * delta_tau);

		// третий столбик матрицы Якоби
		impulse_Sputnik = Sputnik;
		impulse_Sputnik.fly_away_impulse += delta_impulse;
		Nominal_trajectory = Integration_to_Moon_periselene(impulse_Sputnik, h_near_earth, h_outer_space, h_near_moon);
		delta_Be_impulse_plus = Be_target - Nominal_trajectory.back().get_b_plane_parameters()[0]; //посчитал невязки
		delta_Bn_impulse_plus = Bn_target - Nominal_trajectory.back().get_b_plane_parameters()[1];
		delta_t_impulse_plus = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
		impulse_Sputnik.fly_away_impulse -= 2 * delta_impulse;
		Nominal_trajectory = Integration_to_Moon_periselene(impulse_Sputnik, h_near_earth, h_outer_space, h_near_moon);
		delta_Be_impulse_minus = Be_target - Nominal_trajectory.back().get_b_plane_parameters()[0]; //посчитал невязки
		delta_Bn_impulse_minus = Bn_target - Nominal_trajectory.back().get_b_plane_parameters()[1];
		delta_t_impulse_minus = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
		Jacoby.matr[0][2] = (delta_Be_impulse_plus - delta_Be_impulse_minus) / (2 * delta_impulse);
		Jacoby.matr[1][2] = (delta_Bn_impulse_plus - delta_Bn_impulse_minus) / (2 * delta_impulse);
		Jacoby.matr[2][2] = (delta_t_impulse_plus - delta_t_impulse_minus) / (2 * delta_impulse);

		matrix J_obr = Jacoby.Obr_matr();
		Variations = J_obr * Nevyazki;

		Sputnik.RAAN -= koeff * Variations[0];
		Sputnik.delta_t_depart -= koeff * Variations[1];
		Sputnik.fly_away_impulse -= koeff * Variations[2];
		Sputnik.from_Kepl_to_Dec();
		Nominal_trajectory = Integration_to_Moon_periselene(Sputnik, h_near_earth, h_outer_space, h_near_moon);
		if (Nominal_trajectory.back().main_center_flag == "Moon")
		{
			Nevyazki[0] = Be_target - Nominal_trajectory.back().get_b_plane_parameters()[0]; //посчитал невязки
			Nevyazki[1] = Bn_target - Nominal_trajectory.back().get_b_plane_parameters()[1];
			Nevyazki[2] = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
		}
		if (Nominal_trajectory.back().main_center_flag == "Earth")
		{
			Nevyazki[0] = Be_target - Nominal_trajectory.back().get_b_plane_parameters()[0]; //посчитал невязки
			Nevyazki[1] = Bn_target - Nominal_trajectory.back().get_b_plane_parameters()[1];
			Nevyazki[2] = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
		}
		// Считаю модуль вектора невязок для оценивания сходимости, для этого его надо обезразмерить
		// Вектор прицельной дальности поделю на радиус прицельной орбиты ИСЛ, а время на период этой орбиты
		delta_F = sqrt(Nevyazki[0] * Nevyazki[0] / (R_periselene * R_periselene) + Nevyazki[1] * Nevyazki[1] / (R_periselene * R_periselene)
			+ Nevyazki[2] * Nevyazki[2] / (Period_target * Period_target));
		if (koeff >= 0.03)
		{
			if (delta_F < delta_F_last)
			{
				koeff *= 2;
				if (koeff > 1) { koeff = 1; }
				k++;
				fout << k << '\t' << Nevyazki[0] << '\t' << Nevyazki[1] << '\t' << Nevyazki[2] << '\t' << delta_F << '\t' << koeff << '\n';
			}
			else
			{
				koeff /= 2;
				Sputnik = Sputnik_last;
				delta_F = delta_F_last;
				Nevyazki = Nevyazki_last;
				k++;
			}
		}
		else
		{
			if (delta_F < delta_F_last)
			{
				koeff *= 2;
				if (koeff > 1) { koeff = 1; }
				k++;
				fout << k << '\t' << Nevyazki[0] << '\t' << Nevyazki[1] << '\t' << Nevyazki[2] << '\t' << delta_F << '\t' << koeff << '\n';
			}
			else { k++; }
		}
	}
}

void Newtons_method_final_step(Spacecraft& Sputnik, double flight_time, double h_near_earth, double h_outer_space, double h_near_moon, double Orbit_altitude)
{
	Spacecraft Nominal_Sputnik, RAAN_Sputnik, tau_Sputnik, impulse_Sputnik, Sputnik_last;
	std::vector <Spacecraft> Nominal_trajectory; //номинальная траектория (и другие)
	// у спутника уже должны быть первые приближения

	std::vector <double> Nevyazki(3), Nevyazki_last(3), Variations(3);
	matrix Jacoby(3, 3); // матрица Якоби
	double delta_R_RAAN_plus, delta_inc_RAAN_plus, delta_t_RAAN_plus,
		delta_R_RAAN_minus, delta_inc_RAAN_minus, delta_t_RAAN_minus;
	double delta_R_tau_plus, delta_inc_tau_plus, delta_t_tau_plus,
		delta_R_tau_minus, delta_inc_tau_minus, delta_t_tau_minus;
	double delta_R_impulse_plus, delta_inc_impulse_plus, delta_t_impulse_plus,
		delta_R_impulse_minus, delta_inc_impulse_minus, delta_t_impulse_minus;

	double R_target = Orbit_altitude + 1737.1, inc_target = 90 * Pi / 180, time_target = flight_time * 24 * 60 * 60; // задаю условия, в которые нужно попасть
	double Period_target = 2 * Pi * sqrt(pow(R_target, 3) / mu_Moon);
	double delta_RAAN = 0.00000001, delta_tau = 0.01 /*в секундах*/, delta_impulse = 0.0000001; // шаги для производных
	double delta_F, delta_F_last;
	double koeff = 1; // коэффициент для прибавления приращений
	int k = 0; // номер итерации
	std::vector <double> SVSK_crds(6);

	std::ofstream fout("Newton_final_step.txt");
	fout << "Номер итерации" << '\t' << "Дельта R" << '\t' << "Дельта i" << '\t' << "Дельта time" << '\t' << "Дельта F" << '\t' << "Коэффициент K" << '\n';
	Nominal_Sputnik = Sputnik;
	Nominal_trajectory = Integration_to_Moon_periselene(Nominal_Sputnik, h_near_earth, h_outer_space, h_near_moon);
	if (Nominal_trajectory.back().main_center_flag == "Moon")
	{
		Nevyazki[0] = R_target - Nominal_trajectory.back().r_KA; //посчитал невязки
		Nevyazki[1] = inc_target - Nominal_trajectory.back().SVSK_inclination();
		Nevyazki[2] = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
	}
	if (Nominal_trajectory.back().main_center_flag == "Earth")
	{
		Nevyazki[0] = R_target - Nominal_trajectory.back().distance_to_Moon(); //посчитал невязки
		Nevyazki[1] = inc_target - Nominal_trajectory.back().SVSK_inclination();
		Nevyazki[2] = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
	}
	// Считаю модуль вектора невязок для оценивания сходимости, для этого его надо обезразмерить
	// Невязку по радиусу делю на радиус прицельной орбиты ИСЛ, невязку наклонения к прицельному наклонению, а время на период этой орбиты
	delta_F = sqrt(Nevyazki[0] * Nevyazki[0] / (R_target * R_target) + Nevyazki[1] * Nevyazki[1] / (inc_target * inc_target)
		+ Nevyazki[2] * Nevyazki[2] / (Period_target * Period_target));

	fout << k << '\t' << Nevyazki[0] << '\t' << Nevyazki[1] << '\t' << Nevyazki[2] << '\t' << delta_F << '\t' << koeff << '\n';

	while (delta_F >= 0.0001) //(delta_F >= 0.0001)
	{
		delta_F_last = delta_F;
		Sputnik_last = Sputnik;
		Nevyazki_last = Nevyazki;
		// первый столбец матрицы Якоби (пока считаю сразу приращение по РААН, а не по времени старта)
		RAAN_Sputnik = Sputnik;
		RAAN_Sputnik.RAAN += delta_RAAN;
		RAAN_Sputnik.from_Kepl_to_Dec();
		Nominal_trajectory = Integration_to_Moon_periselene(RAAN_Sputnik, h_near_earth, h_outer_space, h_near_moon);
		delta_R_RAAN_plus = R_target - Nominal_trajectory.back().r_KA; //посчитал невязки
		delta_inc_RAAN_plus = inc_target - Nominal_trajectory.back().SVSK_inclination();
		delta_t_RAAN_plus = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
		RAAN_Sputnik.RAAN -= 2 * delta_RAAN;
		RAAN_Sputnik.from_Kepl_to_Dec();
		Nominal_trajectory = Integration_to_Moon_periselene(RAAN_Sputnik, h_near_earth, h_outer_space, h_near_moon);
		delta_R_RAAN_minus = R_target - Nominal_trajectory.back().r_KA; //посчитал невязки
		delta_inc_RAAN_minus = inc_target - Nominal_trajectory.back().SVSK_inclination();
		delta_t_RAAN_minus = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
		Jacoby.matr[0][0] = (delta_R_RAAN_plus - delta_R_RAAN_minus) / (2 * delta_RAAN);
		Jacoby.matr[1][0] = (delta_inc_RAAN_plus - delta_inc_RAAN_minus) / (2 * delta_RAAN);
		Jacoby.matr[2][0] = (delta_t_RAAN_plus - delta_t_RAAN_minus) / (2 * delta_RAAN);

		// второй столбец матрицы Якоби
		tau_Sputnik = Sputnik;
		tau_Sputnik.delta_t_depart += delta_tau;
		Nominal_trajectory = Integration_to_Moon_periselene(tau_Sputnik, h_near_earth, h_outer_space, h_near_moon);
		delta_R_tau_plus = R_target - Nominal_trajectory.back().r_KA; //посчитал невязки
		delta_inc_tau_plus = inc_target - Nominal_trajectory.back().SVSK_inclination();
		delta_t_tau_plus = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
		tau_Sputnik.delta_t_depart -= 2 * delta_tau;
		Nominal_trajectory = Integration_to_Moon_periselene(tau_Sputnik, h_near_earth, h_outer_space, h_near_moon);
		delta_R_tau_minus = R_target - Nominal_trajectory.back().r_KA; //посчитал невязки
		delta_inc_tau_minus = inc_target - Nominal_trajectory.back().SVSK_inclination();
		delta_t_tau_minus = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
		Jacoby.matr[0][1] = (delta_R_tau_plus - delta_R_tau_minus) / (2 * delta_tau);
		Jacoby.matr[1][1] = (delta_inc_tau_plus - delta_inc_tau_minus) / (2 * delta_tau);
		Jacoby.matr[2][1] = (delta_t_tau_plus - delta_t_tau_minus) / (2 * delta_tau);

		// третий столбик матрицы Якоби
		impulse_Sputnik = Sputnik;
		impulse_Sputnik.fly_away_impulse += delta_impulse;
		Nominal_trajectory = Integration_to_Moon_periselene(impulse_Sputnik, h_near_earth, h_outer_space, h_near_moon);
		delta_R_impulse_plus = R_target - Nominal_trajectory.back().r_KA; //посчитал невязки
		delta_inc_impulse_plus = inc_target - Nominal_trajectory.back().SVSK_inclination();
		delta_t_impulse_plus = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
		impulse_Sputnik.fly_away_impulse -= 2 * delta_impulse;
		Nominal_trajectory = Integration_to_Moon_periselene(impulse_Sputnik, h_near_earth, h_outer_space, h_near_moon);
		delta_R_impulse_minus = R_target - Nominal_trajectory.back().r_KA; //посчитал невязки
		delta_inc_impulse_minus = inc_target - Nominal_trajectory.back().SVSK_inclination();
		delta_t_impulse_minus = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
		Jacoby.matr[0][2] = (delta_R_impulse_plus - delta_R_impulse_minus) / (2 * delta_impulse);
		Jacoby.matr[1][2] = (delta_inc_impulse_plus - delta_inc_impulse_minus) / (2 * delta_impulse);
		Jacoby.matr[2][2] = (delta_t_impulse_plus - delta_t_impulse_minus) / (2 * delta_impulse);

		matrix J_obr = Jacoby.Obr_matr();
		Variations = J_obr * Nevyazki;

		Sputnik.RAAN -= koeff * Variations[0];
		Sputnik.delta_t_depart -= koeff * Variations[1];
		Sputnik.fly_away_impulse -= koeff * Variations[2];
		Sputnik.from_Kepl_to_Dec();
		Nominal_trajectory = Integration_to_Moon_periselene(Sputnik, h_near_earth, h_outer_space, h_near_moon);
		if (Nominal_trajectory.back().main_center_flag == "Moon")
		{
			Nevyazki[0] = R_target - Nominal_trajectory.back().r_KA; //посчитал невязки
			Nevyazki[1] = inc_target - Nominal_trajectory.back().SVSK_inclination();
			Nevyazki[2] = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
		}
		if (Nominal_trajectory.back().main_center_flag == "Earth")
		{
			Nevyazki[0] = R_target - Nominal_trajectory.back().r_KA; //посчитал невязки
			Nevyazki[1] = inc_target - Nominal_trajectory.back().SVSK_inclination();
			Nevyazki[2] = time_target - (Nominal_trajectory.back().time - Nominal_trajectory.back().delta_t_depart);
			std::cout << "ERROR";
		}
		// Считаю модуль вектора невязок для оценивания сходимости, для этого его надо обезразмерить
		// Вектор прицельной дальности поделю на радиус прицельной орбиты ИСЛ, а время на период этой орбиты
		delta_F = sqrt(Nevyazki[0] * Nevyazki[0] / (R_target * R_target) + Nevyazki[1] * Nevyazki[1] / (inc_target * inc_target)
			+ Nevyazki[2] * Nevyazki[2] / (Period_target * Period_target));
		if (koeff >= 0.03)
		{
			if (delta_F < delta_F_last)
			{
				koeff *= 2;
				if (koeff > 1) { koeff = 1; }
				k++;
				fout << k << '\t' << Nevyazki[0] << '\t' << Nevyazki[1] << '\t' << Nevyazki[2] << '\t' << delta_F << '\t' << koeff << '\n';
			}
			else
			{
				koeff /= 2;
				Sputnik = Sputnik_last;
				delta_F = delta_F_last;
				Nevyazki = Nevyazki_last;
				k++;
			}
		}
		else
		{
			if (delta_F < delta_F_last)
			{
				koeff *= 2;
				if (koeff > 1) { koeff = 1; }
				k++;
				fout << k << '\t' << Nevyazki[0] << '\t' << Nevyazki[1] << '\t' << Nevyazki[2] << '\t' << delta_F << '\t' << koeff << '\n';
			}
			else { k++; }
		}
	}
}

std::vector <Spacecraft> Orbiting_Moon(Spacecraft Sputnik, double altitude, double integration_time, double h_near_moon) //торможение и движение у Луны
{
	std::vector <Spacecraft> SC_Condition;
	SC_Condition.push_back(Sputnik);
	double thrust_impulse = sqrt(mu_Moon / (1737.1 + altitude)) - Sputnik.V_KA;
	Sputnik.Collinear_impulse(thrust_impulse);
	Sputnik.from_Dec_to_Kepl();
	Sputnik.time = 0;
	SC_Condition.push_back(Sputnik);

	while (Sputnik.time <= integration_time - h_near_moon)
	{
		Sputnik.Runge_Kutta(h_near_moon);
		SC_Condition.push_back(Sputnik);
	}
	return SC_Condition;
}

void NRHO_parameters(double synodic_resonance, double r_peri, double& r_apo, double& Orbit_Altitude_peri, double& Orbit_Altitude_apo, double& Orbit_period, double& V_peri, double& V_apo)
{
	Orbit_period = Moon_synodic_period / synodic_resonance * 86400; // Орбитальный период в сек
	double a = cbrt(mu_Moon * pow(Orbit_period / (2 * Pi), 2));
	r_apo = 2 * a - r_peri;
	Orbit_Altitude_peri = r_peri - 1737.1;
	Orbit_Altitude_apo = r_apo - 1737.1;
	V_peri = sqrt(mu_Moon * (2 / r_peri - 1 / a));
	V_apo = sqrt(mu_Moon * (2 / r_apo - 1 / a));
}

double Precise_NRHO_impulse_aposelene(Spacecraft Sputnik, double NRHO_r_apo, double V_peri0, double h)
{
	std::vector <Spacecraft> SC_Condition;
	double thrust_impulse = 0, r_aposelene = 0;
	double X0 = 0, Y0 = 0, X1 = 0, Y1 = 0, X2 = 0, Y2 = 0;

	thrust_impulse = V_peri0 - Sputnik.V_KA;
	Sputnik.Collinear_impulse(thrust_impulse);
	Sputnik.from_Dec_to_Kepl();
	SC_Condition = Integration_to_Moon_aposelene(Sputnik, h);
	Sputnik.Collinear_impulse(-thrust_impulse); // возвращаю спутник в начальное положение
	Sputnik.from_Dec_to_Kepl();
	r_aposelene = SC_Condition.back().distance_to_Moon();
	X1 = thrust_impulse;
	Y1 = r_aposelene - NRHO_r_apo;

	thrust_impulse += 0.01;
	Sputnik.Collinear_impulse(thrust_impulse);
	Sputnik.from_Dec_to_Kepl();
	SC_Condition = Integration_to_Moon_aposelene(Sputnik, h);
	Sputnik.Collinear_impulse(-thrust_impulse); // возвращаю спутник в начальное положение
	Sputnik.from_Dec_to_Kepl();
	r_aposelene = SC_Condition.back().distance_to_Moon();
	X2 = thrust_impulse;
	Y2 = r_aposelene - NRHO_r_apo;
	do
	{
		X0 = Y1 / (Y1 - Y2) * (X2 - X1) + X1;
		thrust_impulse = X0;
		Sputnik.Collinear_impulse(thrust_impulse);
		Sputnik.from_Dec_to_Kepl();
		SC_Condition = Integration_to_Moon_aposelene(Sputnik, h);
		Sputnik.Collinear_impulse(-thrust_impulse); // возвращаю спутник в начальное положение
		Sputnik.from_Dec_to_Kepl();
		r_aposelene = SC_Condition.back().distance_to_Moon();
		Y0 = r_aposelene - NRHO_r_apo;
		if (abs(Y1) > abs(Y2))
		{
			X1 = X0;
			Y1 = Y0;
		}
		else
		{
			X2 = X0;
			Y2 = Y0;
		}
	} while (abs(Y0) > pow(10, -1) and (abs(X1 - X2) > 0.00001));
	return thrust_impulse;
}

double Precise_NRHO_impulse_period(Spacecraft Sputnik, double NRHO_period, double V_peri0, double h)
{
	std::vector <Spacecraft> SC_Condition;
	double thrust_impulse = 0, time = 0;
	double X0 = 0, Y0 = 0, X1 = 0, Y1 = 0, X2 = 0, Y2 = 0;

	thrust_impulse = V_peri0 - Sputnik.V_KA;
	Sputnik.Collinear_impulse(thrust_impulse);
	Sputnik.from_Dec_to_Kepl();
	SC_Condition = Integration_full_period(Sputnik, h);
	Sputnik.Collinear_impulse(-thrust_impulse); // возвращаю спутник в начальное положение
	Sputnik.from_Dec_to_Kepl();
	time = SC_Condition.back().time;
	X1 = thrust_impulse;
	Y1 = time - NRHO_period;

	thrust_impulse -= 0.002;
	Sputnik.Collinear_impulse(thrust_impulse);
	Sputnik.from_Dec_to_Kepl();
	SC_Condition = Integration_full_period(Sputnik, h);
	Sputnik.Collinear_impulse(-thrust_impulse); // возвращаю спутник в начальное положение
	Sputnik.from_Dec_to_Kepl();
	time = SC_Condition.back().time;
	X2 = thrust_impulse;
	Y2 = time - NRHO_period;
	do
	{
		X0 = Y1 / (Y1 - Y2) * (X2 - X1) + X1;
		thrust_impulse = X0;
		Sputnik.Collinear_impulse(thrust_impulse);
		Sputnik.from_Dec_to_Kepl();
		SC_Condition = Integration_full_period(Sputnik, h);
		Sputnik.Collinear_impulse(-thrust_impulse); // возвращаю спутник в начальное положение
		Sputnik.from_Dec_to_Kepl();
		time = SC_Condition.back().time;
		Y0 = time - NRHO_period;
		if (abs(Y1) > abs(Y2))
		{
			X1 = X0;
			Y1 = Y0;
		}
		else
		{
			X2 = X0;
			Y2 = Y0;
		}
	} while (abs(Y0) > pow(10, 1) and (abs(X1 - X2) > 0.00001));
	return thrust_impulse;
}

std::vector <Spacecraft> NRHO_flight(Spacecraft Sputnik, double NRHO_V, double NRHO_Orbit_period, double integration_time_periods, double h_near_moon) //торможение и полёт на NRHO от Апоселения
{
	std::vector <Spacecraft> SC_Condition;
	SC_Condition.push_back(Sputnik);
	double thrust_impulse = NRHO_V - Sputnik.V_KA;
	Sputnik.Collinear_impulse(thrust_impulse);
	Sputnik.from_Dec_to_Kepl();
	Sputnik.time = 0;
	SC_Condition.push_back(Sputnik);

	double integration_time = integration_time_periods * NRHO_Orbit_period;

	while (Sputnik.time <= integration_time - h_near_moon)
	{
		Sputnik.Runge_Kutta(h_near_moon);
		SC_Condition.push_back(Sputnik);
	}
	return SC_Condition;
}

std::vector <Spacecraft> Integration_R(Spacecraft Sputnik, double R_target, double h)
{
	std::vector <Spacecraft> SC_Condition;
	SC_Condition.push_back(Sputnik);
	while (Sputnik.r_KA <= R_target)
	{
		Sputnik.Runge_Kutta(h);
		SC_Condition.push_back(Sputnik);
	}
	return SC_Condition;
}

std::vector <Spacecraft> Integration_back(Spacecraft Sputnik, double integration_time, double h) // здесь h положительное
{
	std::vector <Spacecraft> SC_Condition;
	SC_Condition.push_back(Sputnik);
	while (Sputnik.time >= -integration_time + h)
	{
		Sputnik.Runge_Kutta(-h);
		SC_Condition.push_back(Sputnik);
	}
	return SC_Condition;
}

void Check_integration_step(Spacecraft Sputnik, double integration_time, double h, double& R_diff, double& V_diff)
{
	Spacecraft Sputnik_goes_back;
	std::vector <double> Difference; Difference.resize(6);

	std::vector <Spacecraft> SC_Condition = Integration(Sputnik, integration_time, h);
	Sputnik_goes_back.set_parameters(SC_Condition.back().get_parameters());
	std::vector <Spacecraft> SC_Condition_to_0 = Integration_back(Sputnik_goes_back, integration_time, h);
	for (size_t i = 0; i < 6; i++)
	{
		Difference[i] = Sputnik.get_parameters()[i] - SC_Condition_to_0.back().get_parameters()[i];
		std::cout << "Difference[" << i << "] = " << Difference[i] << '\n';
	}
	R_diff = sqrt(Difference[0] * Difference[0] + Difference[1] * Difference[1] + Difference[2] * Difference[2]);
	V_diff = sqrt(Difference[3] * Difference[3] + Difference[4] * Difference[4] + Difference[5] * Difference[5]);
	std::cout << "R_diff = " << R_diff << '\n' << "V_diff = " << V_diff << '\n';
}

void Check_integration_step_hyper(Spacecraft Sputnik, double R, double h, double& R_diff, double& V_diff)
{
	Spacecraft Sputnik_goes_back;
	std::vector <double> Difference; Difference.resize(6);

	std::vector <Spacecraft> SC_Condition = Integration_R(Sputnik, R, h);
	Sputnik_goes_back.set_parameters(SC_Condition.back().get_parameters());
	std::vector <Spacecraft> SC_Condition_to_0 = Integration_back(Sputnik_goes_back, SC_Condition.back().time, h);
	for (size_t i = 0; i < 6; i++)
	{
		Difference[i] = Sputnik.get_parameters()[i] - SC_Condition_to_0.back().get_parameters()[i];
		std::cout << "Difference[" << i << "] = " << Difference[i] << '\n';
	}
	R_diff = sqrt(Difference[0] * Difference[0] + Difference[1] * Difference[1] + Difference[2] * Difference[2]);
	V_diff = sqrt(Difference[3] * Difference[3] + Difference[4] * Difference[4] + Difference[5] * Difference[5]);
	std::cout << "R_diff = " << R_diff << '\n' << "V_diff = " << V_diff << '\n';
}

void Searching_for_optimal_int_step(Spacecraft Sputnik)
{
	std::ofstream fout1("IntStep_Diff.txt");
	fout1 << "h" << '\t' << "Difference" << '\n';
	double R_diff, V_diff;
	for (int h = 5; h <= 100; h += 15)
	{
		Check_integration_step(Sputnik, 100 * Sputnik.Period, h, R_diff, V_diff);
		fout1 << h << '\t' << R_diff << '\n';
	}
}

void Searching_for_optimal_int_step_hyper(Spacecraft Sputnik)
{
	std::ofstream fout1("IntStep_Diff.txt");
	fout1 << "h" << '\t' << "Difference" << '\n';
	double R_diff, V_diff;
	for (int h = 5; h <= 100; h += 15)
	{
		Check_integration_step_hyper(Sputnik, 100, h, R_diff, V_diff);
		fout1 << h << '\t' << R_diff << '\n';
	}
}

void Write_to_Files(std::vector <Spacecraft> a, int b)
{
	std::ofstream fout1("Integ.txt");
	std::ofstream fout2("Kepler.txt");
	std::ofstream fout3("BPlane.txt");
	std::ofstream fout4("SVSK.txt");
	std::ofstream fout5("Integ_GreenwichCS.txt");
	fout1 << "Date" << '\t' << "Time" << '\t' << "X" << '\t' << "Y" << '\t' << "Z" << '\t' << "r_KA" << '\t' << "Расстояние до Луны" << '\t' << "Vx" << '\t' << "Vy" << '\t' << "Vz" << '\t' << "V_KA" << '\t' << "Главный центр" << '\t' << "Угол наклона тр." << '\n';
	fout2 << "Date" << '\t' << "Time" << '\t' << "a" << '\t' << "ecc" << '\t' << "inc" << '\t' << "RAAN" << '\t' << "omega" << '\t' << "tau" << '\t' << "theta" << '\t' << "Главный центр" << '\n';
	fout3 << "Date" << '\t' << "Time" << '\t' << "V_inf" << '\t' << "B_ksi" << '\t' << "B_etta" << '\n';
	fout4 << "Date" << '\t' << "Time" << '\t' << "X" << '\t' << "Y" << '\t' << "Z" << '\t' << "R" << '\t' << "Vx" << '\t' << "Vy" << '\t' << "Vz" << '\t' << "inc" << '\n';
	fout5 << "Date" << '\t' << "Time" << '\t' << "X" << '\t' << "Y" << '\t' << "Z" << '\t' << "r_KA" << '\t' << "Vx" << '\t' << "Vy" << '\t' << "Vz" << '\t' << "V_KA" << '\t' << "Угол наклона тр." << '\n';
	fout1 << "fly_away_impulse = " << a[0].fly_away_impulse << '\n';
	for (size_t i = 0; i < a.size(); i++)
	{
		if ((i % b == 0) or (i == a.size() - 1) or (a[i].time == a[i].delta_t_depart) or ((a[i].time - a[i].delta_t_depart > 0) and (a[i].time - a[i].delta_t_depart < 3000) and (i % 4 == 0)))
		{
			a[i].printParameters(fout1);
			a[i].printKepler(fout2);
			if (a[i].distance_to_Moon() < 75000)
			{
				a[i].printBPlane(fout3);
				a[i].printSVSK(fout4);
			}
			a[i].printGCS(fout5);
		}
	}
}

void Write_to_Files_Moon_orbit(std::vector <Spacecraft> a, int b)
{
	std::ofstream fout1("Integ_Moon_orbit.txt");
	std::ofstream fout2("Integ_Moon_orbit_Kepler.txt");
	std::ofstream fout3("Integ_Moon_orbit_SVSK.txt");
	std::ofstream fout4("Integ_Moon_orbit_SISK.txt");
	fout1 << "Date" << '\t' << "Time" << '\t' << "X" << '\t' << "Y" << '\t' << "Z" << '\t' << "r_KA" << '\t' << "Расстояние до Луны" << '\t' << "Vx" << '\t' << "Vy" << '\t' << "Vz" << '\t' << "V_KA" << '\t' << "Главный центр" << '\n';
	fout2 << "Date" << '\t' << "Time" << '\t' << "a" << '\t' << "ecc" << '\t' << "inc" << '\t' << "RAAN" << '\t' << "omega" << '\t' << "tau" << '\t' << "theta" << '\t' << "Главный центр" << '\n';
	fout3 << "Date" << '\t' << "Time" << '\t' << "X" << '\t' << "Y" << '\t' << "Z" << '\t' << "R" << '\t' << "Vx" << '\t' << "Vy" << '\t' << "Vz" << '\t' << "inc" << '\n';
	fout4 << "Date" << '\t' << "Time" << '\t' << "X" << '\t' << "Y" << '\t' << "Z" << '\t' << "R" << '\t' << "Vx" << '\t' << "Vy" << '\t' << "Vz" << '\t' << "inc" << '\n';
	for (size_t i = 0; i < a.size(); i++)
	{
		if ((i % b == 0) or (i == a.size() - 1) or (i == 1))
		{
			a[i].printParameters_relatively_to_Moon(fout1);
			a[i].printKepler_relatively_to_Moon(fout2);
			a[i].printSVSK_relatively_to_Moon(fout3);
			a[i].printSISK_relatively_to_Moon(fout4);
		}
	}
}

void Write_to_Files_initial_cond(std::vector <Spacecraft> a, std::vector <Spacecraft> b)
{
	std::ofstream fout("Initial_conditions.txt"); // файл для н.у. и условий во время подлёта к Луне
	std::cout << std::setprecision(8);
	fout << "Все параматеры представлены в СК J2000" << '\n';
	fout << "Date" << '\t' << "Time" << '\t' << "X" << '\t' << "Y" << '\t' << "Z" << '\t' << "r_KA" << '\t' << "Расстояние до Луны" << '\t' << "Vx" << '\t' << "Vy" << '\t' << "Vz" << '\t' << "V_KA" << '\t' << "Главный центр" << '\n';
	fout << "Начало интегрирования на опорной орбите:" << '\n';
	a[0].printParameters(fout);
	fout << "fly_away_impulse = " << a[0].fly_away_impulse << '\t' << "delta_t_depart = " << a[0].delta_t_depart << '\n';
	fout << "Момент подачи разгонного импульса:" << '\n';
	for (size_t i = 0; i < a.size(); i++)
	{
		if (a[i].time == a[i].delta_t_depart)
		{
			a[i].printParameters(fout);
			a[i + 1].printParameters(fout);
			break;
		}
	}
	fout << "Окончание интегрирования у Луны:" << '\n';
	b[0].printParameters(fout);
	fout << "Сразу после подачи тормозного импульса:" << '\n';
	b[1].printParameters(fout);
}

void Correct_new_vector(Spacecraft& Sputnik)
{
	Sputnik.from_Dec_to_Kepl();
	//Sputnik.RAAN = 0 * Pi / 180;
	//Sputnik.from_Kepl_Moon_to_Dec_J2000();
	Sputnik.omega = 90 * Pi / 180;
	Sputnik.from_Kepl_Moon_to_Dec_J2000();
}

void Impulses_to_hundred(double h0, double& impulses_sum)
{
	double V0 = sqrt(mu_Moon / (h0 + 1737.1));
	double V1 = sqrt(mu_Moon / 1837.1);
	double a_hohman = (100 + 1737.1 + 1737.1 + h0) / 2;
	double imp1 = V0 - sqrt(mu_Moon * (2 / (h0 + 1737.1) - 1 / a_hohman));
	double imp2 = sqrt(mu_Moon * (2 / (100 + 1737.1) - 1 / a_hohman)) - V1;
	impulses_sum = imp1 + imp2;
}

std::vector <Spacecraft> Integration_to_height_hundred(Spacecraft Sputnik, double h)
{
	// н.у. у Спутника уже должны быть прописаны нужные
	std::vector <Spacecraft> SC_Condition;
	SC_Condition.push_back(Sputnik);

	double height = Sputnik.r_KA - 6371;
	while (height >= 100)
	{
		Sputnik.Runge_Kutta(h);
		height = Sputnik.r_KA - 6371;
		SC_Condition.push_back(Sputnik);
	}

	SC_Condition.pop_back(); //удалил последний элемент массива, так как там спутник уже ниже 100 ушёл
	double h_new = h;

	Sputnik.Runge_Kutta(-2 * h_new);
	h_new /= 2;
	do
	{
		do
		{
			Sputnik.Runge_Kutta(h_new);
			height = Sputnik.r_KA - 6371;
		} while (height > 100);
		Sputnik.Runge_Kutta(-2 * h_new);
		height = Sputnik.r_KA - 6371;
		h_new /= 4;
	} while (abs(height - 100) > pow(10, -5) and (h_new > 0.05));
	SC_Condition.push_back(Sputnik);
	return SC_Condition;

}

double Precise_impulse_Dima(Spacecraft Sputnik_Dima, double target_trajectory_slope, double h)
{
	std::vector <Spacecraft> SC_Condition;
	double thrust_impulse = 0, trajectory_slope = 0;
	double X0 = 0, Y0 = 0, X1 = 0, Y1 = 0, X2 = 0, Y2 = 0;

	thrust_impulse = -0.12;
	Sputnik_Dima.Collinear_impulse(thrust_impulse);
	Sputnik_Dima.from_Dec_to_Kepl();
	SC_Condition = Integration_to_height_hundred(Sputnik_Dima, h);
	Sputnik_Dima.Collinear_impulse(-thrust_impulse); // возвращаю спутник в начальное положение
	Sputnik_Dima.from_Dec_to_Kepl();
	trajectory_slope = SC_Condition.back().trajectory_slope;
	X1 = thrust_impulse;
	Y1 = trajectory_slope - target_trajectory_slope;

	thrust_impulse -= 0.004;
	Sputnik_Dima.Collinear_impulse(thrust_impulse);
	Sputnik_Dima.from_Dec_to_Kepl();
	SC_Condition = Integration_to_height_hundred(Sputnik_Dima, h);
	Sputnik_Dima.Collinear_impulse(-thrust_impulse); // возвращаю спутник в начальное положение
	Sputnik_Dima.from_Dec_to_Kepl();
	trajectory_slope = SC_Condition.back().trajectory_slope;
	X2 = thrust_impulse;
	Y2 = trajectory_slope - target_trajectory_slope;
	do
	{
		X0 = Y1 / (Y1 - Y2) * (X2 - X1) + X1;
		thrust_impulse = X0;
		Sputnik_Dima.Collinear_impulse(thrust_impulse);
		Sputnik_Dima.from_Dec_to_Kepl();
		SC_Condition = Integration_to_height_hundred(Sputnik_Dima, h);
		Sputnik_Dima.Collinear_impulse(-thrust_impulse); // возвращаю спутник в начальное положение
		Sputnik_Dima.from_Dec_to_Kepl();
		trajectory_slope = SC_Condition.back().trajectory_slope;
		Y0 = trajectory_slope - target_trajectory_slope;
		if (abs(Y1) > abs(Y2))
		{
			X1 = X0;
			Y1 = Y0;
		}
		else
		{
			X2 = X0;
			Y2 = X0;
		}
	} while (abs(Y0) > pow(10, -4));
	return thrust_impulse;
}

int main()
{
	furnsh_c("naif0012.tls.pc.txt"); // подключил файл с leapseconds
	furnsh_c("de430.bsp"); // подключил эфемериды
	furnsh_c("APOPHIS_2020_2029_2.bsp"); // эфемериды Апофиса
	
	std::string str_approach_date = "2021-Mar-21 12:00:00"; // std::string str_approach_date = "2021-Feb-01 12:00:00"; // планируемая дата попадания в Луну
	Date_time_ms approach_date;
	approach_date.set_date_time(str_approach_date, 0);


	/// ДЕЛАЮ ДИМЕ >>
	/*Spacecraft Sputnik_Dima("Earth", "2021-Feb-01 00:00:00", 000);
	Sputnik_Dima.set_earth_compression(0);
	Sputnik_Dima.inc = 51.6*Pi/180;
	Sputnik_Dima.from_Kepl_to_Dec();
	Sputnik_Dima.from_Dec_to_Kepl();

	//double thrust_imp = Precise_impulse_Dima(Sputnik_Dima, 1.5 * Pi / 180, 1);
	Sputnik_Dima.Collinear_impulse(-0.12);
	Sputnik_Dima.from_Dec_to_Kepl();
	std::vector <Spacecraft> Sputnik_Dima_Par = Integration(Sputnik_Dima, 3600, 1);
	Write_to_Files(Sputnik_Dima_Par, 2);*/
	///ДЕЛАЮ ДИМЕ ^^


	/////////////////////////////////////////// // для начала интегрирования на NRHO 9:2 в перигее
	/*Spacecraft Sputnik("Moon", "2021-Feb-01 00:58:05", 657);
	std::vector <double> New_vector(6);
	New_vector[0] = -703.474;
	New_vector[1] = 899.338;
	New_vector[2] = 2926.11;
	New_vector[3] = -0.429074;
	New_vector[4] = 1.83426;
	New_vector[5] = -0.667872;
	Sputnik.set_parameters(New_vector);
	Correct_new_vector(Sputnik);

	Sputnik.set_earth_compression(0);
	Sputnik.set_moon_compression(0);
	Sputnik.set_sun_pressure(0);
	std::vector <bool> Planet_pertur = { 1,1,0,0,0,0 };
	Sputnik.set_planet_perturbations(Planet_pertur);


	double Synodic_resonance = 4.5, NRHO_r_peri = 3141;
	double NRHO_Orbit_Altitude_peri, NRHO_Orbit_Altitude_apo;
	double NRHO_r_apo, NRHO_Orbit_period, NRHO_V_peri, NRHO_V_apo;

	NRHO_parameters(Synodic_resonance, NRHO_r_peri, NRHO_r_apo, NRHO_Orbit_Altitude_peri, NRHO_Orbit_Altitude_apo, NRHO_Orbit_period, NRHO_V_peri, NRHO_V_apo);
	double Precise_impulse = Precise_NRHO_impulse_period(Sputnik, Moon_synodic_period / Synodic_resonance * 86400, NRHO_V_peri, 20); // импульс по периоду
	//double Precise_impulse = Precise_NRHO_impulse_aposelene(Sputnik, 71000, NRHO_V_peri, 20); // импульс по апоселению

	//std::vector <Spacecraft> Sputnik_Moon = NRHO_flight(Sputnik_Par.back(), NRHO_V_peri, NRHO_Orbit_period, 2, 20);
	//std::vector <Spacecraft> Sputnik_Moon = NRHO_flight(Sputnik, NRHO_V_peri, NRHO_Orbit_period, 2, 20);
	std::vector <Spacecraft> Sputnik_Moon = NRHO_flight(Sputnik, Sputnik.V_KA + Precise_impulse, NRHO_Orbit_period, 2, 20);
	Write_to_Files_Moon_orbit(Sputnik_Moon, 20);
	*/
	///////////////////////////////////////////

	/////////////////////////////////////////// // для начала интегрирования на NRHO 4:1 в перигее
	/*Spacecraft Sputnik("Moon", "2021-Feb-01 00:57:57", 593);
	std::vector <double> New_vector(6);
	New_vector[0] = -1033.59;
	New_vector[1] = 659.881;
	New_vector[2] = 5464.09;
	New_vector[3] = -0.413681;
	New_vector[4] = 1.54586;
	New_vector[5] = -0.265661;
	Sputnik.set_parameters(New_vector);
	Correct_new_vector(Sputnik);

	Sputnik.set_earth_compression(0);
	Sputnik.set_moon_compression(0);
	Sputnik.set_sun_pressure(0);
	std::vector <bool> Planet_pertur = { 0,1,0,0,0,0 };
	Sputnik.set_planet_perturbations(Planet_pertur);


	double Synodic_resonance = 4, NRHO_r_peri = 5600;
	double NRHO_Orbit_Altitude_peri, NRHO_Orbit_Altitude_apo;
	double NRHO_r_apo, NRHO_Orbit_period, NRHO_V_peri, NRHO_V_apo;

	NRHO_parameters(Synodic_resonance, NRHO_r_peri, NRHO_r_apo, NRHO_Orbit_Altitude_peri, NRHO_Orbit_Altitude_apo, NRHO_Orbit_period, NRHO_V_peri, NRHO_V_apo);
	//double Precise_impulse = Precise_NRHO_impulse_period(Sputnik, Moon_synodic_period / Synodic_resonance * 86400, NRHO_V_peri, 20); // импульс по периоду
	//double Precise_impulse = Precise_NRHO_impulse_aposelene(Sputnik, 71000, NRHO_V_peri, 20); // импульс по апоселению

	//std::vector <Spacecraft> Sputnik_Moon = NRHO_flight(Sputnik_Par.back(), NRHO_V_peri, NRHO_Orbit_period, 2, 20);
	std::vector <Spacecraft> Sputnik_Moon = NRHO_flight(Sputnik, NRHO_V_peri, NRHO_Orbit_period, 2, 20);
	//std::vector <Spacecraft> Sputnik_Moon = NRHO_flight(Sputnik, Sputnik.V_KA + Precise_impulse, NRHO_Orbit_period, 2, 20);
	Write_to_Files_Moon_orbit(Sputnik_Moon, 20);*/

	///////////////////////////////////////////


	/////////////////////////////////////////// // для расчёта перелётов методом Ньютона на круговые
	/*Spacecraft Sputnik("Earth", "2021-Feb-01 00:58:05", 657);
	double Circular_Orbit_Altitude = 100;

	Sputnik.set_earth_compression(1);
	Sputnik.set_moon_compression(1);
	Sputnik.set_sun_pressure(1);
	std::vector <bool> Planet_pertur = { 1,1,1,0,0,0 };
	Sputnik.set_planet_perturbations(Planet_pertur);


	//set_parameters_Moon_flight(Sputnik, 3.5, approach_date);
	Newtons_method_coords(Sputnik, 3.5, approach_date, 5, 40, 40); //первый этап задачи попадание в центр Луны
	Newtons_method_b_plane(Sputnik, 3.5, 5, 40, 40, Circular_Orbit_Altitude); //второй этап задачи попадание в параметры Картинной плоскости
	Newtons_method_final_step(Sputnik, 3.5, 5, 40, 40, Circular_Orbit_Altitude); //третий этап задачи попадание на оконч орбиту

	std::vector <Spacecraft> Sputnik_Par = Integration_to_Moon(Sputnik, 3.5, 5, 40, 40);
	Write_to_Files(Sputnik_Par, 20);
	std::vector <Spacecraft> Sputnik_Moon = Orbiting_Moon(Sputnik_Par.back(), Circular_Orbit_Altitude, 42000, 40);
	Write_to_Files_Moon_orbit(Sputnik_Moon, 20);
	Write_to_Files_initial_cond(Sputnik_Par, Sputnik_Moon);*/
	///////////////////////////////////////////

	/////////////////////////////////////////// // для расчёта перелётов методом Ньютона на NRHO 4:1
	/*Spacecraft Sputnik("Earth", "2021-Feb-01 00:58:05", 657);
	Sputnik.set_earth_compression(1);
	Sputnik.set_moon_compression(0);
	Sputnik.set_sun_pressure(1);
	std::vector <bool> Planet_pertur = { 1,1,1,0,0,0 };
	Sputnik.set_planet_perturbations(Planet_pertur);

	double Synodic_resonance = 4, NRHO_r_peri = 5600;
	double NRHO_Orbit_Altitude_peri, NRHO_Orbit_Altitude_apo;
	double NRHO_r_apo, NRHO_Orbit_period, NRHO_V_peri, NRHO_V_apo;
	NRHO_parameters(Synodic_resonance, NRHO_r_peri, NRHO_r_apo, NRHO_Orbit_Altitude_peri, NRHO_Orbit_Altitude_apo, NRHO_Orbit_period, NRHO_V_peri, NRHO_V_apo);

	//set_parameters_Moon_flight(Sputnik, 3.5, approach_date);
	// изменил время перелёта с 3.5 на 5.33
	Newtons_method_coords(Sputnik, 4, approach_date, 5, 40, 40); //первый этап задачи попадание в центр Луны
	Newtons_method_b_plane(Sputnik, 4, 5, 40, 40, NRHO_Orbit_Altitude_apo); //второй этап задачи попадание в параметры Картинной плоскости
	Newtons_method_final_step(Sputnik, 4, 5, 40, 40, NRHO_Orbit_Altitude_apo); //третий этап задачи попадание на оконч орбиту

	std::vector <Spacecraft> Sputnik_Par = Integration_to_Moon(Sputnik, 3.5, 5, 40, 40);
	Write_to_Files(Sputnik_Par, 20);

	Sputnik_Par.back().set_earth_compression(0); // подлетели к Луне, отключил все возмущения
	Sputnik_Par.back().set_moon_compression(0);
	Sputnik_Par.back().set_sun_pressure(0);
	Planet_pertur = { 0,1,0,0,0,0 };
	Sputnik_Par.back().set_planet_perturbations(Planet_pertur);

	std::vector <Spacecraft> Sputnik_Moon = NRHO_flight(Sputnik_Par.back(), NRHO_V_apo, NRHO_Orbit_period, 2, 20);
	Write_to_Files_Moon_orbit(Sputnik_Moon, 20);
	Write_to_Files_initial_cond(Sputnik_Par, Sputnik_Moon);
	*/
	///////////////////////////////////////////

	/////////////////////////////////////////// // перелёт от Земли к Луне уже по заданному
	/*Spacecraft Sputnik("Earth", "2021-Jan-28 12:58:35", 980);
	std::vector <double> New_vector(6);
	New_vector[0] = 6562.69;
	New_vector[1] = 7.9545;
	New_vector[2] = -329.747;
	New_vector[3] = 0.42469;
	New_vector[4] = 6.73611;
	New_vector[5] = 8.59272;
	Sputnik.set_parameters(New_vector);

	Sputnik.set_earth_compression(0);
	Sputnik.set_moon_compression(0);
	Sputnik.set_sun_pressure(0);
	std::vector <bool> Planet_pertur = { 1,0,0,0,0,0 };
	Sputnik.set_planet_perturbations(Planet_pertur);

	std::vector <Spacecraft> Sputnik_Par = Integration(Sputnik, 302400, 40);
	Write_to_Files(Sputnik_Par, 20);*/
	///////////////////////////////////////////

	////////////////////////////////////////// ПРОВЕРКА для Dr. Alfriend
	Spacecraft Sputnik("Moon", "2021-Feb-01 00:57:57", 593);
	std::vector <double> New_vector(6);
	New_vector[0] = 10000;
	New_vector[1] = 0;
	New_vector[2] = 0;
	New_vector[3] = 0;
	New_vector[4] = sqrt(mu_Moon/10000);
	New_vector[5] = 0;
	Sputnik.set_parameters(New_vector);

	std::vector <double> New_j2000_crds = Sputnik.SVSK_to_J2000_6x6();
	double r_j2000 = sqrt(New_j2000_crds[0] * New_j2000_crds[0] + New_j2000_crds[1] * New_j2000_crds[1] + New_j2000_crds[2] * New_j2000_crds[2]);
	double v_j2000 = sqrt(New_j2000_crds[3] * New_j2000_crds[3] + New_j2000_crds[4] * New_j2000_crds[4] + New_j2000_crds[5] * New_j2000_crds[5]);
	Sputnik.a = 10000;
	Sputnik.tau = 0;
	Sputnik.ecc = 0;
	Sputnik.omega = 0 * Pi / 180;
	Sputnik.RAAN = 270*Pi/180;
	Sputnik.inc = 6.68*Pi/180;

	Sputnik.from_Kepl_Moon_to_Dec_J2000();
	Sputnik.from_Dec_to_Kepl();

	Sputnik.set_earth_compression(0);
	Sputnik.set_moon_compression(0);
	Sputnik.set_sun_pressure(0);
	std::vector <bool> Planet_pertur = { 1,1,0,0,0,0 };
	Sputnik.set_planet_perturbations(Planet_pertur);

	std::vector <Spacecraft> Sputnik_Par = Integration(Sputnik, 2375065, 120);
	Write_to_Files(Sputnik_Par, 5);
	Write_to_Files_Moon_orbit(Sputnik_Par, 5);
	//////////////////////////////////////////
}

