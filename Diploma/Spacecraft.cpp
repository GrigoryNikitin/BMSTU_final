#include "pch.h"
#include "Spacecraft.h"

#pragma comment(lib, "../libs/cspice/lib/cspice.lib")
extern "C" {
#include "../libs/cspice/include/SpiceCK.h" // must be before SpiceZdf.h and SpiceZpr.h
#include "../libs/cspice/include/SpiceZdf.h"
#include "../libs/cspice/include/SpiceZpr.h"
}

Spacecraft::Spacecraft()
{
	X_KA = X0;
	Y_KA = Y0;
	Z_KA = Z0;
	Vx_KA = Vx0;
	Vy_KA = Vy0;
	Vz_KA = Vz0;
	r_KA = sqrt(X_KA * X_KA + Y_KA * Y_KA + Z_KA * Z_KA);
	V_KA = sqrt(Vx_KA * Vx_KA + Vy_KA * Vy_KA + Vz_KA * Vz_KA);
	//Date_time.Date = { 0 }; // по умолчанию ставится 0, надо функцией set_Date менять сразу
	//Date_time.ms = 0;
	//пока что варик только при начальных у Земли, поэтому сразу мю Земли все такое
	earth_compression_flag = 0; // по умолчанию ставлю 0
	moon_compression_flag = 0; // по умолчанию ставлю 0
	sun_pressure_flag = 0; // по умолчанию ставлю 0
	Grav_parameter = mu_Earth;
	main_center_flag = "Earth";
	time = 0;
	from_Dec_to_Kepl();
	if (a > 0) { Period = 2 * Pi * sqrt(a * a * a / Grav_parameter); }
	else { Period = -1; }

	planet_perturbations.resize(6);
	for (size_t i = 0; i < 6; i++)
	{
		planet_perturbations[i] = 0; // по умолчанию все 0
	}
	ksi_bplane.resize(3);
	etta_bplane.resize(3);
	zeta_bplane.resize(3);
}

Spacecraft::Spacecraft(std::string main_center, std::string data, int ms)
{
	Date_time.set_date_time(data, ms);
	X_KA = X0;
	Y_KA = Y0;
	Z_KA = Z0;
	Vx_KA = Vx0;
	Vy_KA = Vy0;
	Vz_KA = Vz0;
	r_KA = sqrt(X_KA * X_KA + Y_KA * Y_KA + Z_KA * Z_KA);
	V_KA = sqrt(Vx_KA * Vx_KA + Vy_KA * Vy_KA + Vz_KA * Vz_KA);
	//Date_time.Date = { 0 }; // по умолчанию ставится 0, надо функцией set_Date менять сразу
	//Date_time.ms = 0;
	earth_compression_flag = 0; // по умолчанию ставлю 0
	moon_compression_flag = 0; // по умолчанию ставлю 0
	sun_pressure_flag = 0; // по умолчанию ставлю 0
	if (main_center == "Earth")
	{
		Grav_parameter = mu_Earth;
		main_center_flag = "Earth";
	}
	if (main_center == "Moon")
	{
		Grav_parameter = mu_Moon;
		main_center_flag = "Moon";
	}
	if (main_center == "Sun")
	{
		Grav_parameter = mu_Sun;
		main_center_flag = "Sun";
	}
	if (main_center == "Apophis")
	{
		Grav_parameter = mu_Apophis;
		main_center_flag = "Apophis";
	}
	time = 0;
	from_Dec_to_Kepl();
	if (a > 0) { Period = 2 * Pi * sqrt(a * a * a / Grav_parameter); }
	else { Period = -1; }

	planet_perturbations.resize(6);
	for (size_t i = 0; i < 6; i++)
	{
		planet_perturbations[i] = 0; // по умолчанию все 0
	}
	ksi_bplane.resize(3);
	etta_bplane.resize(3);
	zeta_bplane.resize(3);
}

Spacecraft::~Spacecraft() {}

void Spacecraft::printParameters(std::ofstream& fout)
{
	fout << Date_time.get_string() << '\t' << time << '\t' << X_KA << '\t' << Y_KA << '\t' << Z_KA << '\t' << r_KA << '\t' << distance_to_Moon() << '\t' << Vx_KA << '\t' <<
		Vy_KA << '\t' << Vz_KA << '\t' << V_KA << '\t' << main_center_flag << '\t' << trajectory_slope * 180 / Pi << '\n';

}
void Spacecraft::printKepler(std::ofstream& fout)
{
	fout << Date_time.get_string() << '\t' << time << '\t' << a << '\t' << ecc << '\t' << inc * toDeg << '\t' << RAAN * toDeg << '\t' << omega * toDeg << '\t' << tau << '\t' << theta * toDeg << '\t' << main_center_flag << '\n';
}
void Spacecraft::printBPlane(std::ofstream& fout)
{
	/*char Date_char[32];
	strftime(Date_char, 32, "%Y-%b-%d %H:%M:%S", &Date); // из tm в массив из char	*/
	fout << Date_time.get_string() << '\t' << time << '\t' << V_inf << '\t' << B_ksi << '\t' << B_etta << '\n';
}
void Spacecraft::printSVSK(std::ofstream& fout)
{
	std::vector <double> SVSK_crds(6);
	SVSK_crds = J2000_to_SVSK_6x6();
	fout << Date_time.get_string() << '\t' << time << '\t' << SVSK_crds[0] << '\t' << SVSK_crds[1] << '\t' << SVSK_crds[2] << '\t' << distance_to_Moon() << '\t' <<
		SVSK_crds[3] << '\t' << SVSK_crds[4] << '\t' << SVSK_crds[5] << '\t' << SVSK_inclination() * 180 / Pi << '\n';
}

void Spacecraft::printGCS(std::ofstream& fout)
{
	std::vector <double> GCS_crds(6);
	GCS_crds = J2000_to_GCS_6x6();
	fout << Date_time.get_string() << '\t' << time << '\t' << GCS_crds[0] << '\t' << GCS_crds[1] << '\t' << GCS_crds[2] << '\t' << r_KA <<
		'\t' << GCS_crds[3] << '\t' << GCS_crds[4] << '\t' << GCS_crds[5] << '\t' << V_KA << '\t' << trajectory_slope * 180 / Pi << '\n';
}

void Spacecraft::printParameters_relatively_to_Moon(std::ofstream& fout)
{
	if (main_center_flag == "Moon")
	{
		fout << Date_time.get_string() << '\t' << time << '\t' << X_KA << '\t' << Y_KA << '\t' << Z_KA << '\t' << r_KA << '\t' << distance_to_Moon() << '\t' << Vx_KA << '\t' <<
			Vy_KA << '\t' << Vz_KA << '\t' << V_KA << '\t' << main_center_flag << '\n';
	}
	if (main_center_flag == "Earth")
	{
		std::vector <double> Moon_state(6); // вектор параметров Луны
		double X_M, Y_M, Z_M, Vx_M, Vy_M, Vz_M;
		double X_M0, Y_M0, Z_M0, Vx_M0, Vy_M0, Vz_M0;

		Moon_state = get_planet_parameters(Date_time.get_string().c_str(), "MOON", "J2000", "EARTH");
		X_M = Moon_state[0]; // вектор Земля - Луна
		Y_M = Moon_state[1];
		Z_M = Moon_state[2];
		Vx_M = Moon_state[3]; // вектор Земля - Луна
		Vy_M = Moon_state[4];
		Vz_M = Moon_state[5];
		X_M0 = X_KA - X_M; // считаем вектор от Луны на КА
		Y_M0 = Y_KA - Y_M;
		Z_M0 = Z_KA - Z_M;
		Vx_M0 = Vx_KA - Vx_M; // считаем вектор от Луны на КА
		Vy_M0 = Vy_KA - Vy_M;
		Vz_M0 = Vz_KA - Vz_M;
		fout << Date_time.get_string() << '\t' << time << '\t' << X_M0 << '\t' << Y_M0 << '\t' << Z_M0 << '\t' << r_KA << '\t' << distance_to_Moon() << '\t' << Vx_M0 << '\t' <<
			Vy_M0 << '\t' << Vz_M0 << '\t' << sqrt(pow(Vx_M0, 2) + pow(Vy_M0, 2) + pow(Vz_M0, 2)) << '\t' << main_center_flag << '\n';
	}

}
void Spacecraft::printKepler_relatively_to_Moon(std::ofstream& fout)
{
	if (main_center_flag == "Moon")
	{
		fout << Date_time.get_string() << '\t' << time << '\t' << a << '\t' << ecc << '\t' << inc * toDeg << '\t' << RAAN * toDeg << '\t' << omega * toDeg << '\t' << tau << '\t' << theta * toDeg << '\t' << main_center_flag << '\n';
	}
	if (main_center_flag == "Earth")
	{
		double a_, ecc_, inc_, RAAN_, omega_, tau_, theta_;
		from_Dec_to_Kepl("Moon", a_, ecc_, inc_, RAAN_, omega_, tau_, theta_);
		fout << Date_time.get_string() << '\t' << time << '\t' << a_ << '\t' << ecc_ << '\t' << inc_ * toDeg << '\t' << RAAN_ * toDeg << '\t' << omega_ * toDeg << '\t' << tau_ << '\t' << theta_ * toDeg << '\t' << main_center_flag << '\n';
	}
}
void Spacecraft::printSVSK_relatively_to_Moon(std::ofstream& fout)
{
	std::vector <double> SVSK_crds(6);
	SVSK_crds = J2000_to_SVSK_6x6();
	fout << Date_time.get_string() << '\t' << time << '\t' << SVSK_crds[0] << '\t' << SVSK_crds[1] << '\t' << SVSK_crds[2] << '\t' << distance_to_Moon() << '\t' <<
		SVSK_crds[3] << '\t' << SVSK_crds[4] << '\t' << SVSK_crds[5] << '\t' << SVSK_inclination() * 180 / Pi << '\n';
}

void Spacecraft::printSISK_relatively_to_Moon(std::ofstream& fout)
{
	std::vector <double> SISK_crds(6);
	SISK_crds = J2000_to_SISK_6x6();
	fout << Date_time.get_string() << '\t' << time << '\t' << SISK_crds[0] << '\t' << SISK_crds[1] << '\t' << SISK_crds[2] << '\t' << distance_to_Moon() << '\t' <<
		SISK_crds[3] << '\t' << SISK_crds[4] << '\t' << SISK_crds[5] << '\t' << SVSK_inclination() * 180 / Pi << '\n';
}

std::vector <double> Spacecraft::get_parameters()
{
	std::vector <double> a;
	a.resize(6);
	a[0] = X_KA;
	a[1] = Y_KA;
	a[2] = Z_KA;
	a[3] = Vx_KA;
	a[4] = Vy_KA;
	a[5] = Vz_KA;
	return a;
}

std::vector <double> Spacecraft::get_b_plane_parameters()
{
	std::vector <double> a;
	a.resize(2);
	a[0] = B_ksi;
	a[1] = B_etta;
	return a;
}

void Spacecraft::set_parameters(std::vector <double> Parameters)
{
	X_KA = Parameters[0];
	Y_KA = Parameters[1];
	Z_KA = Parameters[2];
	Vx_KA = Parameters[3];
	Vy_KA = Parameters[4];
	Vz_KA = Parameters[5];
	r_KA = sqrt(X_KA * X_KA + Y_KA * Y_KA + Z_KA * Z_KA);
	V_KA = sqrt(Vx_KA * Vx_KA + Vy_KA * Vy_KA + Vz_KA * Vz_KA);
	from_Dec_to_Kepl();
	if (a > 0) { Period = 2 * Pi * sqrt(a * a * a / Grav_parameter); }
	else { Period = -1; }
}

void Spacecraft::set_Date(std::string a, int ms)
{
	Date_time.set_date_time(a, ms);
}

void Spacecraft::set_earth_compression(int a)
{
	if (a == 0) { earth_compression_flag = 0; }
	else { earth_compression_flag = 1; }
}

void Spacecraft::set_moon_compression(int a)
{
	if (a == 0) { moon_compression_flag = 0; }
	else { moon_compression_flag = 1; }
}

void Spacecraft::set_sun_pressure(int a)
{
	if (a == 0) { sun_pressure_flag = 0; }
	else { sun_pressure_flag = 1; }
}

void Spacecraft::set_planet_perturbations(std::vector <bool> a)
{
	planet_perturbations = a;
}

void Spacecraft::from_Dec_to_Kepl()
{
	double X = 0, Y = 0, Z = 0, Vx = 0, Vy = 0, Vz = 0;
	r_KA = sqrt(X_KA * X_KA + Y_KA * Y_KA + Z_KA * Z_KA);
	if (main_center_flag == "Earth")
	{
		X = X_KA;
		Y = Y_KA;
		Z = Z_KA;
		Vx = Vx_KA;
		Vy = Vy_KA;
		Vz = Vz_KA;
	}
	if (main_center_flag == "Moon")
	{
		std::vector <double> SVSK_crds(6); // теперь считаю по SISK, а не SVSK (ИЛИ НЕТ!)
		SVSK_crds = J2000_to_SVSK_6x6();
		X = SVSK_crds[0];
		Y = SVSK_crds[1];
		Z = SVSK_crds[2];
		Vx = SVSK_crds[3];
		Vy = SVSK_crds[4];
		Vz = SVSK_crds[5];
	}

	double C_x = Y * Vz - Z * Vy;
	double C_y = -X * Vz + Z * Vx;
	double C_z = X * Vy - Y * Vx;
	double C = sqrt(C_x * C_x + C_y * C_y + C_z * C_z);
	trajectory_slope = Pi / 2 - asin(C / (r_KA * V_KA));
	inc = acos(C_z / C);    //inc КА
	if (inc == 0) { inc = 0.0000000001; }

	if (inc < 0.000000001) { RAAN = 0.0000000001; }      // RAAN КА
	else { RAAN = atan2(C_x, -C_y); }
	if (RAAN < 0) { RAAN += 2 * Pi; }

	double f_x = -Grav_parameter * X / r_KA - Vz * (Z * Vx - X * Vz) + Vy * (X * Vy - Y * Vx);
	double f_y = -Grav_parameter * Y / r_KA - Vx * (X * Vy - Y * Vx) + Vz * (Y * Vz - Z * Vy);
	double f_z = -Grav_parameter * Z / r_KA - Vy * (Y * Vz - Z * Vy) + Vx * (Z * Vx - X * Vz);
	double Lap = sqrt(f_x * f_x + f_y * f_y + f_z * f_z);

	p = C * C / Grav_parameter; //параметр орбиты
	ecc = Lap / Grav_parameter;         // Эксцентриситет орбиты КА
	a = p / (1 - ecc * ecc);           //Большая полуось орбиты КА
	if (ecc <= pow(10, -7))
	{
		omega = 0;
	}
	else
	{
		omega = atan2(f_z / (Lap * sin(inc)), (f_x * cos(RAAN) + f_y * sin(RAAN)) / Lap); //Аргумент перицентра КА
		if (omega < 0) { omega += 2 * Pi; }
	}

	u = atan2(Z / (r_KA * sin(inc)), (X * cos(RAAN) + Y * sin(RAAN)) / r_KA);   //Аргумент широты КА
	if (u < 0) { u += 2 * Pi; }

	if (ecc <= pow(10, -7))
	{
		theta = u - omega;
		if (theta < 0) { theta += 2 * Pi; }
	}
	else
	{
		theta = atan2((X * Vx + Y * Vy + Z * Vz) / (r_KA * sqrt(Grav_parameter / p) * ecc), (p / r_KA - 1) / ecc);
		if (theta < 0) { theta += 2 * Pi; }
	}
	if (ecc < 1)
	{
		EA = 2 * atan(sqrt((1 - ecc) / (1 + ecc)) * tan(0.5 * theta));
		if (EA < 0) { EA += 2 * Pi; }
		MA = EA - ecc * sin(EA);
	}
	else
	{
		H = atanh(sqrt(ecc * ecc - 1) * sin(theta) / (ecc + cos(theta)));
		MA = ecc * sinh(H) - H;
	}
	/*if (time == 0) { tau = 0; }
	else { tau = time - sqrt(abs(a*a*a) / Grav_parameter)*MA; }*/
	tau = time - sqrt(abs(a * a * a) / Grav_parameter) * MA;
}

void Spacecraft::from_Dec_to_Kepl(std::string chosen_center_flag, double& a, double& ecc, double& inc, double& RAAN, double& omega, double& tau, double& theta)
{
	double X = 0, Y = 0, Z = 0, Vx = 0, Vy = 0, Vz = 0, r_KA = 0, Grav_parameter = 0;
	if (chosen_center_flag == "Earth")
	{
		Grav_parameter = mu_Earth;
		if (main_center_flag == "Earth")
		{
			X = X_KA;
			Y = Y_KA;
			Z = Z_KA;
			Vx = Vx_KA;
			Vy = Vy_KA;
			Vz = Vz_KA;
		}
		if (main_center_flag == "Moon")
		{
			std::vector <double> Earth_state(6); // вектор параметров Земли относительно Луны
			double X_E, Y_E, Z_E, Vx_E, Vy_E, Vz_E;
			double X_E0, Y_E0, Z_E0, Vx_E0, Vy_E0, Vz_E0;

			Earth_state = get_planet_parameters(Date_time.get_string().c_str(), "EARTH", "J2000", "MOON");
			X_E = Earth_state[0]; // вектор Луна - Земля
			Y_E = Earth_state[1];
			Z_E = Earth_state[2];
			Vx_E = Earth_state[3]; // вектор Луна - Земля
			Vy_E = Earth_state[4];
			Vz_E = Earth_state[5];

			X = X_KA - X_E; // считаем вектор от Земли на КА
			Y = Y_KA - Y_E;
			Z = Z_KA - Z_E;
			Vx = Vx_KA - Vx_E; // считаем вектор от Земли на КА
			Vy = Vy_KA - Vy_E;
			Vz = Vz_KA - Vz_E;
		}
	}
	if (chosen_center_flag == "Moon")
	{
		Grav_parameter = mu_Moon;
		std::vector <double> SVSK_crds(6); // теперь считаю по SISK, а не SVSK
		SVSK_crds = J2000_to_SVSK_6x6();
		X = SVSK_crds[0];
		Y = SVSK_crds[1];
		Z = SVSK_crds[2];
		Vx = SVSK_crds[3];
		Vy = SVSK_crds[4];
		Vz = SVSK_crds[5];
	}
	r_KA = sqrt(X * X + Y * Y + Z * Z);

	double C_x = Y * Vz - Z * Vy;
	double C_y = -X * Vz + Z * Vx;
	double C_z = X * Vy - Y * Vx;
	double C = sqrt(C_x * C_x + C_y * C_y + C_z * C_z);
	inc = acos(C_z / C);    //inc КА
	if (inc == 0) { inc = 0.0000000001; }

	if (inc < 0.000000001) { RAAN = 0.0000000001; }      // RAAN КА
	else { RAAN = atan2(C_x, -C_y); }
	if (RAAN < 0) { RAAN += 2 * Pi; }

	double f_x = -Grav_parameter * X / r_KA - Vz * (Z * Vx - X * Vz) + Vy * (X * Vy - Y * Vx);
	double f_y = -Grav_parameter * Y / r_KA - Vx * (X * Vy - Y * Vx) + Vz * (Y * Vz - Z * Vy);
	double f_z = -Grav_parameter * Z / r_KA - Vy * (Y * Vz - Z * Vy) + Vx * (Z * Vx - X * Vz);
	double Lap = sqrt(f_x * f_x + f_y * f_y + f_z * f_z);

	p = C * C / Grav_parameter; //параметр орбиты
	ecc = Lap / Grav_parameter;         // Эксцентриситет орбиты КА
	a = p / (1 - ecc * ecc);           //Большая полуось орбиты КА
	if (ecc <= pow(10, -7))
	{
		omega = 0;
	}
	else
	{
		omega = atan2(f_z / (Lap * sin(inc)), (f_x * cos(RAAN) + f_y * sin(RAAN)) / Lap); //Аргумент перицентра КА
		if (omega < 0) { omega += 2 * Pi; }
	}

	u = atan2(Z / (r_KA * sin(inc)), (X * cos(RAAN) + Y * sin(RAAN)) / r_KA);   //Аргумент широты КА
	if (u < 0) { u += 2 * Pi; }

	if (ecc <= pow(10, -7))
	{
		theta = u - omega;
		if (theta < 0) { theta += 2 * Pi; }
	}
	else
	{
		theta = atan2((X * Vx + Y * Vy + Z * Vz) / (r_KA * sqrt(Grav_parameter / p) * ecc), (p / r_KA - 1) / ecc);
		if (theta < 0) { theta += 2 * Pi; }
	}
	if (ecc < 1)
	{
		EA = 2 * atan(sqrt((1 - ecc) / (1 + ecc)) * tan(0.5 * theta));
		if (EA < 0) { EA += 2 * Pi; }
		MA = EA - ecc * sin(EA);
	}
	else
	{
		H = atanh(sqrt(ecc * ecc - 1) * sin(theta) / (ecc + cos(theta)));
		MA = ecc * sinh(H) - H;
	}
	/*if (time == 0) { tau = 0; }
	else { tau = time - sqrt(abs(a*a*a) / Grav_parameter)*MA; }*/
	tau = time - sqrt(abs(a * a * a) / Grav_parameter) * MA;
}

void Spacecraft::from_Kepl_to_Dec()
{
	double EA_, EA_next, H_, H_next;
	MA = sqrt(Grav_parameter / (abs(a * a * a))) * (time - tau);
	if (ecc < 1)
	{
		EA_next = MA + ecc * sin(MA);
		do
		{
			EA_ = EA_next;
			EA_next = EA_ - (EA_ - ecc * sin(EA_) - MA) / (1 - ecc * cos(EA_));
		} while (abs(EA_next - EA_) >= pow(10, -9));
		EA = EA_next;
		theta = atan2(sqrt(1 - ecc * ecc) * sin(EA) / (1 - ecc * cos(EA)), (cos(EA) - ecc) / (1 - ecc * cos(EA)));
		if (theta < 0) { theta += 2 * Pi; }
	}
	else
	{
		H_next = -MA + ecc * sinh(MA);
		do
		{
			H_ = H_next;
			H_next = H_ - (ecc * sinh(H_) - H_ - MA) / (ecc * cosh(H_) - 1);
		} while (abs(H_next - H_) >= pow(10, -9));
		H = H_next;
		theta = atan2(sqrt(ecc * ecc - 1) * sinh(H) / (ecc * cosh(H) - 1), (ecc - cosh(H)) / (ecc * cosh(H) - 1));
		if (theta < 0) { theta += 2 * Pi; }
	}

	p = a * (1 - ecc * ecc);
	double r_KA = p / (1 + ecc * cos(theta));
	double V_r = sqrt(Grav_parameter / p) * ecc * sin(theta);
	double V_m = sqrt(Grav_parameter / p) * (1 + ecc * cos(theta));
	u = omega + theta;
	X_KA = (cos(RAAN) * cos(u) - sin(RAAN) * sin(u) * cos(inc)) * r_KA;
	Y_KA = (sin(RAAN) * cos(u) + cos(RAAN) * sin(u) * cos(inc)) * r_KA;
	Z_KA = (sin(u) * sin(inc)) * r_KA;
	Vx_KA = (cos(RAAN) * cos(u) - sin(RAAN) * sin(u) * cos(inc)) * V_r + (-cos(RAAN) * sin(u) - sin(RAAN) * cos(u) * cos(inc)) * V_m;
	Vy_KA = (sin(RAAN) * cos(u) + cos(RAAN) * sin(u) * cos(inc)) * V_r + (-sin(RAAN) * sin(u) + cos(RAAN) * cos(u) * cos(inc)) * V_m;
	Vz_KA = (sin(u) * sin(inc)) * V_r + (cos(u) * sin(inc)) * V_m;
	this->r_KA = sqrt(X_KA * X_KA + Y_KA * Y_KA + Z_KA * Z_KA);
	this->V_KA = sqrt(Vx_KA * Vx_KA + Vy_KA * Vy_KA + Vz_KA * Vz_KA);
}

void Spacecraft::from_Kepl_Moon_to_Dec_J2000()
{
	double EA_, EA_next, H_, H_next;
	MA = sqrt(Grav_parameter / (abs(a * a * a))) * (time - tau);
	if (ecc < 1)
	{
		EA_next = MA + ecc * sin(MA);
		do
		{
			EA_ = EA_next;
			EA_next = EA_ - (EA_ - ecc * sin(EA_) - MA) / (1 - ecc * cos(EA_));
		} while (abs(EA_next - EA_) >= pow(10, -9));
		EA = EA_next;
		theta = atan2(sqrt(1 - ecc * ecc) * sin(EA) / (1 - ecc * cos(EA)), (cos(EA) - ecc) / (1 - ecc * cos(EA)));
		if (theta < 0) { theta += 2 * Pi; }
	}
	else
	{
		H_next = -MA + ecc * sinh(MA);
		do
		{
			H_ = H_next;
			H_next = H_ - (ecc * sinh(H_) - H_ - MA) / (ecc * cosh(H_) - 1);
		} while (abs(H_next - H_) >= pow(10, -9));
		H = H_next;
		theta = atan2(sqrt(ecc * ecc - 1) * sinh(H) / (ecc * cosh(H) - 1), (ecc - cosh(H)) / (ecc * cosh(H) - 1));
		if (theta < 0) { theta += 2 * Pi; }
	}

	p = a * (1 - ecc * ecc);
	double r_KA = p / (1 + ecc * cos(theta));
	double V_r = sqrt(Grav_parameter / p) * ecc * sin(theta);
	double V_m = sqrt(Grav_parameter / p) * (1 + ecc * cos(theta));
	u = omega + theta;
	X_KA = (cos(RAAN) * cos(u) - sin(RAAN) * sin(u) * cos(inc)) * r_KA;
	Y_KA = (sin(RAAN) * cos(u) + cos(RAAN) * sin(u) * cos(inc)) * r_KA;
	Z_KA = (sin(u) * sin(inc)) * r_KA;
	Vx_KA = (cos(RAAN) * cos(u) - sin(RAAN) * sin(u) * cos(inc)) * V_r + (-cos(RAAN) * sin(u) - sin(RAAN) * cos(u) * cos(inc)) * V_m;
	Vy_KA = (sin(RAAN) * cos(u) + cos(RAAN) * sin(u) * cos(inc)) * V_r + (-sin(RAAN) * sin(u) + cos(RAAN) * cos(u) * cos(inc)) * V_m;
	Vz_KA = (sin(u) * sin(inc)) * V_r + (cos(u) * sin(inc)) * V_m;
	// ПОЛУЧИЛ КООРДИНАТЫ В SVSK, ДАЛЬШЕ ПЕРЕВОЖУ ИХ В J2000
	std::vector <double> J2000_crds(6);
	J2000_crds = SVSK_to_J2000_6x6();
	X_KA = J2000_crds[0];
	Y_KA = J2000_crds[1];
	Z_KA = J2000_crds[2];
	Vx_KA = J2000_crds[3];
	Vy_KA = J2000_crds[4];
	Vz_KA = J2000_crds[5];

	this->r_KA = sqrt(X_KA * X_KA + Y_KA * Y_KA + Z_KA * Z_KA);
	this->V_KA = sqrt(Vx_KA * Vx_KA + Vy_KA * Vy_KA + Vz_KA * Vz_KA);
}

std::vector <double> Spacecraft::J2000_to_SISK_6x6()
{
	std::vector <double> SISK_coords(6);
	double delta_t, T, alpha, beta;
	double Omega_L;
	double e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13;
	double a, b;
	double X, Y, Z, Vx, Vy, Vz;

	if (main_center_flag == "Moon")
	{
		X = X_KA;
		Y = Y_KA;
		Z = Z_KA;
		Vx = Vx_KA;
		Vy = Vy_KA;
		Vz = Vz_KA;
	}
	if (main_center_flag == "Earth")
	{
		std::vector <double> Moon_state(6); // вектор параметров Луны
		double X_M, Y_M, Z_M, Vx_M, Vy_M, Vz_M;

		Moon_state = get_planet_parameters(Date_time.get_string().c_str(), "MOON", "J2000", "EARTH");
		X_M = Moon_state[0]; // вектор Земля - Луна
		Y_M = Moon_state[1];
		Z_M = Moon_state[2];
		Vx_M = Moon_state[3]; // вектор Земля - Луна
		Vy_M = Moon_state[4];
		Vz_M = Moon_state[5];

		X = X_KA - X_M; //пересчитываю координаты
		Y = Y_KA - Y_M;
		Z = Z_KA - Z_M;
		Vx = Vx_KA - Vx_M; // считаем вектор от Луны на КА
		Vy = Vy_KA - Vy_M;
		Vz = Vz_KA - Vz_M;
	}

	Omega_L = 13.17635815; // град/сутки
	str2et_c(Date_time.get_string().c_str(), &delta_t); // s, seconds past J2000 epoch
	delta_t /= (60 * 60 * 24); // Разница между датой и J2000 в сутках
	T = delta_t / 36525;

	e1 = (125.045 - 0.0529921 * delta_t) * Pi / 180;  //Сразу все в радианы перевожу
	e2 = (250.089 - 0.1059842 * delta_t) * Pi / 180;
	e3 = (260.008 + 13.0120009 * delta_t) * Pi / 180;
	e4 = (176.625 + 13.3407154 * delta_t) * Pi / 180;
	e5 = (357.529 + 0.9856003 * delta_t) * Pi / 180;
	e6 = (311.589 + 26.4057084 * delta_t) * Pi / 180;
	e7 = (134.963 + 13.0649930 * delta_t) * Pi / 180;
	e8 = (276.617 + 0.3287146 * delta_t) * Pi / 180;
	e9 = (34.226 + 1.7484877 * delta_t) * Pi / 180;
	e10 = (15.134 - 0.1589763 * delta_t) * Pi / 180;
	e11 = (119.743 + 0.0036096 * delta_t) * Pi / 180;
	e12 = (239.961 + 0.1643573 * delta_t) * Pi / 180;
	e13 = (25.053 + 12.9590088 * delta_t) * Pi / 180;

	alpha = 269.9949 + 0.0031 * T - 3.8787 * sin(e1) - 0.1204 * sin(e2) + 0.07 * sin(e3) -
		0.0172 * sin(e4) + 0.0072 * sin(e6) - 0.0052 * sin(e10) + 0.0043 * sin(e13); // в градусах
	beta = 66.5392 + 0.013 * T + 1.5419 * cos(e1) + 0.0239 * cos(e2) - 0.0278 * cos(e3) +
		0.0068 * cos(e4) - 0.0029 * cos(e6) + 0.0009 * cos(e7) + 0.0008 * cos(e10) - 0.0009 * cos(e13); // в градусах

	a = alpha * Pi / 180;
	b = beta * Pi / 180;
	// перевел их в радианы

	SISK_coords[0] = cos(Pi / 2 + a) * X + sin(Pi / 2 + a) * Y + 0 * Z;
	SISK_coords[1] = -sin(Pi / 2 + a) * cos(Pi / 2 - b) * X + cos(Pi / 2 + a) * cos(Pi / 2 - b) * Y + sin(Pi / 2 - b) * Z;
	SISK_coords[2] = sin(Pi / 2 + a) * sin(Pi / 2 - b) * X - cos(Pi / 2 + a) * sin(Pi / 2 - b) * Y + cos(Pi / 2 - b) * Z;

	SISK_coords[3] = cos(Pi / 2 + a) * Vx + sin(Pi / 2 + a) * Vy + 0 * Vz;
	SISK_coords[4] = -sin(Pi / 2 + a) * cos(Pi / 2 - b) * Vx + cos(Pi / 2 + a) * cos(Pi / 2 - b) * Vy + sin(Pi / 2 - b) * Vz;
	SISK_coords[5] = sin(Pi / 2 + a) * sin(Pi / 2 - b) * Vx - cos(Pi / 2 + a) * sin(Pi / 2 - b) * Vy + cos(Pi / 2 - b) * Vz;
	return SISK_coords;
}

std::vector <double> Spacecraft::J2000_to_GCS_6x6()
{
	std::vector <double> GCS_coords(6);
	double X, Y, Z, Vx, Vy, Vz;

	if (main_center_flag == "Earth")
	{
		X = X_KA;
		Y = Y_KA;
		Z = Z_KA;
		Vx = Vx_KA;
		Vy = Vy_KA;
		Vz = Vz_KA;
	}
	if (main_center_flag == "Moon")
	{
		std::vector <double> Earth_state(6); // вектор параметров Земли
		double X_E, Y_E, Z_E, Vx_E, Vy_E, Vz_E;

		Earth_state = get_planet_parameters(Date_time.get_string().c_str(), "EARTH", "J2000", "MOON");
		X_E = Earth_state[0]; // вектор Луна - Земля
		Y_E = Earth_state[1];
		Z_E = Earth_state[2];
		Vx_E = Earth_state[3]; // вектор Луна - Земля
		Vy_E = Earth_state[4];
		Vz_E = Earth_state[5];

		X = X_KA - X_E; //пересчитываю координаты
		Y = Y_KA - Y_E;
		Z = Z_KA - Z_E;
		Vx = Vx_KA - Vx_E; // считаем вектор от Луны на КА
		Vy = Vy_KA - Vy_E;
		Vz = Vz_KA - Vz_E;
	}

	double alpha = 0;
	alpha = alpfa0(Date_time.Date);  //альфа в часах
	alpha = alpha * 15;           //в градусах
	alpha = alpha * toRad;      // в радианах

	GCS_coords[0] = cos(alpha) * X + sin(alpha) * Y + 0 * Z;
	GCS_coords[1] = -sin(alpha) * X + cos(alpha) * Y + 0 * Z;
	GCS_coords[2] = 0 * X + 0 * Y + 1 * Z;

	GCS_coords[3] = cos(alpha) * Vx + sin(alpha) * Vy + 0 * Vz;
	GCS_coords[4] = -sin(alpha) * Vx + cos(alpha) * Vy + 0 * Vz;
	GCS_coords[5] = 0 * Vx + 0 * Vy + 1 * Vz;
	return GCS_coords;
}

matrix Spacecraft::J2000_to_SVSK_matrix_3x3(Date_time_ms Date)
{
	matrix A(3, 3);
	double delta_t, T, alpha, beta, omega_delta_t;
	double Omega_L;
	double e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13;
	double a, b, om;
	double X, Y, Z;

	if (main_center_flag == "Moon")
	{
		X = X_KA;
		Y = Y_KA;
		Z = Z_KA;
	}
	if (main_center_flag == "Earth")
	{
		std::vector <double> Moon_state(6); // вектор параметров Луны
		double X_M, Y_M, Z_M, Vx_M, Vy_M, Vz_M;

		Moon_state = get_planet_parameters(Date_time.get_string().c_str(), "MOON", "J2000", "EARTH");
		X_M = Moon_state[0]; // вектор Земля - Луна
		Y_M = Moon_state[1];
		Z_M = Moon_state[2];

		X = X_KA - X_M; //пересчитываю координаты
		Y = Y_KA - Y_M;
		Z = Z_KA - Z_M;
	}

	Omega_L = 13.17635815; // град/сутки
	str2et_c(Date.get_string().c_str(), &delta_t); // s, seconds past J2000 epoch
	delta_t /= (60 * 60 * 24); // Разница между датой и J2000 в сутках
	T = delta_t / 36525;

	e1 = (125.045 - 0.0529921 * delta_t) * Pi / 180;  //Сразу все в радианы перевожу
	e2 = (250.089 - 0.1059842 * delta_t) * Pi / 180;
	e3 = (260.008 + 13.0120009 * delta_t) * Pi / 180;
	e4 = (176.625 + 13.3407154 * delta_t) * Pi / 180;
	e5 = (357.529 + 0.9856003 * delta_t) * Pi / 180;
	e6 = (311.589 + 26.4057084 * delta_t) * Pi / 180;
	e7 = (134.963 + 13.0649930 * delta_t) * Pi / 180;
	e8 = (276.617 + 0.3287146 * delta_t) * Pi / 180;
	e9 = (34.226 + 1.7484877 * delta_t) * Pi / 180;
	e10 = (15.134 - 0.1589763 * delta_t) * Pi / 180;
	e11 = (119.743 + 0.0036096 * delta_t) * Pi / 180;
	e12 = (239.961 + 0.1643573 * delta_t) * Pi / 180;
	e13 = (25.053 + 12.9590088 * delta_t) * Pi / 180;

	alpha = 269.9949 + 0.0031 * T - 3.8787 * sin(e1) - 0.1204 * sin(e2) + 0.07 * sin(e3) -
		0.0172 * sin(e4) + 0.0072 * sin(e6) - 0.0052 * sin(e10) + 0.0043 * sin(e13); // в градусах
	beta = 66.5392 + 0.013 * T + 1.5419 * cos(e1) + 0.0239 * cos(e2) - 0.0278 * cos(e3) +
		0.0068 * cos(e4) - 0.0029 * cos(e6) + 0.0009 * cos(e7) + 0.0008 * cos(e10) - 0.0009 * cos(e13); // в градусах
	omega_delta_t = 38.3213 + Omega_L * delta_t - 1.4 * pow(10, -12) * delta_t * delta_t + 3.561 * sin(e1) +
		0.1208 * sin(e2) - 0.0642 * sin(e3) + 0.0158 * sin(e4) + 0.0252 * sin(e5) - 0.0066 * sin(e6) - 0.0047 * sin(e7) -
		0.0046 * sin(e8) + 0.0028 * sin(e9) + 0.0052 * sin(e10) + 0.004 * sin(e11) + 0.0019 * sin(e12) - 0.0044 * sin(e13); //в градусах

	a = alpha * Pi / 180;
	b = beta * Pi / 180;
	om = omega_delta_t * Pi / 180;// перевел их в радианы

	A.matr[0][0] = -sin(a) * cos(om) - cos(a) * sin(b) * sin(om);
	A.matr[0][1] = cos(a) * cos(om) - sin(a) * sin(b) * sin(om);
	A.matr[0][2] = sin(om) * cos(b);
	A.matr[1][0] = sin(a) * sin(om) - cos(a) * sin(b) * cos(om);
	A.matr[1][1] = -cos(a) * sin(om) - sin(a) * sin(b) * cos(om);
	A.matr[1][2] = cos(om) * cos(b);
	A.matr[2][0] = cos(a) * cos(b);
	A.matr[2][1] = sin(a) * cos(b);
	A.matr[2][2] = sin(b);
	return A;
}

std::vector <double> Spacecraft::J2000_to_SVSK_3x3()
{
	std::vector <double> SVSK_Vector;
	std::vector <double> SVSK_coords(3);
	std::vector <double> SVSK_velocities(3);

	std::vector <double> J2000_coords(3);
	std::vector <double> J2000_velocities(3);
	J2000_coords[0] = X_KA;
	J2000_coords[1] = Y_KA;
	J2000_coords[2] = Z_KA;
	J2000_velocities[0] = Vx_KA;
	J2000_velocities[1] = Vy_KA;
	J2000_velocities[2] = Vz_KA;

	SVSK_coords = J2000_to_SVSK_matrix_3x3(Date_time) * J2000_coords;
	SVSK_velocities = J2000_to_SVSK_matrix_3x3(Date_time) * J2000_velocities;

	SVSK_Vector.reserve(6); // preallocate memory
	SVSK_Vector.insert(SVSK_Vector.end(), SVSK_coords.begin(), SVSK_coords.end());
	SVSK_Vector.insert(SVSK_Vector.end(), SVSK_velocities.begin(), SVSK_velocities.end());
	return SVSK_Vector;
}

std::vector <double> Spacecraft::SVSK_to_J2000_3x3(std::vector <double> SVSK, Date_time_ms Date_time)
{
	std::vector <double> SVSK_coords(3);
	std::vector <double> SVSK_velocities(3);
	std::vector <double> J2000_coords(3);
	std::vector <double> J2000_velocities(3);
	std::vector <double> J2000_Vector;

	SVSK_coords[0] = SVSK[0];
	SVSK_coords[1] = SVSK[1];
	SVSK_coords[2] = SVSK[2];
	SVSK_velocities[0] = SVSK[3];
	SVSK_velocities[1] = SVSK[4];
	SVSK_velocities[2] = SVSK[5];

	matrix SVSK_J2000_3x3 = J2000_to_SVSK_matrix_3x3(Date_time).Transpose();

	J2000_coords = SVSK_J2000_3x3 * SVSK_coords;
	J2000_velocities = SVSK_J2000_3x3 * SVSK_velocities;

	J2000_Vector.reserve(6); // preallocate memory
	J2000_Vector.insert(J2000_Vector.end(), J2000_coords.begin(), J2000_coords.end());
	J2000_Vector.insert(J2000_Vector.end(), J2000_velocities.begin(), J2000_velocities.end());
	return J2000_Vector;
}

std::vector <double> Spacecraft::J2000_to_SVSK_6x6()
{
	std::vector <double> SVSK_coords(6);
	std::vector <double> J2000_coords(6);
	if (main_center_flag == "Moon")
	{
		J2000_coords[0] = X_KA;
		J2000_coords[1] = Y_KA;
		J2000_coords[2] = Z_KA;
		J2000_coords[3] = Vx_KA;
		J2000_coords[4] = Vy_KA;
		J2000_coords[5] = Vz_KA;
	}
	if (main_center_flag == "Earth")
	{
		std::vector <double> Moon_state(6); // вектор параметров Луны
		double X_M, Y_M, Z_M, Vx_M, Vy_M, Vz_M;
		double X_M0, Y_M0, Z_M0, Vx_M0, Vy_M0, Vz_M0;

		Moon_state = get_planet_parameters(Date_time.get_string().c_str(), "MOON", "J2000", "EARTH");
		X_M = Moon_state[0]; // вектор Земля - Луна
		Y_M = Moon_state[1];
		Z_M = Moon_state[2];
		Vx_M = Moon_state[3]; // вектор Земля - Луна
		Vy_M = Moon_state[4];
		Vz_M = Moon_state[5];
		X_M0 = X_KA - X_M; // считаем вектор от Луны на КА
		Y_M0 = Y_KA - Y_M;
		Z_M0 = Z_KA - Z_M;
		Vx_M0 = Vx_KA - Vx_M; // считаем вектор от Луны на КА
		Vy_M0 = Vy_KA - Vy_M;
		Vz_M0 = Vz_KA - Vz_M;

		J2000_coords[0] = X_M0;
		J2000_coords[1] = Y_M0;
		J2000_coords[2] = Z_M0;
		J2000_coords[3] = Vx_M0;
		J2000_coords[4] = Vy_M0;
		J2000_coords[5] = Vz_M0;
	}
	matrix A_6x6(6, 6);
	matrix SVSK_J2000_3x3 = J2000_to_SVSK_matrix_3x3(Date_time).Transpose();

	Date_time_ms Date_time_plus = Date_time;
	Date_time_plus.ms = Date_time_plus.ms + 500;
	Date_time_plus.update_date_time(); // прибавил 0.5 сек

	Date_time_ms Date_time_minus = Date_time;
	Date_time_minus.ms = Date_time_minus.ms - 500;
	Date_time_minus.update_date_time(); // отнял 0.5 сек

	matrix Helping = J2000_to_SVSK_matrix_3x3(Date_time_plus).Transpose() - J2000_to_SVSK_matrix_3x3(Date_time_minus).Transpose();
	A_6x6.matr[0][0] = SVSK_J2000_3x3.matr[0][0];
	A_6x6.matr[0][1] = SVSK_J2000_3x3.matr[0][1];
	A_6x6.matr[0][2] = SVSK_J2000_3x3.matr[0][2];
	A_6x6.matr[0][3] = 0;
	A_6x6.matr[0][4] = 0;
	A_6x6.matr[0][5] = 0;
	A_6x6.matr[1][0] = SVSK_J2000_3x3.matr[1][0];
	A_6x6.matr[1][1] = SVSK_J2000_3x3.matr[1][1];
	A_6x6.matr[1][2] = SVSK_J2000_3x3.matr[1][2];
	A_6x6.matr[1][3] = 0;
	A_6x6.matr[1][4] = 0;
	A_6x6.matr[1][5] = 0;
	A_6x6.matr[2][0] = SVSK_J2000_3x3.matr[2][0];
	A_6x6.matr[2][1] = SVSK_J2000_3x3.matr[2][1];
	A_6x6.matr[2][2] = SVSK_J2000_3x3.matr[2][2];
	A_6x6.matr[2][3] = 0;
	A_6x6.matr[2][4] = 0;
	A_6x6.matr[2][5] = 0;
	A_6x6.matr[3][0] = Helping.matr[0][0];
	A_6x6.matr[3][1] = Helping.matr[0][1];
	A_6x6.matr[3][2] = Helping.matr[0][2];
	A_6x6.matr[3][3] = SVSK_J2000_3x3.matr[0][0];
	A_6x6.matr[3][4] = SVSK_J2000_3x3.matr[0][1];
	A_6x6.matr[3][5] = SVSK_J2000_3x3.matr[0][2];
	A_6x6.matr[4][0] = Helping.matr[1][0];
	A_6x6.matr[4][1] = Helping.matr[1][1];
	A_6x6.matr[4][2] = Helping.matr[1][2];
	A_6x6.matr[4][3] = SVSK_J2000_3x3.matr[1][0];
	A_6x6.matr[4][4] = SVSK_J2000_3x3.matr[1][1];
	A_6x6.matr[4][5] = SVSK_J2000_3x3.matr[1][2];
	A_6x6.matr[5][0] = Helping.matr[2][0];
	A_6x6.matr[5][1] = Helping.matr[2][1];
	A_6x6.matr[5][2] = Helping.matr[2][2];
	A_6x6.matr[5][3] = SVSK_J2000_3x3.matr[2][0];
	A_6x6.matr[5][4] = SVSK_J2000_3x3.matr[2][1];
	A_6x6.matr[5][5] = SVSK_J2000_3x3.matr[2][2];

	matrix J2000_SVSK_6x6 = A_6x6.Transpose();
	SVSK_coords = J2000_SVSK_6x6 * J2000_coords;

	return SVSK_coords;
}

std::vector <double> Spacecraft::SVSK_to_J2000_6x6()
{
	std::vector <double> SVSK_coords(6);
	std::vector <double> J2000_coords(6);
	if (main_center_flag == "Moon")
	{
		SVSK_coords[0] = X_KA;
		SVSK_coords[1] = Y_KA;
		SVSK_coords[2] = Z_KA;
		SVSK_coords[3] = Vx_KA;
		SVSK_coords[4] = Vy_KA;
		SVSK_coords[5] = Vz_KA;
	}
	if (main_center_flag == "Earth")
	{
		std::vector <double> Moon_state(6); // вектор параметров Луны
		double X_M, Y_M, Z_M, Vx_M, Vy_M, Vz_M;
		double X_M0, Y_M0, Z_M0, Vx_M0, Vy_M0, Vz_M0;

		Moon_state = get_planet_parameters(Date_time.get_string().c_str(), "MOON", "J2000", "EARTH");
		X_M = Moon_state[0]; // вектор Земля - Луна
		Y_M = Moon_state[1];
		Z_M = Moon_state[2];
		Vx_M = Moon_state[3]; // вектор Земля - Луна
		Vy_M = Moon_state[4];
		Vz_M = Moon_state[5];
		X_M0 = X_KA - X_M; // считаем вектор от Луны на КА
		Y_M0 = Y_KA - Y_M;
		Z_M0 = Z_KA - Z_M;
		Vx_M0 = Vx_KA - Vx_M; // считаем вектор от Луны на КА
		Vy_M0 = Vy_KA - Vy_M;
		Vz_M0 = Vz_KA - Vz_M;

		SVSK_coords[0] = X_M0;
		SVSK_coords[1] = Y_M0;
		SVSK_coords[2] = Z_M0;
		SVSK_coords[3] = Vx_M0;
		SVSK_coords[4] = Vy_M0;
		SVSK_coords[5] = Vz_M0;
	}
	matrix A_6x6(6, 6);
	matrix SVSK_J2000_3x3 = J2000_to_SVSK_matrix_3x3(Date_time).Transpose();

	Date_time_ms Date_time_plus = Date_time;
	Date_time_plus.ms = Date_time_plus.ms + 500;
	Date_time_plus.update_date_time(); // прибавил 0.5 сек

	Date_time_ms Date_time_minus = Date_time;
	Date_time_minus.ms = Date_time_minus.ms - 500;
	Date_time_minus.update_date_time(); // отнял 0.5 сек

	matrix Helping = J2000_to_SVSK_matrix_3x3(Date_time_plus).Transpose() - J2000_to_SVSK_matrix_3x3(Date_time_minus).Transpose();
	A_6x6.matr[0][0] = SVSK_J2000_3x3.matr[0][0];
	A_6x6.matr[0][1] = SVSK_J2000_3x3.matr[0][1];
	A_6x6.matr[0][2] = SVSK_J2000_3x3.matr[0][2];
	A_6x6.matr[0][3] = 0;
	A_6x6.matr[0][4] = 0;
	A_6x6.matr[0][5] = 0;
	A_6x6.matr[1][0] = SVSK_J2000_3x3.matr[1][0];
	A_6x6.matr[1][1] = SVSK_J2000_3x3.matr[1][1];
	A_6x6.matr[1][2] = SVSK_J2000_3x3.matr[1][2];
	A_6x6.matr[1][3] = 0;
	A_6x6.matr[1][4] = 0;
	A_6x6.matr[1][5] = 0;
	A_6x6.matr[2][0] = SVSK_J2000_3x3.matr[2][0];
	A_6x6.matr[2][1] = SVSK_J2000_3x3.matr[2][1];
	A_6x6.matr[2][2] = SVSK_J2000_3x3.matr[2][2];
	A_6x6.matr[2][3] = 0;
	A_6x6.matr[2][4] = 0;
	A_6x6.matr[2][5] = 0;
	A_6x6.matr[3][0] = Helping.matr[0][0];
	A_6x6.matr[3][1] = Helping.matr[0][1];
	A_6x6.matr[3][2] = Helping.matr[0][2];
	A_6x6.matr[3][3] = SVSK_J2000_3x3.matr[0][0];
	A_6x6.matr[3][4] = SVSK_J2000_3x3.matr[0][1];
	A_6x6.matr[3][5] = SVSK_J2000_3x3.matr[0][2];
	A_6x6.matr[4][0] = Helping.matr[1][0];
	A_6x6.matr[4][1] = Helping.matr[1][1];
	A_6x6.matr[4][2] = Helping.matr[1][2];
	A_6x6.matr[4][3] = SVSK_J2000_3x3.matr[1][0];
	A_6x6.matr[4][4] = SVSK_J2000_3x3.matr[1][1];
	A_6x6.matr[4][5] = SVSK_J2000_3x3.matr[1][2];
	A_6x6.matr[5][0] = Helping.matr[2][0];
	A_6x6.matr[5][1] = Helping.matr[2][1];
	A_6x6.matr[5][2] = Helping.matr[2][2];
	A_6x6.matr[5][3] = SVSK_J2000_3x3.matr[2][0];
	A_6x6.matr[5][4] = SVSK_J2000_3x3.matr[2][1];
	A_6x6.matr[5][5] = SVSK_J2000_3x3.matr[2][2];

	matrix SVSK_J2000_6x6 = A_6x6;
	J2000_coords = SVSK_J2000_6x6 * SVSK_coords;

	return J2000_coords;
}

std::vector <double> Spacecraft::SVSK_to_J2000_6x6(std::vector <double> SVSK, Date_time_ms Date_time)
{
	std::vector <double> SVSK_coords(6);
	std::vector <double> J2000_coords(6);

	SVSK_coords[0] = SVSK[0];
	SVSK_coords[1] = SVSK[1];
	SVSK_coords[2] = SVSK[2];
	SVSK_coords[3] = SVSK[3];
	SVSK_coords[4] = SVSK[4];
	SVSK_coords[5] = SVSK[5];

	matrix A_6x6(6, 6);
	matrix SVSK_J2000_3x3 = J2000_to_SVSK_matrix_3x3(Date_time).Transpose();

	Date_time_ms Date_time_plus = Date_time;
	Date_time_plus.ms = Date_time_plus.ms + 500;
	Date_time_plus.update_date_time(); // прибавил 0.5 сек

	Date_time_ms Date_time_minus = Date_time;
	Date_time_minus.ms = Date_time_minus.ms - 500;
	Date_time_minus.update_date_time(); // отнял 0.5 сек

	matrix Helping = J2000_to_SVSK_matrix_3x3(Date_time_plus).Transpose() - J2000_to_SVSK_matrix_3x3(Date_time_minus).Transpose();
	A_6x6.matr[0][0] = SVSK_J2000_3x3.matr[0][0];
	A_6x6.matr[0][1] = SVSK_J2000_3x3.matr[0][1];
	A_6x6.matr[0][2] = SVSK_J2000_3x3.matr[0][2];
	A_6x6.matr[0][3] = 0;
	A_6x6.matr[0][4] = 0;
	A_6x6.matr[0][5] = 0;
	A_6x6.matr[1][0] = SVSK_J2000_3x3.matr[1][0];
	A_6x6.matr[1][1] = SVSK_J2000_3x3.matr[1][1];
	A_6x6.matr[1][2] = SVSK_J2000_3x3.matr[1][2];
	A_6x6.matr[1][3] = 0;
	A_6x6.matr[1][4] = 0;
	A_6x6.matr[1][5] = 0;
	A_6x6.matr[2][0] = SVSK_J2000_3x3.matr[2][0];
	A_6x6.matr[2][1] = SVSK_J2000_3x3.matr[2][1];
	A_6x6.matr[2][2] = SVSK_J2000_3x3.matr[2][2];
	A_6x6.matr[2][3] = 0;
	A_6x6.matr[2][4] = 0;
	A_6x6.matr[2][5] = 0;
	A_6x6.matr[3][0] = Helping.matr[0][0];
	A_6x6.matr[3][1] = Helping.matr[0][1];
	A_6x6.matr[3][2] = Helping.matr[0][2];
	A_6x6.matr[3][3] = SVSK_J2000_3x3.matr[0][0];
	A_6x6.matr[3][4] = SVSK_J2000_3x3.matr[0][1];
	A_6x6.matr[3][5] = SVSK_J2000_3x3.matr[0][2];
	A_6x6.matr[4][0] = Helping.matr[1][0];
	A_6x6.matr[4][1] = Helping.matr[1][1];
	A_6x6.matr[4][2] = Helping.matr[1][2];
	A_6x6.matr[4][3] = SVSK_J2000_3x3.matr[1][0];
	A_6x6.matr[4][4] = SVSK_J2000_3x3.matr[1][1];
	A_6x6.matr[4][5] = SVSK_J2000_3x3.matr[1][2];
	A_6x6.matr[5][0] = Helping.matr[2][0];
	A_6x6.matr[5][1] = Helping.matr[2][1];
	A_6x6.matr[5][2] = Helping.matr[2][2];
	A_6x6.matr[5][3] = SVSK_J2000_3x3.matr[2][0];
	A_6x6.matr[5][4] = SVSK_J2000_3x3.matr[2][1];
	A_6x6.matr[5][5] = SVSK_J2000_3x3.matr[2][2];

	matrix SVSK_J2000_6x6 = A_6x6;
	J2000_coords = SVSK_J2000_6x6 * SVSK_coords;

	return J2000_coords;
}

double Spacecraft::SVSK_inclination()
{
	std::vector <double> SVSK_coords(6);
	SVSK_coords = J2000_to_SVSK_6x6();

	double X = SVSK_coords[0];
	double Y = SVSK_coords[1];
	double Z = SVSK_coords[2];
	double Vx = SVSK_coords[3];
	double Vy = SVSK_coords[4];
	double Vz = SVSK_coords[5];

	double C_x = Y * Vz - Z * Vy;
	double C_y = -X * Vz + Z * Vx;
	double C_z = X * Vy - Y * Vx;
	double C = sqrt(C_x * C_x + C_y * C_y + C_z * C_z);
	inc = acos(C_z / C);    //inc КА
	if (inc == 0) { inc = 0.0000000001; }
	return inc;
}

std::vector <double> Spacecraft::get_planet_parameters(const char Date_time[], const char target_body[], const char reference_frame[],
	const char observer_body[])
{
	SpiceDouble   et_s;
	SpiceDouble   state[6];
	SpiceDouble   lt_s;
	str2et_c(Date_time, &et_s); // https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/str2et_c.html
	spkezr_c(target_body,// target body name or NAIF ID, https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
		et_s,		// s, seconds past J2000 epoch, https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/time.html
		reference_frame,	// reference frame, https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html
		"NONE",
		observer_body,	// observer body name or NAIF ID, https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
		state,		// output parameter, body state, which is [X km, Y km, Z km, Vx km/s, Vy km/s, Vz km/s]
		&lt_s);		// s, the one-way light time between the observer and target

	std::vector <double> Parameters(6);
	for (size_t i = 0; i < 6; i++)
	{
		Parameters[i] = state[i];
	}
	return Parameters;
	/*std::cout << "X = " << state[0] << " km" << std::endl
		<< "Y = " << state[1] << " km" << std::endl
		<< "Z = " << state[2] << " km" << std::endl
		<< "Vx = " << state[3] << " km/s" << std::endl
		<< "Vy = " << state[4] << " km/s" << std::endl
		<< "Vz = " << state[5] << " km/s" << std::endl;*/
}

void Spacecraft::Sun_perturbations(double& aS_x, double& aS_y, double& aS_z)
{
	std::vector <double> Sun_state(6); // вектор параметров Солнца
	double X_S, Y_S, Z_S, r_S, qs, fs;
	double X_S0, Y_S0, Z_S0, r_S0;


	if (main_center_flag == "Earth")
	{
		Sun_state = get_planet_parameters(Date_time.get_string().c_str(), "SUN", "J2000", "EARTH");
		X_S = Sun_state[0]; // вектор Земля - Солнце
		Y_S = Sun_state[1];
		Z_S = Sun_state[2];
	}
	if (main_center_flag == "Moon")
	{
		Sun_state = get_planet_parameters(Date_time.get_string().c_str(), "SUN", "J2000", "MOON");
		X_S = Sun_state[0]; // вектор Луна - Солнце
		Y_S = Sun_state[1];
		Z_S = Sun_state[2];
	}
	if (main_center_flag == "Apophis")
	{
		Sun_state = get_planet_parameters(Date_time.get_string().c_str(), "SUN", "J2000", "-40999");
		X_S = Sun_state[0]; // вектор Апофис - Солнце
		Y_S = Sun_state[1];
		Z_S = Sun_state[2];
	}
	r_S = sqrt(X_S * X_S + Y_S * Y_S + Z_S * Z_S);
	X_S0 = X_KA - X_S; // считаем вектор от Солнца на КА
	Y_S0 = Y_KA - Y_S;
	Z_S0 = Z_KA - Z_S;
	r_S0 = sqrt(X_S0 * X_S0 + Y_S0 * Y_S0 + Z_S0 * Z_S0);
	qs = X_KA / r_S * (X_KA / r_S - 2 * X_S / r_S) + Y_KA / r_S * (Y_KA / r_S - 2 * Y_S / r_S) + Z_KA / r_S * (Z_KA / r_S - 2 * Z_S / r_S);
	fs = (3 + 3 * qs + qs * qs) / (pow(1 + qs, 3 / 2) + 1) * qs;
	aS_x = -mu_Sun / (r_S0 * r_S0 * r_S0) * (X_KA + fs * X_S); // считаем возмущающее ускорение от Солнца
	aS_y = -mu_Sun / (r_S0 * r_S0 * r_S0) * (Y_KA + fs * Y_S);
	aS_z = -mu_Sun / (r_S0 * r_S0 * r_S0) * (Z_KA + fs * Z_S);
}

void Spacecraft::Moon_perturbations(double& aM_x, double& aM_y, double& aM_z)
{
	std::vector <double> Moon_state(6); // вектор параметров Луны
	double X_M, Y_M, Z_M, r_M, qm, fm;
	double X_M0, Y_M0, Z_M0, r_M0;

	if (main_center_flag == "Earth")
	{
		Moon_state = get_planet_parameters(Date_time.get_string().c_str(), "MOON", "J2000", "EARTH");
		X_M = Moon_state[0]; // вектор Земля - Луна
		Y_M = Moon_state[1];
		Z_M = Moon_state[2];
	}
	r_M = sqrt(X_M * X_M + Y_M * Y_M + Z_M * Z_M);
	X_M0 = X_KA - X_M; // считаем вектор от Луны на КА
	Y_M0 = Y_KA - Y_M;
	Z_M0 = Z_KA - Z_M;
	r_M0 = sqrt(X_M0 * X_M0 + Y_M0 * Y_M0 + Z_M0 * Z_M0);
	qm = X_KA / r_M * (X_KA / r_M - 2 * X_M / r_M) + Y_KA / r_M * (Y_KA / r_M - 2 * Y_M / r_M) + Z_KA / r_M * (Z_KA / r_M - 2 * Z_M / r_M);
	fm = (3 + 3 * qm + qm * qm) / (pow(1 + qm, 3 / 2) + 1) * qm;
	aM_x = -mu_Moon / (r_M0 * r_M0 * r_M0) * (X_KA + fm * X_M); // считаем возмущающее ускорение от Луны
	aM_y = -mu_Moon / (r_M0 * r_M0 * r_M0) * (Y_KA + fm * Y_M);
	aM_z = -mu_Moon / (r_M0 * r_M0 * r_M0) * (Z_KA + fm * Z_M);
}

void Spacecraft::Earth_perturbations(double& aE_x, double& aE_y, double& aE_z)
{
	std::vector <double> Earth_state(6); // вектор параметров Земли относительно Луны
	double X_E, Y_E, Z_E, r_E, qe, fe;
	double X_E0, Y_E0, Z_E0, r_E0;

	if (main_center_flag == "Moon")
	{
		Earth_state = get_planet_parameters(Date_time.get_string().c_str(), "EARTH", "J2000", "MOON");
		X_E = Earth_state[0]; // вектор Луна - Земля
		Y_E = Earth_state[1];
		Z_E = Earth_state[2];
	}
	if (main_center_flag == "Apophis")
	{
		Earth_state = get_planet_parameters(Date_time.get_string().c_str(), "EARTH", "J2000", "-40999");
		X_E = Earth_state[0]; // вектор Апофис - Земля
		Y_E = Earth_state[1];
		Z_E = Earth_state[2];
	}
	r_E = sqrt(X_E * X_E + Y_E * Y_E + Z_E * Z_E);
	X_E0 = X_KA - X_E; // считаем вектор от Земли на КА
	Y_E0 = Y_KA - Y_E;
	Z_E0 = Z_KA - Z_E;
	r_E0 = sqrt(X_E0 * X_E0 + Y_E0 * Y_E0 + Z_E0 * Z_E0);
	qe = X_KA / r_E * (X_KA / r_E - 2 * X_E / r_E) + Y_KA / r_E * (Y_KA / r_E - 2 * Y_E / r_E) + Z_KA / r_E * (Z_KA / r_E - 2 * Z_E / r_E);
	fe = (3 + 3 * qe + qe * qe) / (pow(1 + qe, 3 / 2) + 1) * qe;
	aE_x = -mu_Earth / (r_E0 * r_E0 * r_E0) * (X_KA + fe * X_E); // считаем возмущающее ускорение от Земли
	aE_y = -mu_Earth / (r_E0 * r_E0 * r_E0) * (Y_KA + fe * Y_E);
	aE_z = -mu_Earth / (r_E0 * r_E0 * r_E0) * (Z_KA + fe * Z_E);

}

void Spacecraft::Jupiter_perturbations(double& aJup_x, double& aJup_y, double& aJup_z)
{
	std::vector <double> Jupiter_state(6); // вектор параметров Юпитера относительно Солнца
	double X_Jup, Y_Jup, Z_Jup, r_Jup, qJup, fJup;
	double X_Jup0, Y_Jup0, Z_Jup0, r_Jup0;

	if (main_center_flag == "Apophis")
	{
		Jupiter_state = get_planet_parameters(Date_time.get_string().c_str(), "JUPITER Barycenter", "J2000", "-40999");
		X_Jup = Jupiter_state[0]; // вектор Апофис - Юпитер
		Y_Jup = Jupiter_state[1];
		Z_Jup = Jupiter_state[2];
	}
	r_Jup = sqrt(X_Jup * X_Jup + Y_Jup * Y_Jup + Z_Jup * Z_Jup);
	X_Jup0 = X_KA - X_Jup; // считаем вектор от Юпитера на КА
	Y_Jup0 = Y_KA - Y_Jup;
	Z_Jup0 = Z_KA - Z_Jup;
	r_Jup0 = sqrt(X_Jup0 * X_Jup0 + Y_Jup0 * Y_Jup0 + Z_Jup0 * Z_Jup0);
	qJup = X_KA / r_Jup * (X_KA / r_Jup - 2 * X_Jup / r_Jup) + Y_KA / r_Jup * (Y_KA / r_Jup - 2 * Y_Jup / r_Jup) + Z_KA / r_Jup * (Z_KA / r_Jup - 2 * Z_Jup / r_Jup);
	fJup = (3 + 3 * qJup + qJup * qJup) / (pow(1 + qJup, 3 / 2) + 1) * qJup;
	aJup_x = -mu_Jupiter / (r_Jup0 * r_Jup0 * r_Jup0) * (X_KA + fJup * X_Jup); // считаем возмущающее ускорение от Юпитера
	aJup_y = -mu_Jupiter / (r_Jup0 * r_Jup0 * r_Jup0) * (Y_KA + fJup * Y_Jup);
	aJup_z = -mu_Jupiter / (r_Jup0 * r_Jup0 * r_Jup0) * (Z_KA + fJup * Z_Jup);
}

void Spacecraft::Venus_perturbations(double& aVen_x, double& aVen_y, double& aVen_z)
{
	std::vector <double> Venus_state(6); // вектор параметров Венеры относительно Солнца
	double X_Ven, Y_Ven, Z_Ven, r_Ven, qVen, fVen;
	double X_Ven0, Y_Ven0, Z_Ven0, r_Ven0;

	if (main_center_flag == "Apophis")
	{
		Venus_state = get_planet_parameters(Date_time.get_string().c_str(), "VENUS", "J2000", "-40999");
		X_Ven = Venus_state[0]; // вектор Апофис - Венера
		Y_Ven = Venus_state[1];
		Z_Ven = Venus_state[2];
	}
	r_Ven = sqrt(X_Ven * X_Ven + Y_Ven * Y_Ven + Z_Ven * Z_Ven);
	X_Ven0 = X_KA - X_Ven; // считаем вектор от Венеры на КА
	Y_Ven0 = Y_KA - Y_Ven;
	Z_Ven0 = Z_KA - Z_Ven;
	r_Ven0 = sqrt(X_Ven0 * X_Ven0 + Y_Ven0 * Y_Ven0 + Z_Ven0 * Z_Ven0);
	qVen = X_KA / r_Ven * (X_KA / r_Ven - 2 * X_Ven / r_Ven) + Y_KA / r_Ven * (Y_KA / r_Ven - 2 * Y_Ven / r_Ven) + Z_KA / r_Ven * (Z_KA / r_Ven - 2 * Z_Ven / r_Ven);
	fVen = (3 + 3 * qVen + qVen * qVen) / (pow(1 + qVen, 3 / 2) + 1) * qVen;
	aVen_x = -mu_Venus / (r_Ven0 * r_Ven0 * r_Ven0) * (X_KA + fVen * X_Ven); // считаем возмущающее ускорение от Венеры
	aVen_y = -mu_Venus / (r_Ven0 * r_Ven0 * r_Ven0) * (Y_KA + fVen * Y_Ven);
	aVen_z = -mu_Venus / (r_Ven0 * r_Ven0 * r_Ven0) * (Z_KA + fVen * Z_Ven);
}

void Spacecraft::Mars_perturbations(double& aMar_x, double& aMar_y, double& aMar_z)
{
	std::vector <double> Mars_state(6); // вектор параметров Марса относительно Солнца
	double X_Mar, Y_Mar, Z_Mar, r_Mar, qMar, fMar;
	double X_Mar0, Y_Mar0, Z_Mar0, r_Mar0;

	if (main_center_flag == "Apophis")
	{
		Mars_state = get_planet_parameters(Date_time.get_string().c_str(), "MARS Barycenter", "J2000", "-40999");
		X_Mar = Mars_state[0]; // вектор Апофис - Марс
		Y_Mar = Mars_state[1];
		Z_Mar = Mars_state[2];
	}
	r_Mar = sqrt(X_Mar * X_Mar + Y_Mar * Y_Mar + Z_Mar * Z_Mar);
	X_Mar0 = X_KA - X_Mar; // считаем вектор от Марса на КА
	Y_Mar0 = Y_KA - Y_Mar;
	Z_Mar0 = Z_KA - Z_Mar;
	r_Mar0 = sqrt(X_Mar0 * X_Mar0 + Y_Mar0 * Y_Mar0 + Z_Mar0 * Z_Mar0);
	qMar = X_KA / r_Mar * (X_KA / r_Mar - 2 * X_Mar / r_Mar) + Y_KA / r_Mar * (Y_KA / r_Mar - 2 * Y_Mar / r_Mar) + Z_KA / r_Mar * (Z_KA / r_Mar - 2 * Z_Mar / r_Mar);
	fMar = (3 + 3 * qMar + qMar * qMar) / (pow(1 + qMar, 3 / 2) + 1) * qMar;
	aMar_x = -mu_Mars / (r_Mar0 * r_Mar0 * r_Mar0) * (X_KA + fMar * X_Mar); // считаем возмущающее ускорение от Марса
	aMar_y = -mu_Mars / (r_Mar0 * r_Mar0 * r_Mar0) * (Y_KA + fMar * Y_Mar);
	aMar_z = -mu_Mars / (r_Mar0 * r_Mar0 * r_Mar0) * (Z_KA + fMar * Z_Mar);
}

void Spacecraft::Sun_pressure(double& aSP_x, double& aSP_y, double& aSP_z)
{
	std::vector <double> Sun_state(6); // вектор параметров Солнца
	double X_S, Y_S, Z_S, r_S;
	double X_S0, Y_S0, Z_S0, r_S0;
	const double q0 = 1361.0 / 299792458.0; // давление на расстоянии 1 а.е.
	double q1;
	double tau_x, tau_y, tau_z;

	if (main_center_flag == "Earth")
	{
		Sun_state = get_planet_parameters(Date_time.get_string().c_str(), "SUN", "J2000", "EARTH");
		X_S = Sun_state[0]; // вектор Земля - Солнце
		Y_S = Sun_state[1];
		Z_S = Sun_state[2];
	}
	if (main_center_flag == "Moon")
	{
		Sun_state = get_planet_parameters(Date_time.get_string().c_str(), "SUN", "J2000", "MOON");
		X_S = Sun_state[0]; // вектор Луна - Солнце
		Y_S = Sun_state[1];
		Z_S = Sun_state[2];
	}
	if (main_center_flag == "Apophis")
	{
		Sun_state = get_planet_parameters(Date_time.get_string().c_str(), "SUN", "J2000", "-40999");
		X_S = Sun_state[0]; // вектор Апофис - Солнце
		Y_S = Sun_state[1];
		Z_S = Sun_state[2];
	}
	r_S = sqrt(X_S * X_S + Y_S * Y_S + Z_S * Z_S);
	X_S0 = X_KA - X_S; // считаем вектор от Солнца на КА
	Y_S0 = Y_KA - Y_S;
	Z_S0 = Z_KA - Z_S;
	r_S0 = sqrt(X_S0 * X_S0 + Y_S0 * Y_S0 + Z_S0 * Z_S0); //расстояние от КА до Солнца
	q1 = q0 * pow(Ae / r_S0, 2);
	tau_x = X_S0 / r_S0;
	tau_y = Y_S0 / r_S0;
	tau_z = Z_S0 / r_S0;

	aSP_x = koeff_otraj * q1 * Area * tau_x / mass / pow(10, 3); //перевел в км/c^2 (порядок 10^-11)
	aSP_y = koeff_otraj * q1 * Area * tau_y / mass / pow(10, 3);
	aSP_z = koeff_otraj * q1 * Area * tau_z / mass / pow(10, 3);
}

double Spacecraft::distance_to_Moon()
{
	if (main_center_flag == "Moon")
	{
		return sqrt(X_KA * X_KA + Y_KA * Y_KA + Z_KA * Z_KA);
	}
	else if (main_center_flag == "Earth")
	{
		std::vector <double> Moon_state(6); // вектор параметров Луны
		double X_M, Y_M, Z_M;
		double X_M0, Y_M0, Z_M0;

		Moon_state = get_planet_parameters(Date_time.get_string().c_str(), "MOON", "J2000", "EARTH");
		X_M = Moon_state[0]; // вектор Земля - Луна
		Y_M = Moon_state[1];
		Z_M = Moon_state[2];
		X_M0 = X_KA - X_M; // считаем вектор от Луны на КА
		Y_M0 = Y_KA - Y_M;
		Z_M0 = Z_KA - Z_M;
		return sqrt(X_M0 * X_M0 + Y_M0 * Y_M0 + Z_M0 * Z_M0);
	}
	else { return -1; }
}

double Spacecraft::distance_to_Earth()
{
	if (main_center_flag == "Earth")
	{
		return sqrt(X_KA * X_KA + Y_KA * Y_KA + Z_KA * Z_KA);
	}
	else if (main_center_flag == "Moon")
	{
		std::vector <double> Earth_state(6); // вектор параметров Земли относительно Луны
		double X_E, Y_E, Z_E;
		double X_E0, Y_E0, Z_E0;

		Earth_state = get_planet_parameters(Date_time.get_string().c_str(), "EARTH", "J2000", "MOON");
		X_E = Earth_state[0]; // вектор Луна - Земля
		Y_E = Earth_state[1];
		Z_E = Earth_state[2];
		X_E0 = X_KA - X_E; // считаем вектор от Земли на КА
		Y_E0 = Y_KA - Y_E;
		Z_E0 = Z_KA - Z_E;
		return sqrt(X_E0 * X_E0 + Y_E0 * Y_E0 + Z_E0 * Z_E0);
	}
	else { return -1; }
}

void Spacecraft::define_main_center()
{
	if (distance_to_Moon() <= 66000) //66000
	{
		if (main_center_flag == "Moon")
		{
			main_center_flag = "Moon";
		}
		if (main_center_flag == "Earth")
		{
			std::vector <double> Moon_state(6); // вектор параметров Луны
			double X_M, Y_M, Z_M, Vx_M, Vy_M, Vz_M;

			Moon_state = get_planet_parameters(Date_time.get_string().c_str(), "MOON", "J2000", "EARTH");
			X_M = Moon_state[0]; // вектор Земля - Луна
			Y_M = Moon_state[1];
			Z_M = Moon_state[2];
			Vx_M = Moon_state[3];
			Vy_M = Moon_state[4];
			Vz_M = Moon_state[5];

			X_KA = X_KA - X_M; //пересчитываю координаты и скорости
			Y_KA = Y_KA - Y_M;
			Z_KA = Z_KA - Z_M;
			Vx_KA = Vx_KA - Vx_M;
			Vy_KA = Vy_KA - Vy_M;
			Vz_KA = Vz_KA - Vz_M;
			// пересчитываю остальные параметры
			r_KA = sqrt(X_KA * X_KA + Y_KA * Y_KA + Z_KA * Z_KA);
			V_KA = sqrt(Vx_KA * Vx_KA + Vy_KA * Vy_KA + Vz_KA * Vz_KA);
			Grav_parameter = mu_Moon;
			main_center_flag = "Moon"; //переключил флаг
		}
	}
	else
	{
		if (main_center_flag == "Moon")
		{
			std::vector <double> Earth_state(6); // вектор параметров Земли относительно Луны
			double X_E, Y_E, Z_E, Vx_E, Vy_E, Vz_E;

			Earth_state = get_planet_parameters(Date_time.get_string().c_str(), "EARTH", "J2000", "MOON");
			X_E = Earth_state[0]; // вектор Луна - Земля
			Y_E = Earth_state[1];
			Z_E = Earth_state[2];
			Vx_E = Earth_state[3];
			Vy_E = Earth_state[4];
			Vz_E = Earth_state[5];

			X_KA = X_KA - X_E; //пересчитываю координаты и скорости
			Y_KA = Y_KA - Y_E;
			Z_KA = Z_KA - Z_E;
			Vx_KA = Vx_KA - Vx_E;
			Vy_KA = Vy_KA - Vy_E;
			Vz_KA = Vz_KA - Vz_E;
			// пересчитываю остальные параметры
			r_KA = sqrt(X_KA * X_KA + Y_KA * Y_KA + Z_KA * Z_KA);
			V_KA = sqrt(Vx_KA * Vx_KA + Vy_KA * Vy_KA + Vz_KA * Vz_KA);
			Grav_parameter = mu_Earth;
			main_center_flag = "Earth"; //переключил флаг
		}
		if (main_center_flag == "Earth")
		{
			main_center_flag = "Earth";
		}
	}
}

double Spacecraft::alpfa0(tm D1)
{
	double day_new, god_new, a_jd, y_jd, m_jd, JDN, JD, T, R_new, chas1_new, GST;
	char Year[8], Month[8], Day[8], Hour[8], Minute[8], Second[8];
	strftime(Year, 8, "%Y", &D1);
	strftime(Month, 8, "%m", &D1);
	strftime(Day, 8, "%d", &D1);
	strftime(Hour, 8, "%H", &D1);
	strftime(Minute, 8, "%M", &D1);
	strftime(Second, 8, "%S", &D1);

	day_new = atof(Day);
	god_new = atof(Year);
	a_jd = trunc((14 - atof(Month)) / 12);
	y_jd = god_new + 4800 - a_jd;
	m_jd = atof(Month) + 12 * a_jd - 3;
	JDN = day_new + trunc((153 * m_jd + 2) / 5) + 365 * y_jd + trunc(y_jd / 4) - trunc(y_jd / 100) + trunc(y_jd / 400) - 32045;
	JD = JDN + (atof(Hour) - 12) / 24 + atof(Minute) / 1440 + atof(Second) / 86400;
	T = (JD - 2451545) / 36525;
	R_new = 6.697374558 + 2400.051336 * T + 0.000025862 * T * T;
	while (R_new < 0)
	{
		R_new += 24;
	}
	while (R_new > 24)
	{
		R_new -= 24;
	}
	chas1_new = atof(Hour) + atof(Minute) / 60 + atof(Second) / 3600;
	GST = chas1_new * 1.002737909 + R_new;
	if (GST > 24) { GST -= 24; }
	if (GST < 0) { GST += 24; }
	return GST;
}

void Spacecraft::Finding_RAAN0_north(Date_time_ms Start_Date, double fi, double lambda, double& RAAN0_north)
{
	double sinA, x, alpha;
	sinA = abs(cos(inc0)) / cos(fi);
	x = asin(abs(sin(fi)) * sinA / sin(inc0));
	if (x < 0) { x = Pi - asin(abs(sin(fi)) * sinA / sin(inc0)); }
	alpha = alpfa0(Start_Date.Date);  //альфа в часах
	alpha = alpha * 15;           //в градусах
	alpha = alpha * toRad;      // в радианах

	RAAN0_north = alpha + lambda - x;
	while (RAAN0_north < 0) { RAAN0_north += 2 * Pi; }
	while (RAAN0_north > 2 * Pi) { RAAN0_north -= 2 * Pi; }
}

void Spacecraft::Finding_RAAN0_south(Date_time_ms Start_Date, double fi, double lambda, double& RAAN0_south)
{
	double sinA, x, alpha;
	sinA = abs(cos(inc0)) / cos(fi);
	x = asin(abs(sin(fi)) * sinA / sin(inc0));
	if (x < 0) { x = Pi - asin(abs(sin(fi)) * sinA / sin(inc0)); }
	alpha = alpfa0(Start_Date.Date);  //альфа в часах
	alpha = alpha * 15;           //в градусах
	alpha = alpha * toRad;      // в радианах

	RAAN0_south = alpha + lambda + x + Pi;
	while (RAAN0_south < 0) { RAAN0_south += 2 * Pi; }
	while (RAAN0_south > 2 * Pi) { RAAN0_south -= 2 * Pi; }
}

void Spacecraft::Moon_flight_without_perturbations_north(double flight_time /*время перелёта в сутках*/, Date_time_ms& approach_date)
{
	double fractional_flight_time, integer_flight_time;
	double X_targ, Y_targ, Z_targ, R_targ;
	double X_depart, Y_depart, Z_depart;
	double delta_a;
	double sin_u_targ, cos_u_targ, u_targ_plus, u_targ_minus, RAAN0_north, theta_targ, flight_time_1, E_anomaly_targ;
	//double bob;
	std::vector <double> Moon_state(6); // вектор параметров Луны

		// approach_date выбирать в начале суток, потом будет перебор вперед по минуте
		// flight_time в сутках
		// launch time в секундах
	inc = inc0;
	inc1 = inc0;
	double delta_t_sec = 0, integer_delta_t_sec = 0, fractional_delta_t_sec = 0;
	double koeff = 1, number = 0;
	do
	{
		number++;
		//if (number == 6) { koeff *=0.75; } // если уже долго топчется, то делаю шаги в два раза меньше
		approach_date.Date.tm_sec += integer_delta_t_sec * koeff; // 
		approach_date.update_date_time();
		approach_date.ms += 1000 * fractional_delta_t_sec * koeff;
		approach_date.update_date_time();
		// ниже поиск RAAN орбиты для попадания в Луну для заданного наклонения
		Moon_state = get_planet_parameters(approach_date.get_string().c_str(), "MOON", "J2000", "EARTH");
		X_targ = Moon_state[0]; // вектор Земля - Луна
		Y_targ = Moon_state[1];
		Z_targ = Moon_state[2];
		R_targ = sqrt(X_targ * X_targ + Y_targ * Y_targ + Z_targ * Z_targ);
		sin_u_targ = Z_targ / (R_targ * sin(inc));
		RAAN0 = 0.0000001;
		do
		{
			RAAN0 = RAAN0 + 0.000001;
			cos_u_targ = (X_targ / R_targ + sin(RAAN0) * sin_u_targ * cos(inc)) / (cos(RAAN0));
			if (cos_u_targ > 1)
			{
				delta_a = 1001;
				continue;
			}
			delta_a = abs(R_targ * (sin(RAAN0) * cos_u_targ + cos(RAAN0) * sin_u_targ * cos(inc)) - Y_targ);
		} while ((delta_a > 2) and (RAAN0 < 2 * Pi));

		// ниже поиск RAAN орбиты при старте с Земли в рассчитанное время. ДВУ должны совпасть
		Start_date = approach_date;
		fractional_flight_time = modf(flight_time, &integer_flight_time); // разделил время перелёта на целую и дробную части
		Start_date.Date.tm_mday -= integer_flight_time; // вычел время перелёта в сутках (целую часть)
		Start_date.update_date_time();
		Start_date.Date.tm_sec -= (fractional_flight_time * (24 * 60 * 60) + launch_time); // вычел время вывода на орбиту в секундах и дробную часть времени перелёта
		Start_date.update_date_time(); // получил время старта без учета движения по опорной орбите
		// исправить расчёт даты старта РН с Земли (учитывать еще движение по опорной орбите)!!!!!
		Finding_RAAN0_north(Start_date, fi, lambda, RAAN0_north); //долгота и широта Восточного
		if (RAAN0 > RAAN0_north)
		{
			if ((RAAN0 - RAAN0_north) < Pi)
			{
				delta_t_sec = (RAAN0 - RAAN0_north) / omega_earth;
				fractional_delta_t_sec = modf(delta_t_sec, &integer_delta_t_sec); // разделил дельта т на целую и дробную части
			}
			else
			{
				delta_t_sec = -(2 * Pi + RAAN0_north - RAAN0) / omega_earth;
				fractional_delta_t_sec = modf(delta_t_sec, &integer_delta_t_sec); // разделил дельта т на целую и дробную части
			}
		}
		else
		{
			if ((RAAN0_north - RAAN0) < Pi)
			{
				delta_t_sec = -(RAAN0_north - RAAN0) / omega_earth;
				fractional_delta_t_sec = modf(delta_t_sec, &integer_delta_t_sec); // разделил дельта т на целую и дробную части
			}
			else
			{
				delta_t_sec = (2 * Pi + RAAN0 - RAAN0_north) / omega_earth;
				fractional_delta_t_sec = modf(delta_t_sec, &integer_delta_t_sec); // разделил дельта т на целую и дробную части
			}
		}
	} while (abs(delta_t_sec) >= 0.3);

	u_targ_plus = acos(cos_u_targ);
	u_targ_minus = 2 * Pi - acos(cos_u_targ);
	time_t tt = mktime(&approach_date.Date); // перевожу в секунды, чтобы вычесть
	time_t tstr = tt - flight_time * (24 * 60 * 60);
	localtime_s(&tau1_date.Date, &tstr); // тау в дате 
	tau1_date.update_date_time();

	theta_targ = 1.57;      // начинаю перебор сразу с 90 град
	flight_time_1 = 0;
	while (abs(flight_time - flight_time_1) > 0.001 / (24))
	{
		theta_targ = theta_targ + 0.0001 * toRad;
		a1 = (a0 * (a0 - R_targ * cos(theta_targ))) / (2 * a0 - R_targ * (1 + cos(theta_targ))); //Суханова для гиперболы посмотреть или методу Казаковцева
		if (a1 < 0) { continue; }
		ecc1 = (R_targ - a0) / (a0 - R_targ * cos(theta_targ));
		E_anomaly_targ = 2 * atan(tan(theta_targ / 2) / sqrt((1 + ecc1) / (1 - ecc1)));
		if (E_anomaly_targ < 0) { E_anomaly_targ = E_anomaly_targ + Pi; }
		flight_time_1 = sqrt(a1 * a1 * a1 / mu_Earth) * (E_anomaly_targ - ecc1 * sin(E_anomaly_targ));
		flight_time_1 = flight_time_1 / (60 * 60 * 24); // перевел время полёта в сутки
	}
	omega1 = u_targ_plus - theta_targ; // два варианта u, надо решить с этим!
	if (omega1 < 0) { omega1 = omega1 + 2 * Pi; }
	RAAN1 = RAAN0;

	//om1:= Pi + om1; //не знаю почему, но так

	X_depart = (cos(RAAN1) * cos(omega1) - sin(RAAN1) * sin(omega1) * cos(inc1)) * a0; // определяю координаты точки подачи импульса
	Y_depart = (sin(RAAN1) * cos(omega1) + cos(RAAN1) * sin(omega1) * cos(inc1)) * a0;
	Z_depart = (sin(omega1) * sin(inc1)) * a0;
	u_depart = atan2(Z_depart / (a0 * sin(inc0)), (X_depart * cos(RAAN0) + Y_depart * sin(RAAN0)) / a0); // аргументы широты точки отлета
	if (u_depart < 0) { u_depart += 2 * Pi; }

	double EA_depart = 2 * atan(sqrt((1 - ecc0) / (1 + ecc0)) * tan(0.5 * u_depart));
	if (EA_depart < 0) { EA_depart += 2 * Pi; }
	double MA_depart = EA_depart - ecc0 * sin(EA_depart);

	delta_t_depart = sqrt(abs(a0 * a0 * a0) / mu_Earth) * MA_depart;


	fly_away_impulse = sqrt(2 * mu_Earth * (1 / a0 - 1 / (2 * a1))) - sqrt(mu_Earth / a0);

}

void Spacecraft::Moon_flight_without_perturbations_south(double flight_time /*время перелёта в сутках*/, Date_time_ms& approach_date)
{
	double fractional_flight_time, integer_flight_time;
	double X_targ, Y_targ, Z_targ, R_targ;
	double X_depart, Y_depart, Z_depart;
	double delta_a;
	double sin_u_targ, cos_u_targ, u_targ_plus, u_targ_minus, RAAN0_south, theta_targ, flight_time_1, E_anomaly_targ;
	//double bob;
	std::vector <double> Moon_state(6); // вектор параметров Луны

		// approach_date выбирать в начале суток, потом будет перебор вперед по минуте
		// flight_time в сутках
		// launch time в секундах
	inc = inc0;
	inc1 = inc0;
	double delta_t_sec = 0, integer_delta_t_sec = 0, fractional_delta_t_sec = 0;
	do
	{
		approach_date.Date.tm_sec += integer_delta_t_sec; // 
		approach_date.update_date_time();
		approach_date.ms += 1000 * fractional_delta_t_sec;
		approach_date.update_date_time();
		// ниже поиск RAAN орбиты для попадания в Луну для заданного наклонения
		Moon_state = get_planet_parameters(approach_date.get_string().c_str(), "MOON", "J2000", "EARTH");
		X_targ = Moon_state[0]; // вектор Земля - Луна
		Y_targ = Moon_state[1];
		Z_targ = Moon_state[2];
		R_targ = sqrt(X_targ * X_targ + Y_targ * Y_targ + Z_targ * Z_targ);
		sin_u_targ = Z_targ / (R_targ * sin(inc));
		RAAN0 = 0.0000001;
		do
		{
			RAAN0 = RAAN0 + 0.00001;
			cos_u_targ = (X_targ / R_targ + sin(RAAN0) * sin_u_targ * cos(inc)) / (cos(RAAN0));
			if (cos_u_targ > 1)
			{
				delta_a = 1001;
				continue;
			}
			delta_a = abs(R_targ * (sin(RAAN0) * cos_u_targ + cos(RAAN0) * sin_u_targ * cos(inc)) - Y_targ);
		} while ((delta_a > 10) and (RAAN0 < 2 * Pi));

		// ниже поиск RAAN орбиты при старте с Земли в рассчитанное время. ДВУ должны совпасть
		Start_date = approach_date;
		fractional_flight_time = modf(flight_time, &integer_flight_time); // разделил время перелёта на целую и дробную части
		Start_date.Date.tm_mday -= integer_flight_time; // вычел время перелёта в сутках (целую часть)
		Start_date.update_date_time();
		Start_date.Date.tm_sec -= (fractional_flight_time * (24 * 60 * 60) + launch_time); // вычел время вывода на орбиту в секундах и дробную часть времени перелёта
		Start_date.update_date_time(); // получил время старта без учета движения по опорной орбите
		// исправить расчёт даты старта РН с Земли (учитывать еще движение по опорной орбите)!!!!!
		Finding_RAAN0_south(Start_date, fi, lambda, RAAN0_south); //долгота и широта Восточного
		if (RAAN0 > RAAN0_south)
		{
			if ((RAAN0 - RAAN0_south) < Pi)
			{
				delta_t_sec = (RAAN0 - RAAN0_south) / omega_earth;
				fractional_delta_t_sec = modf(delta_t_sec, &integer_delta_t_sec); // разделил дельта т на целую и дробную части
			}
			else
			{
				delta_t_sec = -(2 * Pi + RAAN0_south - RAAN0) / omega_earth;
				fractional_delta_t_sec = modf(delta_t_sec, &integer_delta_t_sec); // разделил дельта т на целую и дробную части
			}
		}
		else
		{
			if ((RAAN0_south - RAAN0) < Pi)
			{
				delta_t_sec = -(RAAN0_south - RAAN0) / omega_earth;
				fractional_delta_t_sec = modf(delta_t_sec, &integer_delta_t_sec); // разделил дельта т на целую и дробную части
			}
			else
			{
				delta_t_sec = (2 * Pi + RAAN0 - RAAN0_south) / omega_earth;
				fractional_delta_t_sec = modf(delta_t_sec, &integer_delta_t_sec); // разделил дельта т на целую и дробную части
			}
		}
	} while (abs(delta_t_sec) >= 0.25);

	u_targ_plus = acos(cos_u_targ);
	u_targ_minus = 2 * Pi - acos(cos_u_targ);
	time_t tt = mktime(&approach_date.Date); // перевожу в секунды, чтобы вычесть
	time_t tstr = tt - flight_time * (24 * 60 * 60);
	localtime_s(&tau1_date.Date, &tstr); // тау в дате 
	tau1_date.update_date_time();

	theta_targ = 1.57;      // начинаю перебор сразу с 90 град
	flight_time_1 = 0;
	while (abs(flight_time - flight_time_1) > 0.001 / (24))
	{
		theta_targ = theta_targ + 0.0001 * toRad;
		a1 = (a0 * (a0 - R_targ * cos(theta_targ))) / (2 * a0 - R_targ * (1 + cos(theta_targ))); //Суханова для гиперболы посмотреть или методу Казаковцева
		if (a1 < 0) { continue; }
		ecc1 = (R_targ - a0) / (a0 - R_targ * cos(theta_targ));
		E_anomaly_targ = 2 * atan(tan(theta_targ / 2) / sqrt((1 + ecc1) / (1 - ecc1)));
		if (E_anomaly_targ < 0) { E_anomaly_targ = E_anomaly_targ + Pi; }
		flight_time_1 = sqrt(a1 * a1 * a1 / mu_Earth) * (E_anomaly_targ - ecc1 * sin(E_anomaly_targ));
		flight_time_1 = flight_time_1 / (60 * 60 * 24); // перевел время полёта в сутки
	}
	omega1 = u_targ_plus - theta_targ; // два варианта u, надо решить с этим!
	if (omega1 < 0) { omega1 = omega1 + 2 * Pi; }
	RAAN1 = RAAN0;

	//om1:= Pi + om1; //не знаю почему, но так

	X_depart = (cos(RAAN1) * cos(omega1) - sin(RAAN1) * sin(omega1) * cos(inc1)) * a0; // определяю координаты точки подачи импульса
	Y_depart = (sin(RAAN1) * cos(omega1) + cos(RAAN1) * sin(omega1) * cos(inc1)) * a0;
	Z_depart = (sin(omega1) * sin(inc1)) * a0;
	u_depart = atan2(Z_depart / (a0 * sin(inc0)), (X_depart * cos(RAAN0) + Y_depart * sin(RAAN0)) / a0); // аргументы широты точки отлета
	if (u_depart < 0) { u_depart += 2 * Pi; }

	double EA_depart = 2 * atan(sqrt((1 - ecc0) / (1 + ecc0)) * tan(0.5 * u_depart));
	if (EA_depart < 0) { EA_depart += 2 * Pi; }
	double MA_depart = EA_depart - ecc0 * sin(EA_depart);

	delta_t_depart = sqrt(abs(a0 * a0 * a0) / mu_Earth) * MA_depart;


	fly_away_impulse = sqrt(2 * mu_Earth * (1 / a0 - 1 / (2 * a1))) - sqrt(mu_Earth / a0);

}

//void Spacecraft::Moon_flight_without_perturbations(double flight_time /*время перелёта в сутках*/, tm &approach_date)
/*{
	double fractional_flight_time, integral_flight_time;
	double X_targ, Y_targ, Z_targ, R_targ;
	double X_depart, Y_depart, Z_depart;
	double delta_a;
	double sin_u_targ, cos_u_targ, u_targ, RAAN0_1, theta_targ, flight_time_1, E_anomaly_targ;
	//double bob;
	char approach_date_char[32];
	std::vector <double> Moon_state(6); // вектор параметров Луны

		// approach_date выбирать в начале суток, потом будет перебор вперед по минуте
		// flight_time в сутках
		// launch time в секундах
	inc1 = inc0;
	do
	{
		approach_date.tm_sec += 1; // +1
		mktime(&approach_date);
		// ниже поиск RAAN орбиты для попадания в Луну для заданного наклонения
		strftime(approach_date_char, 32, "%Y-%b-%d %H:%M:%S", &approach_date); // из tm в массив из char
		Moon_state = get_planet_parameters(approach_date_char, "MOON", "J2000", "EARTH");
		X_targ = Moon_state[0]; // вектор Земля - Луна
		Y_targ = Moon_state[1];
		Z_targ = Moon_state[2];
		R_targ = sqrt(X_targ * X_targ + Y_targ * Y_targ + Z_targ * Z_targ);
		sin_u_targ = Z_targ / (R_targ*sin(inc));
		RAAN0 = 0.0000001;
		//RAAN0:=3.14;    // если нужно второе решение для ДВУ
		do
		{
			RAAN0 = RAAN0 + 0.00001;
			cos_u_targ = (X_targ / R_targ + sin(RAAN0)*sin_u_targ*cos(inc)) / (cos(RAAN0));
			if (cos_u_targ > 1)
			{
			delta_a = 1001;
				continue;
			}
			delta_a = abs(R_targ*(sin(RAAN0)*cos_u_targ + cos(RAAN0)*sin_u_targ*cos(inc)) - Y_targ);
		} while ((delta_a > 10) and (RAAN0 < 2 * Pi));


			/* if RAAN0 >= 2 * Pi then
				begin
				u_pr : = Pi - arcsin(Z_pr / (R_pr*sin(i1)));
		RAAN0: = 0;
		delta_a: = abs(R_pr*(cos(RAAN0)*cos(u_pr) - sin(RAAN0)*sin(u_pr)*cos(i1)) - X_pr);
		delta_b: = abs(R_pr*(sin(RAAN0)*cos(u_pr) + cos(RAAN0)*sin(u_pr)*cos(i1)) - Y_pr);
			while (((delta_a > 1000) or (delta_b > 1000)) and (RAAN0 <= 2 * Pi)) do
				begin
				RAAN0 : = RAAN0 + 0.001;
		delta_a: = abs(R_pr*(cos(RAAN0)*cos(u_pr) - sin(RAAN0)*sin(u_pr)*cos(i1)) - X_pr);
		delta_b: = abs(R_pr*(sin(RAAN0)*cos(u_pr) + cos(RAAN0)*sin(u_pr)*cos(i1)) - Y_pr);
			end;
			end;  */
			// ниже поиск RAAN орбиты при старте с Земли в рассчитанное время. ДВУ должны совпасть
	/*		Start_date = approach_date;
			fractional_flight_time = modf(flight_time, &integral_flight_time); // разделил время перелёта на целую и дробную части
			Start_date.tm_mday -= integral_flight_time; // вычел время перелёта в сутках (целую часть)
			mktime(&Start_date);
			Start_date.tm_sec -= (fractional_flight_time * (24 * 60 * 60) + launch_time); // вычел время вывода на орбиту в секундах и дробную часть времени перелёта
			mktime(&Start_date); // получил время старта без учета движения по опорной орбите
			// исправить расчёт даты старта РН с Земли (учитывать еще движение по опорной орбите)!!!!!
			Finding_RAAN0(Start_date, fi, lambda, RAAN0_1); //долгота и широта Восточного
		}
		while (abs(RAAN0 - RAAN0_1) >= 0.05*toRad);

		u_targ = acos(cos_u_targ);
		time_t tt = mktime(&approach_date); // перевожу в секунды, чтобы вычесть
		time_t tstr = tt - flight_time * (24 * 60 * 60);
		localtime_s(&tau1_date, &tstr); // тау в дате

		theta_targ = 1.57;      // начинаю перебор сразу с 90 град
		flight_time_1 = 0;
		while (abs(flight_time - flight_time_1) > 0.001 / (24))
		{
			theta_targ = theta_targ + 0.0001*toRad;
			a1 = (a0*(a0 - R_targ * cos(theta_targ))) / (2 * a0 - R_targ * (1 + cos(theta_targ))); //Суханова для гиперболы посмотреть или методу Казаковцева
			if (a1 < 0) { continue; }
			ecc1 = (R_targ - a0) / (a0 - R_targ * cos(theta_targ));
			E_anomaly_targ = 2 * atan(tan(theta_targ / 2) / sqrt((1 + ecc1) / (1 - ecc1)));
			if (E_anomaly_targ < 0) { E_anomaly_targ = E_anomaly_targ + Pi; }
			flight_time_1 = sqrt(a1*a1*a1 / mu_Earth)*(E_anomaly_targ - ecc1 * sin(E_anomaly_targ));
			flight_time_1 = flight_time_1 / (60 * 60 * 24); // перевел время полёта в сутки
		}
		omega1 = u_targ - theta_targ;
		if (omega1 < 0) { omega1 = omega1 + 2 * Pi; }
		RAAN1 = RAAN0;

		//om1:= Pi + om1; //не знаю почему, но так

		X_depart = (cos(RAAN1)*cos(omega1) - sin(RAAN1)*sin(omega1)*cos(inc1))*a0; // определяю координаты точки подачи импульса
		Y_depart = (sin(RAAN1)*cos(omega1) + cos(RAAN1)*sin(omega1)*cos(inc1))*a0;
		Z_depart = (sin(omega1)*sin(inc1))*a0;
		u_depart = atan2(Z_depart / (a0*sin(inc0)), (X_depart*cos(RAAN0) + Y_depart * sin(RAAN0)) / a0); // аргументы широты точки отлета
		if (u_depart < 0) { u_depart += 2 * Pi; }

		double EA_depart = 2 * atan(sqrt((1 - ecc0) / (1 + ecc0))*tan(0.5*u_depart));
		if (EA_depart < 0) { EA_depart += 2 * Pi; }
		double MA_depart = EA_depart - ecc0 * sin(EA_depart);

		delta_t_depart =  sqrt(abs(a0*a0*a0) / mu_Earth)*MA_depart;


		fly_away_impulse = sqrt(2 * mu_Earth*(1 / a0 - 1 / (2 * a1))) - sqrt(mu_Earth / a0);

	}*/

void Spacecraft::calc_B_plane_parameters()
{
	double X_M0, Y_M0, Z_M0, r_M0, Vx_M0, Vy_M0, Vz_M0, V_M0;
	std::vector <double> Moon_state(6); // вектор параметров Луны
	double X_M, Y_M, Z_M, r_M, Vx_M, Vy_M, Vz_M;
	double h;
	Moon_state = get_planet_parameters(Date_time.get_string().c_str(), "MOON", "J2000", "EARTH");
	X_M = Moon_state[0]; // вектор Земля - Луна
	Y_M = Moon_state[1];
	Z_M = Moon_state[2];
	Vx_M = Moon_state[3];
	Vy_M = Moon_state[4];
	Vz_M = Moon_state[5];
	if (main_center_flag == "Earth")
	{
		X_M0 = X_KA - X_M; // считаем вектор от Луны на КА
		Y_M0 = Y_KA - Y_M;
		Z_M0 = Z_KA - Z_M;
		r_M0 = sqrt(X_M0 * X_M0 + Y_M0 * Y_M0 + Z_M0 * Z_M0);
		Vx_M0 = Vx_KA - Vx_M;
		Vy_M0 = Vy_KA - Vy_M;
		Vz_M0 = Vz_KA - Vz_M;
		V_M0 = sqrt(Vx_M0 * Vx_M0 + Vy_M0 * Vy_M0 + Vz_M0 * Vz_M0);
	}
	if (main_center_flag == "Moon")
	{
		X_M0 = X_KA; // считаем вектор от Луны на КА
		Y_M0 = Y_KA;
		Z_M0 = Z_KA;
		r_M0 = sqrt(X_M0 * X_M0 + Y_M0 * Y_M0 + Z_M0 * Z_M0);
		Vx_M0 = Vx_KA;
		Vy_M0 = Vy_KA;
		Vz_M0 = Vz_KA;
		V_M0 = sqrt(Vx_M0 * Vx_M0 + Vy_M0 * Vy_M0 + Vz_M0 * Vz_M0);
	}

	// единичный вектор нормали к орбите КА
	/*n_x = C_x / C;
	n_y = C_y / C;
	n_z = C_z / C;
	n = sqrt(n_x * n_x + n_y * n_y + n_z * n_z);*/
	// вектор эксцентриситета
	/*e_x = (V_M0*V_M0 / mu_Moon - 1 / r_M0)*X_M0 - (X_M0*Vx_M0 + Y_M0 * Vy_M0 + Z_M0 * Vz_M0) / mu_Moon * Vx_M0;
	e_y = (V_M0*V_M0 / mu_Moon - 1 / r_M0)*Y_M0 - (X_M0*Vx_M0 + Y_M0 * Vy_M0 + Z_M0 * Vz_M0) / mu_Moon * Vy_M0;
	e_z = (V_M0*V_M0 / mu_Moon - 1 / r_M0)*Z_M0 - (X_M0*Vx_M0 + Y_M0 * Vy_M0 + Z_M0 * Vz_M0) / mu_Moon * Vz_M0;
	e = sqrt(e_x * e_x + e_y * e_y + e_z * e_z);
	alpha_hyp = acos(1 / e);*/



	if (b_plane_first_time == 1) // если первый раз считаем Картинную плоскость то определяем скорость на бесконечности
	{
		/*V_inf_plus[0] = cos(alpha_hyp)*e_x + sin(alpha_hyp)*(n_y*e_z - n_z * e_y);
		V_inf_plus[1] = cos(alpha_hyp)*e_y + sin(alpha_hyp)*(n_z*e_x - n_x * e_z);
		V_inf_plus[2] = cos(alpha_hyp)*e_z + sin(alpha_hyp)*(n_x*e_y - n_y * e_x);
		double V_inf_plus_abs = sqrt(V_inf_plus[0] * V_inf_plus[0] + V_inf_plus[1] * V_inf_plus[1] + V_inf_plus[2] * V_inf_plus[2]);

		V_inf_plus[0] = V_inf_plus[0] / V_inf_plus_abs * V_inf; // делю на длину вектора и умножаю на величину скорости на бесконечности
		V_inf_plus[1] = V_inf_plus[1] / V_inf_plus_abs * V_inf;
		V_inf_plus[2] = V_inf_plus[2] / V_inf_plus_abs * V_inf;

		V_inf_min[0] = -cos(alpha_hyp)*e_x + sin(alpha_hyp)*(n_y*e_z - n_z * e_y);
		V_inf_min[1] = -cos(alpha_hyp)*e_y + sin(alpha_hyp)*(n_z*e_x - n_x * e_z);
		V_inf_min[2] = -cos(alpha_hyp)*e_z + sin(alpha_hyp)*(n_x*e_y - n_y * e_x);
		double V_inf_min_abs = sqrt(V_inf_min[0] * V_inf_min[0] + V_inf_min[1] * V_inf_min[1] + V_inf_min[2] * V_inf_min[2]);

		V_inf_min[0] = V_inf_min[0] / V_inf_min_abs * V_inf; // делю на длину вектора и умножаю на величину скорости на бесконечности
		V_inf_min[1] = V_inf_min[1] / V_inf_min_abs * V_inf;
		V_inf_min[2] = V_inf_min[2] / V_inf_min_abs * V_inf;*/

		//b_plane_first_time = 0;

		h = V_M0 * V_M0 - 2 * mu_Moon / r_M0; // инт энергии
		V_inf0 = sqrt(h);
		V_inf = sqrt(h);

		double C_x, C_y, C_z, C; // интеграл площадей орбиты относительно Луны!
		C_x = Y_M0 * Vz_M0 - Z_M0 * Vy_M0;
		C_y = Z_M0 * Vx_M0 - X_M0 * Vz_M0;
		C_z = X_M0 * Vy_M0 - Y_M0 * Vx_M0;
		C = sqrt(C_x * C_x + C_y * C_y + C_z * C_z);

		// вектор Лапласа
		double f_x = -Grav_parameter * X_M0 / r_M0 - Vz_M0 * C_y + Vy_M0 * C_z;
		double f_y = -Grav_parameter * Y_M0 / r_M0 - Vx_M0 * C_z + Vz_M0 * C_x;
		double f_z = -Grav_parameter * Z_M0 / r_M0 - Vy_M0 * C_x + Vx_M0 * C_y;
		double Lap = sqrt(f_x * f_x + f_y * f_y + f_z * f_z);

		double d_x, d_y, d_z, d; // вектор d = [c x f]
		d_x = C_y * f_z - C_z * f_y;
		d_y = C_z * f_x - C_x * f_z;
		d_z = C_x * f_y - C_y * f_x;
		d = sqrt(d_x * d_x + d_y * d_y + d_z * d_z);

		double B_x, B_y, B_z, B;// вектор прицельной дальности
		B_x = (C * C * f_x - mu_Moon * d_x / V_inf) / (Lap * Lap);
		B_y = (C * C * f_y - mu_Moon * d_y / V_inf) / (Lap * Lap);
		B_z = (C * C * f_z - mu_Moon * d_z / V_inf) / (Lap * Lap);
		B = sqrt(B_x * B_x + B_y * B_y + B_z * B_z);

		double V_inf_x, V_inf_y, V_inf_z; // вектор скорости на бесконечности
		V_inf_x = V_inf * (mu_Moon * f_x + V_inf * d_x) / (Lap * Lap);
		V_inf_y = V_inf * (mu_Moon * f_y + V_inf * d_y) / (Lap * Lap);
		V_inf_z = V_inf * (mu_Moon * f_z + V_inf * d_z) / (Lap * Lap);

		// определим направления векторов кси и этта, которые являются осями СК в Картинной плоскости
		// сначала определим едничный вектор нормали к орбите Луны
		double C_x_M, C_y_M, C_z_M, C_M; // интеграл площадей орбиты Луны!
		double N_x_M, N_y_M, N_z_M, N_M; // вектор нормали к орбите Луны
		C_x_M = Y_M * Vz_M - Z_M * Vy_M;
		C_y_M = Z_M * Vx_M - X_M * Vz_M;
		C_z_M = X_M * Vy_M - Y_M * Vx_M;
		C_M = sqrt(C_x_M * C_x_M + C_y_M * C_y_M + C_z_M * C_z_M);
		// единичный вектор нормали к орбите Луны
		N_x_M = C_x_M / C_M;
		N_y_M = C_y_M / C_M;
		N_z_M = C_z_M / C_M;
		N_M = sqrt(N_x_M * N_x_M + N_y_M * N_y_M + N_z_M * N_z_M);

		zeta_bplane[0] = -V_inf_x / V_inf;
		zeta_bplane[1] = -V_inf_y / V_inf;
		zeta_bplane[2] = -V_inf_z / V_inf;
		ksi_bplane[0] = N_y_M * zeta_bplane[2] - N_z_M * zeta_bplane[1];
		ksi_bplane[1] = N_z_M * zeta_bplane[0] - N_x_M * zeta_bplane[2];
		ksi_bplane[2] = N_x_M * zeta_bplane[1] - N_y_M * zeta_bplane[0];
		double ksi = sqrt(ksi_bplane[0] * ksi_bplane[0] + ksi_bplane[1] * ksi_bplane[1] + ksi_bplane[2] * ksi_bplane[2]);
		ksi_bplane[0] /= ksi;
		ksi_bplane[1] /= ksi;
		ksi_bplane[2] /= ksi;
		etta_bplane[0] = zeta_bplane[1] * ksi_bplane[2] - zeta_bplane[2] * ksi_bplane[1];
		etta_bplane[1] = zeta_bplane[2] * ksi_bplane[0] - zeta_bplane[0] * ksi_bplane[2];
		etta_bplane[2] = zeta_bplane[0] * ksi_bplane[1] - zeta_bplane[1] * ksi_bplane[0];

		// Посчитаем проекции вектора прицельной дальности на оси
		B_ksi = B_x * ksi_bplane[0] + B_y * ksi_bplane[1] + B_z * ksi_bplane[2];
		B_etta = B_x * etta_bplane[0] + B_y * etta_bplane[1] + B_z * etta_bplane[2];
		B_zeta = B_x * zeta_bplane[0] + B_y * zeta_bplane[1] + B_z * zeta_bplane[2];
		/*ksi_bplane[0] = (V_inf_plus[1] * N_z_M - V_inf_plus[2] * N_y_M) / V_inf; // делю на скорость на бесконечности, чтобы сделать её вектор единичным
		ksi_bplane[1] = (V_inf_plus[2] * N_x_M - V_inf_plus[0] * N_z_M) / V_inf;
		ksi_bplane[2] = (V_inf_plus[0] * N_y_M - V_inf_plus[0] * N_x_M) / V_inf;
		etta_bplane[0] = (V_inf_plus[1] * ksi_bplane[2] - V_inf_plus[2] * ksi_bplane[1]) / V_inf; // делю на скорость на бесконечности, чтобы сделать её вектор единичным;
		etta_bplane[1] = (V_inf_plus[2] * ksi_bplane[0] - V_inf_plus[0] * ksi_bplane[2]) / V_inf;
		etta_bplane[2] = (V_inf_plus[0] * ksi_bplane[1] - V_inf_plus[1] * ksi_bplane[0]) / V_inf;
		zeta_bplane[0] = -V_inf_plus[0] / V_inf;
		zeta_bplane[1] = -V_inf_plus[1] / V_inf;
		zeta_bplane[2] = -V_inf_plus[2] / V_inf;*/
	}



	/*a_hyp = 1 / (2 / r_M0 - V_M0 * V_M0 / mu_Moon);
	b_hyp = abs(a_hyp)*sqrt(e*e - 1);

	//дальше ищем проекции вектора прицельной дальности на эти оси
	double B_plus_x, B_plus_y, B_plus_z, B_min_x, B_min_y, B_min_z;
	B_plus_x = b_hyp * (V_inf_plus[1] * n_z - V_inf_plus[2] * n_y) / V_inf; // делю на скорость на бесконечности, чтобы сделать её вектор единичным;
	B_plus_y = b_hyp * (V_inf_plus[2] * n_x - V_inf_plus[0] * n_z) / V_inf;
	B_plus_z = b_hyp * (V_inf_plus[0] * n_y - V_inf_plus[1] * n_x) / V_inf;
	B_min_x = b_hyp * (V_inf_min[1] * n_z - V_inf_min[2] * n_y) / V_inf;
	B_min_y = b_hyp * (V_inf_min[2] * n_x - V_inf_min[0] * n_z) / V_inf;
	B_min_z = b_hyp * (V_inf_min[0] * n_y - V_inf_min[1] * n_x) / V_inf;

	B_plus_ksi = B_plus_x * ksi_bplane[0] + B_plus_y * ksi_bplane[1] + B_plus_z * ksi_bplane[2];
	B_plus_etta = B_plus_x * etta_bplane[0] + B_plus_y * etta_bplane[1] + B_plus_z * etta_bplane[2];
	B_plus_zeta = B_plus_x * zeta_bplane[0] + B_plus_y * zeta_bplane[1] + B_plus_z * zeta_bplane[2]; // должно быть 0
	B_min_ksi = B_min_x * ksi_bplane[0] + B_min_y * ksi_bplane[1] + B_min_z * ksi_bplane[2];
	B_min_etta = B_min_x * etta_bplane[0] + B_min_y * etta_bplane[1] + B_min_z * etta_bplane[2];
	B_min_zeta = B_min_x * zeta_bplane[0] + B_min_y * zeta_bplane[1] + B_min_z * zeta_bplane[2]; // должно быть 0*/
}

void Spacecraft::Collinear_impulse(double impulse)
{
	Vx_KA = Vx_KA / V_KA * (V_KA + impulse);
	Vy_KA = Vy_KA / V_KA * (V_KA + impulse);
	Vz_KA = Vz_KA / V_KA * (V_KA + impulse);
	V_KA += impulse;
}

void Spacecraft::calc_derivatives()
{
	double X = X_KA;
	double Y = Y_KA;
	double Z = Z_KA;
	double Vx = Vx_KA;
	double Vy = Vy_KA;
	double Vz = Vz_KA;
	double r = r_KA;

	dX = Vx; // производные координат
	dY = Vy;
	dZ = Vz;

	double E = 2.634 * pow(10, 10); // эпсилон для учёта сжатия Земли
	double E_Moon = 1.5 * pow(1738.14, 2) * 4902.8 * 2.033542482111609 * pow(10, -4);  // c22 = 2.240509900845084 * power (10, -5)
	double aS_x, aS_y, aS_z; // возмущения от Солнца
	double aM_x, aM_y, aM_z; // возмущения от Луны
	double aE_x, aE_y, aE_z; // возмущения от Земли
	double aJup_x, aJup_y, aJup_z;
	double aVen_x, aVen_y, aVen_z;
	double aMar_x, aMar_y, aMar_z;
	double delta1_x, delta1_y, delta1_z; // возмущения от сжатия Земли
	double delta1_x_M, delta1_y_M, delta1_z_M; // возмущения от сжатия Луны
	double aSP_x, aSP_y, aSP_z; // возмущения от Солнечного давления

	if (main_center_flag == "Earth")
	{
		if (planet_perturbations[1] == 1)
		{
			Moon_perturbations(aM_x, aM_y, aM_z); //посчитали возмущения от Луны
		}
		else
		{
			aM_x = 0;
			aM_y = 0;
			aM_z = 0;
		}
		if (planet_perturbations[2] == 1)
		{
			Sun_perturbations(aS_x, aS_y, aS_z); //посчитали возмущения от Солнца
		}
		else
		{
			aS_x = 0;
			aS_y = 0;
			aS_z = 0;
		}
		if (planet_perturbations[3] == 1)
		{
			Jupiter_perturbations(aJup_x, aJup_y, aJup_z); // посчитали возмущения от Юпитера
		}
		else
		{
			aJup_x = 0;
			aJup_y = 0;
			aJup_z = 0;
		}
		if (planet_perturbations[4] == 1)
		{
			Venus_perturbations(aVen_x, aVen_y, aVen_z); // посчитали возмущения от Венеры
		}
		else
		{
			aVen_x = 0;
			aVen_y = 0;
			aVen_z = 0;
		}
		if (planet_perturbations[5] == 1)
		{
			Mars_perturbations(aMar_x, aMar_y, aMar_z); // посчитали возмущения от Марса
		}
		else
		{
			aMar_x = 0;
			aMar_y = 0;
			aMar_z = 0;
		}

		if (earth_compression_flag == 1)
		{
			delta1_x = E / pow(r, 4) * (5 * Z * Z / (r * r) - 1) * X / r;
			delta1_y = E / pow(r, 4) * (5 * Z * Z / (r * r) - 1) * Y / r;
			delta1_z = E / pow(r, 4) * (5 * Z * Z / (r * r) - 3) * Z / r;
		}
		else
		{
			delta1_x = 0;
			delta1_y = 0;
			delta1_z = 0;
		}

		if (sun_pressure_flag == 1)
		{
			Sun_pressure(aSP_x, aSP_y, aSP_z);
		}
		else
		{
			aSP_x = 0;
			aSP_y = 0;
			aSP_z = 0;
		}

		dVx = -mu_Earth * X / (r * r * r) + aS_x + aM_x + delta1_x + aSP_x + aJup_x + aVen_x + aMar_x;/* + (Xa_x + Ya_x);*/

		dVy = -mu_Earth * Y / (r * r * r) + aS_y + aM_y + delta1_y + aSP_y + aJup_y + aVen_y + aMar_y;/* + (Xa_y + Ya_y);*/

		dVz = -mu_Earth * Z / (r * r * r) + aS_z + aM_z + delta1_z + aSP_z + aJup_z + aVen_z + aMar_z;/* + (Xa_z + Ya_z);*/
	}
	if (main_center_flag == "Moon")
	{
		/*J2000_To_SVSK(CurDat, x, y, z, X_SVSK, Y_SVSK, Z_SVSK, alpha_svsk, beta_svsk);
		r_SVSK = sqrt(pow(X_SVSK, 2) + pow(Y_SVSK, 2) + pow(Z_SVSK, 2));
		roS = sqrt(pow(XS + ro_x, 2) + pow(YS + ro_y, 2) + pow(ZS + ro_z, 2));
		roS_M = sqrt(pow(XS - XL, 2) + pow(YS - YL, 2) + pow(ZS - ZL, 2));
		Solar_Pressure.Calc_Solar_Pressure(1.1, 3000, XS, YS, ZS, (XS + ro_x), (YS + ro_y), (ZS + ro_z), roS,
			Solar_x, Solar_y, Solar_z);*/
		if (planet_perturbations[0] == 1)
		{
			Earth_perturbations(aE_x, aE_y, aE_z); //посчитали возмущения от Земли
		}
		else
		{
			aE_x = 0;
			aE_y = 0;
			aE_z = 0;
		}
		if (planet_perturbations[2] == 1)
		{
			Sun_perturbations(aS_x, aS_y, aS_z); //посчитали возмущения от Солнца
		}
		else
		{
			aS_x = 0;
			aS_y = 0;
			aS_z = 0;
		}

		if (moon_compression_flag == 1)
		{
			std::vector <double> SVSK_crds(6);
			std::vector <double> SVSK_delta1(6);
			std::vector <double> J2000_delta1(6);
			SVSK_crds = J2000_to_SVSK_3x3();
			double X_SVSK = SVSK_crds[0], Y_SVSK = SVSK_crds[1], Z_SVSK = SVSK_crds[2];

			SVSK_delta1[0] = E_Moon / pow(r, 4) * (5 * Z_SVSK * Z_SVSK / (r * r) - 1) * X_SVSK / r; // рассчитал возмущ ускор в SVSK
			SVSK_delta1[1] = E_Moon / pow(r, 4) * (5 * Z_SVSK * Z_SVSK / (r * r) - 1) * Y_SVSK / r;
			SVSK_delta1[2] = E_Moon / pow(r, 4) * (5 * Z_SVSK * Z_SVSK / (r * r) - 3) * Z_SVSK / r;
			SVSK_delta1[3] = 0;
			SVSK_delta1[4] = 0;
			SVSK_delta1[5] = 0;

			// Теперь перевожу этот вектор в J2000
			J2000_delta1 = SVSK_to_J2000_3x3(SVSK_delta1, Date_time);
			delta1_x_M = J2000_delta1[0];
			delta1_y_M = J2000_delta1[1];
			delta1_z_M = J2000_delta1[2];
		}
		else
		{
			delta1_x_M = 0;
			delta1_y_M = 0;
			delta1_z_M = 0;
		}

		if (sun_pressure_flag == 1)
		{
			Sun_pressure(aSP_x, aSP_y, aSP_z);
		}
		else
		{
			aSP_x = 0;
			aSP_y = 0;
			aSP_z = 0;
		}

		dVx = -mu_Moon * X / (r * r * r) + aS_x + aE_x + aSP_x + delta1_x_M;

		dVy = -mu_Moon * Y / (r * r * r) + aS_y + aE_y + aSP_y + delta1_y_M;

		dVz = -mu_Moon * Z / (r * r * r) + aS_z + aE_z + aSP_z + delta1_z_M;
	}
	if (main_center_flag == "Sun")
	{
		dVx = -mu_Sun * X / (r * r * r);

		dVy = -mu_Sun * Y / (r * r * r);

		dVz = -mu_Sun * Z / (r * r * r);
	}
	if (main_center_flag == "Apophis")
	{
		if (planet_perturbations[0] == 1)
		{
			Earth_perturbations(aE_x, aE_y, aE_z); //посчитали возмущения от Земли
		}
		else
		{
			aE_x = 0;
			aE_y = 0;
			aE_z = 0;
		}
		if (planet_perturbations[1] == 1)
		{
			Moon_perturbations(aM_x, aM_y, aM_z); //посчитали возмущения от Луны (тут не использую)
		}
		else
		{
			aM_x = 0;
			aM_y = 0;
			aM_z = 0;
		}
		if (planet_perturbations[2] == 1)
		{
			Sun_perturbations(aS_x, aS_y, aS_z); //посчитали возмущения от Солнца
		}
		else
		{
			aS_x = 0;
			aS_y = 0;
			aS_z = 0;
		}
		if (planet_perturbations[3] == 1)
		{
			Jupiter_perturbations(aJup_x, aJup_y, aJup_z); // посчитали возмущения от Юпитера
		}
		else
		{
			aJup_x = 0;
			aJup_y = 0;
			aJup_z = 0;
		}
		if (planet_perturbations[4] == 1)
		{
			Venus_perturbations(aVen_x, aVen_y, aVen_z); // посчитали возмущения от Венеры
		}
		else
		{
			aVen_x = 0;
			aVen_y = 0;
			aVen_z = 0;
		}
		if (planet_perturbations[5] == 1)
		{
			Mars_perturbations(aMar_x, aMar_y, aMar_z); // посчитали возмущения от Марса
		}
		else
		{
			aMar_x = 0;
			aMar_y = 0;
			aMar_z = 0;
		}

		if (sun_pressure_flag == 1)
		{
			Sun_pressure(aSP_x, aSP_y, aSP_z);
		}
		else
		{
			aSP_x = 0;
			aSP_y = 0;
			aSP_z = 0;
		}

		dVx = -mu_Apophis * X / (r * r * r) + aS_x + aE_x + aSP_x + aJup_x + aVen_x + aMar_x;

		dVy = -mu_Apophis * Y / (r * r * r) + aS_y + aE_y + aSP_y + aJup_y + aVen_y + aMar_y;

		dVz = -mu_Apophis * Z / (r * r * r) + aS_z + aE_z + aSP_z + aJup_z + aVen_z + aMar_z;
	}
}

void Spacecraft::Runge_Kutta(double h) //Р-К на один шаг
{
	double K_X = 0, K_Y = 0, K_Z = 0, K_Vx = 0, K_Vy = 0, K_Vz = 0;
	double time0 = time;
	Date_time_ms Date_time0 = Date_time;
	double X0_KA = X_KA;
	double Y0_KA = Y_KA;
	double Z0_KA = Z_KA;
	double Vx0_KA = Vx_KA;
	double Vy0_KA = Vy_KA;
	double Vz0_KA = Vz_KA;

	double integer_h, fractional_h;

	for (int k = 1; k <= 4; k++)
	{
		switch (k)
		{
		default: //при k == 1 или 2
			calc_derivatives(); //первый или второй раз производные
			K_X += k * dX; // прибавляем при k == 1 или k == 2
			K_Y += k * dY;
			K_Z += k * dZ;
			K_Vx += k * dVx;
			K_Vy += k * dVy;
			K_Vz += k * dVz;
			X_KA = X0_KA + dX * h / 2;
			Y_KA = Y0_KA + dY * h / 2;
			Z_KA = Z0_KA + dZ * h / 2;
			Vx_KA = Vx0_KA + dVx * h / 2;
			Vy_KA = Vy0_KA + dVy * h / 2;
			Vz_KA = Vz0_KA + dVz * h / 2;
			r_KA = sqrt(X_KA * X_KA + Y_KA * Y_KA + Z_KA * Z_KA);
			V_KA = sqrt(Vx_KA * Vx_KA + Vy_KA * Vy_KA + Vz_KA * Vz_KA);
			time = time0 + h / 2;
			Date_time = Date_time0;
			fractional_h = modf(h, &integer_h);
			Date_time.Date.tm_sec = Date_time.Date.tm_sec + integer_h / 2;
			Date_time.ms = Date_time.ms + 1000 * fractional_h / 2;
			Date_time.update_date_time(); // обновил все поля в структуре tm и в миллисекундах
			break;
		case 3: // k == 3
			calc_derivatives(); // третий раз производные
			K_X += 2 * dX;
			K_Y += 2 * dY;
			K_Z += 2 * dZ;
			K_Vx += 2 * dVx;
			K_Vy += 2 * dVy;
			K_Vz += 2 * dVz;
			X_KA = X0_KA + dX * h;
			Y_KA = Y0_KA + dY * h;
			Z_KA = Z0_KA + dZ * h;
			Vx_KA = Vx0_KA + dVx * h;
			Vy_KA = Vy0_KA + dVy * h;
			Vz_KA = Vz0_KA + dVz * h;
			r_KA = sqrt(X_KA * X_KA + Y_KA * Y_KA + Z_KA * Z_KA);
			V_KA = sqrt(Vx_KA * Vx_KA + Vy_KA * Vy_KA + Vz_KA * Vz_KA);
			time = time0 + h;
			Date_time = Date_time0;
			fractional_h = modf(h, &integer_h);
			Date_time.Date.tm_sec = Date_time.Date.tm_sec + integer_h;
			Date_time.ms = Date_time.ms + 1000 * fractional_h;
			Date_time.update_date_time(); // обновил все поля в структуре tm и в миллисекундах
			break;
		case 4: // k == 4
			calc_derivatives(); // четвертый раз производные
			K_X += dX;
			K_Y += dY;
			K_Z += dZ;
			K_Vx += dVx;
			K_Vy += dVy;
			K_Vz += dVz;
			X_KA = X0_KA + K_X * h / 6;
			Y_KA = Y0_KA + K_Y * h / 6;
			Z_KA = Z0_KA + K_Z * h / 6;
			Vx_KA = Vx0_KA + K_Vx * h / 6;
			Vy_KA = Vy0_KA + K_Vy * h / 6;
			Vz_KA = Vz0_KA + K_Vz * h / 6;
			r_KA = sqrt(X_KA * X_KA + Y_KA * Y_KA + Z_KA * Z_KA);
			V_KA = sqrt(Vx_KA * Vx_KA + Vy_KA * Vy_KA + Vz_KA * Vz_KA);
			time = time0 + h;
			Date_time = Date_time0;
			fractional_h = modf(h, &integer_h);
			Date_time.Date.tm_sec = Date_time.Date.tm_sec + integer_h;
			Date_time.ms = Date_time.ms + 1000 * fractional_h;
			Date_time.update_date_time(); // обновил все поля в структуре tm и в миллисекундах
			break;
		}
	}
	define_main_center(); // определяю, менять или нет главный центр
	from_Dec_to_Kepl(); // пересчитываю все параметры
	if (a > 0) { Period = 2 * Pi * sqrt(a * a * a / Grav_parameter); }
	else { Period = -1; }
	//if (distance_to_Moon() < 66000) { calc_B_plane_parameters(); }
	calc_B_plane_parameters();
}

Spacecraft& Spacecraft::operator=(const Spacecraft& right)
{
	//проверка на самоприсваивание
	if (this == &right)
	{
		return *this;
	}
	X_KA = right.X_KA;
	Y_KA = right.Y_KA;
	Z_KA = right.Z_KA;
	Vx_KA = right.Vx_KA;
	Vy_KA = right.Vy_KA;
	Vz_KA = right.Vz_KA;
	dX = right.dX;
	dY = right.dY;
	dZ = right.dZ;
	dVx = right.dVx;
	dVy = right.dVy;
	dVz = right.dVz;
	Date_time = right.Date_time; //дата
	EA = right.EA;
	MA = right.MA;
	H = right.H;
	Grav_parameter = right.Grav_parameter;
	earth_compression_flag = right.earth_compression_flag;
	moon_compression_flag = right.moon_compression_flag;
	sun_pressure_flag = right.sun_pressure_flag;
	planet_perturbations = right.planet_perturbations;

	main_center_flag = right.main_center_flag;
	ecc = right.ecc; inc = right.inc; RAAN = right.RAAN; p = right.p; a = right.a; omega = right.omega; u = right.u;
	theta = right.theta; tau = right.tau;
	time = right.time;
	Period = right.Period;
	r_KA = right.r_KA;
	V_KA = right.V_KA;
	rad_KA = right.rad_KA; //радиус КА
	mass = right.mass; // масса КА
	Area = right.Area; //площадь для учёта солнечного давления
	koeff_otraj = right.koeff_otraj; // коэффициент для солнечного давления 1 < k < 1.44

	//ниже начальные параметры для перелёта к Луне (с индексом 0 - опорной орбиты; с индексом 1 - перелйтной орбиты; с индексом targ - точки попадания в Луну)
	fi = right.fi, lambda = right.lambda; // координаты космодрома
	inc0 = right.inc0; r0 = right.r0; a0 = right.a0; ecc0 = right.ecc0; RAAN0 = right.RAAN0;
	launch_time = right.launch_time; //время запуска КА с Земли в секундах
	Start_date = right.Start_date; // старт с Земли
	inc1 = right.inc1; a1 = right.a1; ecc1 = right.ecc1; RAAN1 = right.RAAN1; omega1 = right.omega1;
	tau1_date = right.tau1_date; // дата подачи импульса
	u_depart = right.u_depart; /*аргумент широты при подаче импульса*/ delta_t_depart = right.delta_t_depart; /*время от theta=0 до точки отлёта от Земли*/
	fly_away_impulse = right.fly_away_impulse; // отлетный импульс
	V_inf = right.V_inf; // модуль скорости на бесконечности

	return *this;
}