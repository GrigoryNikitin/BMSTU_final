#pragma once
#include "pch.h"


class Spacecraft
{
public:
	Spacecraft();
	Spacecraft(std::string main_center, std::string a, int ms);
	~Spacecraft();
	void printParameters(std::ofstream& fout);
	void printKepler(std::ofstream& fout);
	void printBPlane(std::ofstream& fout);
	void printSVSK(std::ofstream& fout);
	void printGCS(std::ofstream& fout);
	void printParameters_relatively_to_Moon(std::ofstream& fout);
	void printKepler_relatively_to_Moon(std::ofstream& fout);
	void printSVSK_relatively_to_Moon(std::ofstream& fout);
	void printSISK_relatively_to_Moon(std::ofstream& fout);
	std::vector <double> get_parameters();
	std::vector <double> get_b_plane_parameters();
	void set_parameters(std::vector <double> Parameters);
	void set_Date(std::string a, int ms);
	void set_earth_compression(int a);
	void set_moon_compression(int a);
	void set_sun_pressure(int a);
	void set_planet_perturbations(std::vector <bool> a);
	void from_Dec_to_Kepl();
	void from_Dec_to_Kepl(std::string chosen_center_flag, double& a, double& ecc, double& inc, double& RAAN, double& omega, double& tau, double& theta); // перегрузка для удобного использования без учёта главного центра
	void from_Kepl_to_Dec();
	void from_Kepl_Moon_to_Dec_J2000();
	std::vector <double> J2000_to_SISK_6x6();
	std::vector <double> J2000_to_GCS_6x6();
	matrix J2000_to_SVSK_matrix_3x3(Date_time_ms Date);
	std::vector <double> J2000_to_SVSK_3x3();
	std::vector <double> SVSK_to_J2000_3x3(std::vector <double> SVSK, Date_time_ms Date_time); // для расчёта по конкретному вектору и дате
	std::vector <double> J2000_to_SVSK_6x6();
	std::vector <double> SVSK_to_J2000_6x6();
	std::vector <double> SVSK_to_J2000_6x6(std::vector <double> SVSK, Date_time_ms Date_time); // перегрузка для расчёта по конкретному вектору и дате
	double SVSK_inclination();
	std::vector <double> get_planet_parameters(const char Date_time[], const char target_body[], const char reference_frame[], const char observer_body[]);
	void Sun_perturbations(double& aS_x, double& aS_y, double& aS_z);
	void Moon_perturbations(double& aM_x, double& aM_y, double& aM_z);
	void Earth_perturbations(double& aE_x, double& aE_y, double& aE_z);
	void Jupiter_perturbations(double& aJup_x, double& aJup_y, double& aJup_z);
	void Venus_perturbations(double& aVen_x, double& aVen_y, double& aVen_z);
	void Mars_perturbations(double& aMar_x, double& aMar_y, double& aMar_z);
	void Sun_pressure(double& aSP_x, double& aSP_y, double& aSP_z);
	double distance_to_Moon();
	double distance_to_Earth();
	void define_main_center();
	double alpfa0(tm D1);
	void Finding_RAAN0_north(Date_time_ms Start_Date, double fi, double lambda, double& RAAN0_north);
	void Finding_RAAN0_south(Date_time_ms Start_Date, double fi, double lambda, double& RAAN0_south);
	void Moon_flight_without_perturbations_north(double flight_time /*время перелёта в сутках*/, Date_time_ms& approach_date);
	void Moon_flight_without_perturbations_south(double flight_time /*время перелёта в сутках*/, Date_time_ms& approach_date);
	void calc_B_plane_parameters();
	void Collinear_impulse(double impulse);
	void calc_derivatives();
	void Runge_Kutta(double h);

	Spacecraft& operator=(const Spacecraft& right);

	Date_time_ms Date_time; // дата моя с миллисекундами
	std::string main_center_flag;
	double ecc, inc, RAAN, p, a, omega, u, theta, tau;
	double time, Period, r_KA, V_KA, trajectory_slope;
	double rad_KA = 2; //радиус КА
	double mass = 3000; // масса КА
	double Area = Pi * rad_KA * rad_KA; //площадь для учёта солнечного давления
	double koeff_otraj = 1; // коэффициент для солнечного давления 1 < k < 1.44

	//ниже начальные параметры для перелёта к Луне (с индексом 0 - опорной орбиты; с индексом 1 - перелйтной орбиты; с индексом targ - точки попадания в Луну)
	double fi = 51.88 * toRad, lambda = 128.33 * toRad; // координаты космодрома
	double inc0 = 52 * toRad, r0 = 6571, a0 = 6571, ecc0 = 0.0000000001, RAAN0;
	double launch_time = 1200; //время запуска КА с Земли в секундах
	Date_time_ms Start_date; // старт с Земли
	double inc1, a1, ecc1, RAAN1, omega1;
	Date_time_ms tau1_date; // дата подачи импульса
	double u_depart /*аргумент широты при подаче импульса*/, delta_t_depart /*время от theta=0 до точки отлёта от Земли*/;
	double fly_away_impulse; // отлетный импульс
	double V_inf, V_inf0; // модуль скорости на бесконечности

protected:
	double X_KA, Y_KA, Z_KA, Vx_KA, Vy_KA, Vz_KA; //параметры

	double dX, dY, dZ, dVx, dVy, dVz; //производные
	double EA, MA, H;
	double Grav_parameter;
	bool earth_compression_flag, moon_compression_flag, sun_pressure_flag;
	std::vector <bool> planet_perturbations; // массив из 0 и 1, чтобы учитывать возмущ от нужных планет
	/*Земля, Луна, Солнце, Юпитер, Венера, Марс*/

	// Ниже параметры для расчёта Картинной плоскости и сами её параметры
	double n, n_x, n_y, n_z; //вектор нормали к плоскости орбиты КА
	double e, e_x, e_y, e_z; //вектор эксцентриситета
	//double a_hyp, b_hyp, alpha_hyp; //большая и малая полуоси гиперболы и альфа
	//std::vector <double> V_inf_plus, V_inf_min; // вектора скорости на бесконечности относительно Луны
	bool b_plane_first_time = 1;
	std::vector <double> ksi_bplane, etta_bplane, zeta_bplane; // кси и этта оси Картинной плоскости
	//double B_plus_ksi, B_plus_etta, B_plus_zeta, B_min_ksi, B_min_etta, B_min_zeta;
	double B_ksi, B_etta, B_zeta;

};


