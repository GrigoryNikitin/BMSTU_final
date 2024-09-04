#pragma once

const double Pi = 3.14159265359;
const double Ae = 149597870.700; // астр. единица
const double Moon_synodic_period = 29.5306; // Синодический период Луны в сутках
const double mu_Earth = 398600.44; //398600.44
const double omega_earth = 7.2921158553 * pow(10, -5); // скорость вращени¤ «емли рад/с
const double mu_Moon = 4902.80; //4902.80
const double mu_Sun = 132712440018; //132712440018
const double mu_Apophis = 2 * pow(10, -9); // 2*pow(10, -9) мю апофиса
const double mu_Jupiter = 126686534; // 126686534
const double mu_Venus = 324859; // 324859
const double mu_Mars = 42828; // 42828 
const double r_per = 0.5; //0.5
const double r_apo = 0.54; //0.54
const double X0 = 1837.1;//0.5;
const double Y0 = 0;
const double Z0 = 0;
const double Vx0 = 0;
const double Vy0 = 0;
const double Vz0 = 1.633637485; //sqrt(2 * mu_Earth / r_per * r_apo / (r_per + r_apo));
const double h = 20; //20
const double toDeg = 180 / Pi;
const double toRad = Pi / 180;