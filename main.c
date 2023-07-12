/*
* A C code version of gravity turn from MATLAB OrbitalMechanics Examples 
* "Example_11_03"
* Implementation of Runge-Kutta-Fehlberg method is written by ChatGPT
* Equations and MATLAB Code comes from the book Orbital Mechanics For Engineering Students (third edition)
* + Output File 
* + Configurable File
* 07/12/2023
* --- Xenocide ---
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265359
#define MAX_LINE_LENGTH 100

double deg_to_rad(double deg){
	return deg * (PI / 180.0);
}

void read_config_file(const char *filename, double *variables){
	FILE *file = fopen(filename, "r");
	if(file == NULL){
		printf("Configuration file not found.\n");
		return;
	}
	char line[256];
	while (fgets(line, sizeof(line), file)){
		if(line[0] == '#' || line[0] == '\n'){
			continue; // Ignore comment line and empty lines
		}
		char variable[256];
		double value;
		if(sscanf(line, "%s = %lf", variable, &value) == 2){
			if(strcmp(variable, "diam") == 0){
				variables[0] = value;
			}else if(strcmp(variable, "A") == 0){
				variables[1] = value;
			}else if(strcmp(variable, "CD") == 0){
				variables[2] = value;
			}else if(strcmp(variable, "m0") == 0){
				variables[3] = value;
			}else if(strcmp(variable, "n") == 0){
				variables[4] = value;
			}else if(strcmp(variable, "T2W") == 0){
				variables[5] = value;
			}else if(strcmp(variable, "Isp") == 0){
				variables[6] = value;
			}else if(strcmp(variable, "mfinal") == 0){
				variables[7] = value;
			}else if(strcmp(variable, "Thrust") == 0){
				variables[8] = value;
			}else if(strcmp(variable, "m_dot") == 0){
				variables[9] = value;
			}else if(strcmp(variable, "mprop") == 0){
				variables[10] = value;
			}else if(strcmp(variable, "tburn") == 0){
				variables[11] = value;
			}else if(strcmp(variable, "hturn") == 0){
				variables[12] = value;
			}
		}
	}
	fclose(file);
}

void rkf45(double (*rates)(double, double[], double[], const double *), double tspan[], double y0[], double t[], double f[][6], int length, const double *variables){
	double h = (tspan[1] - tspan[0]) / (length - 1);
	double t_current = tspan[0];
	double y[6];
	
	for(int i = 0; i < 6; i++){
		y[i] = y0[i];
		f[0][i] = y0[i];
	}
	
	for(int i = 1; i < length; i++){
		double k1[6], k2[6], k3[6], k4[6], k5[6], k6[6], y_temp[6], error[6];

		// Calculate k1
		rates(t_current, y, k1, variables);

		// Calculate k2
		for(int j = 0; j < 6; j++){
			y_temp[j] = y[j] + h * (k1[j] / 4.0);
		}
		rates(t_current + h / 4.0, y_temp, k2, variables);

		// Calculate k3
		for(int j = 0; j < 6; j++){
			y_temp[j] = y[j] + h * (3.0 * k1[j] / 32.0 + 9.0 * k2[j] / 32.0);
		}
		rates(t_current + 3.0 * h / 8.0, y_temp, k3, variables);

		// Calculate k4
		for(int j = 0; j < 6; j++){
			y_temp[j] = y[j] + h * (1932.0 * k1[j] / 2197.0 - 7200.0 * k2[j] / 2197.0 + 7296.0 * k3[j] / 2197.0);
		}
		rates(t_current + 12.0 * h / 13.0, y_temp, k4, variables);

		// Calculate k5
		for(int j = 0; j < 6; j++){
			y_temp[j] = y[j] + h * (439.0 * k1[j] / 216.0 - 8.0 * k2[j] + 3680.0 * k3[j] / 513.0 - 845.0 * k4[j] / 4104.0);
		}
		rates(t_current + h, y_temp, k5, variables);

		// Calculate k6
		for(int j = 0; j < 6; j++){
			y_temp[j] = y[j] + h * (-8.0 * k1[j] / 27.0 + 2.0 * k2[j] - 3544.0 * k3[j] / 2565.0 + 1859.0 * k4[j] / 4104.0 - 11.0 * k5[j] / 40.0);
		}
		rates(t_current + h / 2.0, y_temp, k6, variables);

		// Solve at t = t_current + h
		for(int j = 0; j < 6; j++){
			y[j] = y[j] + h * (16.0 * k1[j] / 135.0 + 6656.0 * k3[j] / 12825.0 + 28561.0 * k4[j] / 56430.0 - 9.0 * k5[j] / 50.0 + 2.0 * k6[j] / 55.0);
		}
		t_current += h;
		// Store solution at t = t_current
		t[i] = t_current;
		for(int j = 0; j < 6; j++){
			f[i][j] = y[j];
		}
	}
}

double rates(double t, double y[], double dydt[], const double *variables){
	double deg = PI / 180.0;
	double g0 = 9.81;
	double Re = 6378e3;
	double hscale = 7.5e3;
	double rho0 = 1.225;
	double diam = variables[0];
	double A = variables[1];
	double CD = variables[2];
	double m0 = variables[3];
	double n = variables[4];
	double T2W = variables[5];
	double Isp = variables[6];
	double mfinal = variables[7];
	double Thrust = variables[8];
	double m_dot = variables[9];
	double mprop = variables[10];
	double tburn = variables[11];
	double hturn = variables[12];
	double h = y[3];
	double v = y[0];
	double gamma = y[1];
	double m;
	if(t < tburn){
		m = m0 - m_dot * t;
	}else{
		m = m0 - m_dot * tburn;
	}
	double g = g0 / pow(1 + h / Re, 2);
	double rho = rho0 * exp(-h / hscale);
	double D = 0.5 * rho * pow(v, 2) * A * CD;
	double v_dot, gamma_dot, x_dot, h_dot, vG_dot, vD_dot;
	if(h <= hturn){
		gamma_dot = 0;
		v_dot = (Thrust / m) - (D / m) - g;
		x_dot = 0;
		h_dot = v;
		vG_dot = -g;
	}else{
		v_dot = (Thrust / m) - (D / m) - g * sin(gamma);
		gamma_dot = (-1 / v) * (g - pow(v, 2) / (Re + h)) * cos(gamma);
		x_dot = (Re / (Re + h)) * v * cos(gamma);
		h_dot = v * sin(gamma);
		vG_dot = -g * sin(gamma);
	}
	vD_dot = -D / m;
	dydt[0] = v_dot;
	dydt[1] = gamma_dot;
	dydt[2] = x_dot;
	dydt[3] = h_dot;
	dydt[4] = vD_dot;
	dydt[5] = vG_dot;
	return 0;
}

void output(double t[], double f[][6], int length, double gamma0, double hturn, double tburn){
	double deg = PI / 180.0;
	double rho0 = 1.225;
	double hscale = 7.5e3;
	double q[1000];

	FILE *dataFile = fopen("output.txt", "w");
	if(dataFile == NULL){
		printf("output.txt file not found.\n");
		return;
	}// Write variable names and units in the first line of the data file
	fprintf(dataFile, "Time (s), Velocity (m/s), Flight Path Angle (deg), Downrange Distance (km), Altitude (km), Drag Loss (m/s), Gravity Loss (m/s)\n");

	for(int i = 0; i < length; i++){
		double h = f[i][3] * 1e-3;
		double Rho = rho0 * exp(-h * 1000.0 / hscale);
		q[i] = 0.5 * Rho * pow(f[i][0] * 1e-3, 2);
		// Write data to the output file
		fprintf(dataFile, "%.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f\n", t[i], f[i][0], f[i][1] / deg, f[i][2] * 1e-3, f[i][3] * 1e-3, -f[i][4] * 1e-3, -f[i][5] * 1e-3);
	}
	fclose(dataFile);
	printf("\n\n -----------------------------------\n");
	printf("\n Initial flight path angle = %10g deg ", gamma0 / deg);
	printf("\n Pitchover altitude          = %10g m   ", hturn);
	printf("\n Burn time                        = %10g s   ", tburn);
	printf("\n Final speed                     = %10g km/s", f[length - 1][0]/1000);
	printf("\n Final flight path angle   = %10g deg ", f[length - 1][1]*180/PI);
	printf("\n Altitude                            = %10g km  ", f[length - 1][3]/1000);
	printf("\n Downrange distance      = %10g km  ", f[length - 1][2]/1000);
	printf("\n Drag loss                         = %10g km/s", -f[length - 1][4]/1000);
	printf("\n Gravity loss                     = %10g km/s", -f[length - 1][5]/1000);
	printf("\n\n -----------------------------------\n");
}

int main(){
	double deg = PI / 180.0;
	double g0 = 9.81;
	double Re = 6378e3;
	double hscale = 7.5e3;
	double rho0 = 1.225;
	double variables[13]; 
	
	// Read configuration file and update variables
	read_config_file("files/gturn.cfg", variables);
	
	double t0 = 0;
	double tf = variables[11];
	double tspan[2] = {t0, tf};
	double v0 = 0;
	double gamma0 = 89.85 * deg;
	double x0 = 0;
	double h0 = 0;
	double vD0 = 0;
	double vG0 = 0;
	double f0[6] = {v0, gamma0, x0, h0, vD0, vG0};
	double t[1000];									// Adjust the size as needed
	double f[1000][6];								// Adjust the size as needed
	rkf45(rates, tspan, f0, t, f, 1000, variables); // Runge-Kutta-Fehlberg method 'rkf45'
	output(t, f, 1000, gamma0, variables[12], variables[11]);
	return 0;
}
