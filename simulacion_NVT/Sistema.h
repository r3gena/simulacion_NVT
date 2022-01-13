#pragma once
#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include <string>
#include "potenciales.h"
#include "progress_bar.h"

class Sistema
{
public:
	//constantes del sistema
	int N;
	float rho, V, L, rc, rc2;
	double E_t, E_cin, E_p, dphi, d2phi;

	Sistema(int num_k, int N, float rho, float E_t);//constructor de la clase, inicializa variables del sistema y crea los vectores del tama�o indicado
	Sistema(std::string ascii_filename);//extrae los parametros del archivo indicado y busca el archivo binario para exxtraer los vectores
	Sistema();
	~Sistema();//destructor, se asegura de eliminar los vectores para no tener fugas de memoria

	void imprimir_energias();
	void exportar_datos(std::string ascii_filename, std::string bin_filename);//exporta los parametros del sistema (energias, densidad, numero de particulas...) a un archivo de texto
	void exportar_datos(std::string ascii_filename, std::string bin_filename, float dt, int n_pasos);
	void exportar_configuracion(std::string bin_filename);//exporta los vectores r, v y a a un archivo binario

	//algoritmos de evolucion temporal
	void stab_verlet(int n_pasos, float dt, int save_rate, std::string energy_filename, std::string vel_filename, std::string process_name);
	void run_verlet(int n_pasos, float dt, int save_rate, std::string results_file, std::string energy_filename, std::string bin_filename, std::string process_name);
	//revisar, no se si en el run quiero exportar un archivo de energias o si quiero la energia media de cada proceso
	void corregir_energia();

protected:

	void colocar_red(int num_k, double a);

	void velocity_verlet(float dt);

	//posiciones, velocidades y aceleraciones, creados como pionters porque no se pueden definir arrays sin tama�o
	double *rx, *ry, *rz;
	double *vx, *vy, *vz;
	double *ax, *ay, *az;

};
