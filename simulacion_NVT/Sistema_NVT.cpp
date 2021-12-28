#include "Sistema_NVT.h"

//constructores

//destructor

//Algoritmo gear
void Sistema_NVT::gear4(float dt, double deltax[], double deltay[], double deltaz[])
{
	//fijamos la Ec a 0
	this->E_cin = 0;
	double deltas = 0;//diferencia de aceleraciones en el foco

	//primero generamos las previsiones de los vectores
	for (int i = 0; i < N; i++)
	{
		rx[i] += vx[i] * dt + 0.5 * ax[i] * dt * dt + 0.167 * dax[i] * dt * dt * dt + 0.04167 * d2ax[i] * dt * dt * dt * dt;
		ry[i] += vy[i] * dt + 0.5 * ay[i] * dt * dt + 0.167 * day[i] * dt * dt * dt + 0.04167 * d2ay[i] * dt * dt * dt * dt;
		rz[i] += vz[i] * dt + 0.5 * az[i] * dt * dt + 0.167 * daz[i] * dt * dt * dt + 0.04167 * d2az[i] * dt * dt * dt * dt;
		
		vx[i] += ax[i] * dt + 0.5 * dax[i] * dt * dt + 0.167 * d2ax[i] * dt * dt * dt;
		vy[i] += ay[i] * dt + 0.5 * day[i] * dt * dt + 0.167 * d2ay[i] * dt * dt * dt;
		vz[i] += az[i] * dt + 0.5 * daz[i] * dt * dt + 0.167 * d2az[i] * dt * dt * dt;

		ax[i] += dax[i] * dt + 0.5 * d2ax[i] * dt * dt;
		ay[i] += day[i] * dt + 0.5 * d2ay[i] * dt * dt;
		az[i] += daz[i] * dt + 0.5 * d2az[i] * dt * dt;
		//añadimos las previsiones al vector delta
		deltax[i] = ax[i];
		deltay[i] = ay[i];
		deltaz[i] = az[i];

		dax[i] += d2ax[i] * dt;
		day[i] += d2ay[i] * dt;
		daz[i] += d2az[i] * dt;

		//previsión para s
		s += ds * dt + 0.5 * d2s * dt * dt + 0.167 * d3s * dt * dt * dt + 0.04167 * d4s * dt * dt * dt * dt;
		ds += d2s * dt + 0.5 * d3s * dt * dt + 0.167 * d4s * dt * dt * dt;
		d2s += d3s * dt + 0.5 * d4s * dt * dt;
		d3s += d4s * dt;

		deltas = d2s;
	}
	//llamamos al potencial para calcular las aceleraciones que corresponderian a estas posiciones
	lennard_jones(N, rc, rc2, L, V, &E_p, &dphi, &d2phi, rx, ry, rz, ax, ay, az);
	//aceleracion del foco
	d2s = 0;
	for (int i = 0; i < N; i++)
	{
		d2s +=(rx[i] * rx[i] + ry[i] * ry[i] + rz[i] * rz[i]);
	}
	d2s *= s;
	d2s -= T / s;//falta un factor v que no recuerdo a que se corresponde
	d2s *= 1 / MT;
	//diferencia de aceleraciones del foco
	deltas = d2s - deltas;

	//hay que ajustar las aceleraciones con el foco termico
	for (int i = 0; i < N; i++)
	{
		ax[i] *= 1 / (s * s);
		ay[i] *= 1 / (s * s);
		az[i] *= 1 / (s * s);

		ax[i] -= 2 * ds * vx[i] / s;
		ay[i] -= 2 * ds * vy[i] / s;
		az[i] -= 2 * ds * vz[i] / s;
	}

	//calculamos el vector delta para corregir las previsiones
	for (int i = 0; i < N; i++)
	{
		deltax[i] = ax[i] - deltax[i];
		deltay[i] = ay[i] - deltay[i];
		deltaz[i] = az[i] - deltaz[i];
	}
	//por ultimo calculamos las correcciones al resto de variables
	for (int i = 0; i < N; i++)
	{
		rx[i] += 0.079167 * dt * dt * deltax[i];
		ry[i] += 0.079167 * dt * dt * deltay[i];
		rz[i] += 0.079167 * dt * dt * deltaz[i];

		vx[i] += 0.365 * dt * deltax[i];
		vy[i] += 0.365 * dt * deltay[i];
		vz[i] += 0.365 * dt * deltaz[i];
		//recalculamos la energia cinetica
		E_cin += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];

		dax[i] += 1.5 * deltax[i] / dt;
		day[i] += 1.5 * deltay[i] / dt;
		daz[i] += 1.5 * deltaz[i] / dt;

		d2ax[i] += deltax[i] / (dt * dt);
		d2ay[i] += deltay[i] / (dt * dt);
		d2az[i] += deltaz[i] / (dt * dt);

		//correccion de las variables del foco
		s += 0.079167 * dt * dt * deltas;
		ds += 0.365 * dt * deltas;
		d3s += 1.5 * deltas / dt;
		d4s += deltas / (dt * dt);
	}

	//factor 0.5 de la cinetica y hay que multiplicar por s para obtener las velocidades reales
	E_cin *= 0.5*s*s;
}

void Sistema_NVT::stab_gear(int npasos, float dt, int save_rate, std::string energy_filename, std::string vel_filename, std::string process_name)
{
	std::ofstream energy_output(energy_filename, std::ios::out);
	std::ofstream velocity_output(vel_filename, std::ios::binary);
	progress_bar bar(process_name, npasos);
	int array_size = N * sizeof(double);
	//definimos la variable delta que hace falta para el gear, la creo como variable externa para no estar creando y destruyendo 3*N doubles en cada paso del gear
	double *deltax = new double[this->N];
	double *deltay = new double[this->N];
	double *deltaz = new double[this->N];

	for (int j = 0; j < npasos; j++)
	{
		gear4(dt, deltax, deltay, deltaz);

		if (j % save_rate == 0)
		{
			energy_output << std::fixed << std::setw(10) << std::setprecision(6) << E_p << " " << E_cin << std::endl;
			velocity_output.write((char*)vx, array_size);
			velocity_output.write((char*)vy, array_size);
			velocity_output.write((char*)vz, array_size);
		}
		bar.update(j);
	}
	energy_output.close();
	velocity_output.close();

	//borramos las variables auxiliares delta
	delete[] deltax;
	delete[] deltay;
	delete[] deltaz;
}

void Sistema_NVT::run_gear(int n_pasos, float dt, int save_rate, std::string results_file, std::string energy_filename, std::string bin_filename, std::string process_name)
{
	std::ofstream energy_output(energy_filename, std::ios::out);
	std::ofstream bin_output(bin_filename, std::ios::binary);
	progress_bar bar(process_name, n_pasos);
	int array_size = N * sizeof(double);

	//definimos la variable delta que hace falta para el gear
	double* deltax = new double[this->N];
	double* deltay = new double[this->N];
	double* deltaz = new double[this->N];
	
	for (int j = 0; j < n_pasos; j++)
	{
		gear4(dt, deltax, deltay, deltaz);
		//calculamos las energias medias y la media de las derivadas
		

		if (j % save_rate == 0)
		{
			energy_output << std::fixed << std::setw(10) << std::setprecision(6) << E_p << " " << E_cin << " " << dphi << " " << d2phi << std::endl;
			bin_output.write((char*)rx, array_size);
			bin_output.write((char*)ry, array_size);
			bin_output.write((char*)rz, array_size);
			bin_output.write((char*)vx, array_size);
			bin_output.write((char*)vy, array_size);
			bin_output.write((char*)vz, array_size);
		}
		bar.update(j);
	}
	

	//por ultimo usamos esto para calcular los parametros del sistema en unidades reducidas
	
	//A partir de estos resultados calculamos el resto:
	

	//calcular alpha_E por el segundo metodo y comprobar que dan cosas similares

	energy_output.close();
	bin_output.close();

	//escribimos los resultadosa un archivo de forma que sea facil de abrir con numpy e interpretar a simple vista

	//borramos las variables auxiliares delta;
	delete[] deltax;
	delete[] deltay;
	delete[] deltaz;
	
}
