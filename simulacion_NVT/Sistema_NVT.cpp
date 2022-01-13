#include "Sistema_NVT.h"

//constructor
Sistema_NVT::Sistema_NVT(std::string filename)
{
	//abro el archivo
	std::string bin_filename;
	std::ifstream inputfile(filename, std::ios::in);
	inputfile >> this->N >> this->rho >> this->L >> this->rc;
  inputfile >> this->T >> this-> MT;
  inputfile >> bin_filename;
	inputfile.close();

	this->V=L*L*L;
	this->rc2=rc*rc;

	//inicializo las posiciones
	this->rx = new double[N];
	this->ry = new double[N];
	this->rz = new double[N];
	this->vx = new double[N];
	this->vy = new double[N];
	this->vz = new double[N];
	this->ax = new double[N];
	this->ay = new double[N];
	this->az = new double[N];
	this->dax = new double[N];
	this->day = new double[N];
	this->daz = new double[N];
	this->d2ax = new double[N];
	this->d2ay = new double[N];
	this->d2az = new double[N];

	//hay que inicializar tambien las variables del foco
	s=1;
	ds=0, d2s=0, d3s=0, d4s=0;

	std::ifstream binfile(bin_filename, std::ios::binary);
	binfile.read((char*)rx, N*sizeof(double));
	binfile.read((char*)ry, N*sizeof(double));
	binfile.read((char*)rz, N*sizeof(double));
	binfile.read((char*)vx, N*sizeof(double));
	binfile.read((char*)vy, N*sizeof(double));
	binfile.read((char*)vz, N*sizeof(double));
	binfile.read((char*)ax, N*sizeof(double));
	binfile.read((char*)ay, N*sizeof(double));
	binfile.read((char*)az, N*sizeof(double));
	binfile.close();

	//ponemos a 0 las derivadas extra de a
	for(int i=0; i < N; i++)
	{
		dax[i]=0;
		day[i]=0;
		daz[i]=0;
		d2ax[i]=0;
		d2ay[i]=0;
		d2az[i]=0;
	}
	lennard_jones(N, rc, rc2, L, V, &E_p, &dphi, &d2phi, rx, ry, rz, ax, ay, az);
}
//destructor
Sistema_NVT::~Sistema_NVT()
{
	delete[] dax;//creo que llega con poner las variables que no estan en el destructor de la clase padre
	delete[] day;
	delete[] daz;
	delete[] d2ax;
	delete[] d2ay;
	delete[] d2az;
}

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
		//a�adimos las previsiones al vector delta
		deltax[i] = ax[i];
		deltay[i] = ay[i];
		deltaz[i] = az[i];

		dax[i] += d2ax[i] * dt;
		day[i] += d2ay[i] * dt;
		daz[i] += d2az[i] * dt;

	}
	//previsi�n para s
	s += ds * dt + 0.5 * d2s * dt * dt + 0.167 * d3s * dt * dt * dt + 0.04167 * d4s * dt * dt * dt * dt;
	ds += d2s * dt + 0.5 * d3s * dt * dt + 0.167 * d4s * dt * dt * dt;
	d2s += d3s * dt + 0.5 * d4s * dt * dt;
	d3s += d4s * dt;

	deltas = d2s;

	//llamamos al potencial para calcular las aceleraciones que corresponderian a estas posiciones
	lennard_jones(N, rc, rc2, L, V, &E_p, &dphi, &d2phi, rx, ry, rz, ax, ay, az);
	//aceleracion del foco
	d2s = 0;
	for (int i = 0; i < N; i++)
	{
		d2s = d2s + (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
	}
	d2s = d2s * s;
	d2s = d2s - (T * (3*N-3)/ s);//falta un factor v que no recuerdo a que se corresponde
	d2s = d2s / MT;
	//diferencia de aceleraciones del foco
	deltas = d2s - deltas;

	//hay que ajustar las aceleraciones con el foco termico
	for (int i = 0; i < N; i++)
	{
		ax[i] = ax[i] / (s * s);
		ay[i] = ay[i] / (s * s);
		az[i] = az[i] / (s * s);

		ax[i] = ax[i] + 2 * ds * vx[i] / s;
		ay[i] = ay[i] + 2 * ds * vy[i] / s;
		az[i] = az[i] + 2 * ds * vz[i] / s;
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
		this->E_cin = E_cin + vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];

		dax[i] += 1.5 * deltax[i] / dt;
		day[i] += 1.5 * deltay[i] / dt;
		daz[i] += 1.5 * deltaz[i] / dt;

		d2ax[i] += deltax[i] / (dt * dt);
		d2ay[i] += deltay[i] / (dt * dt);
		d2az[i] += deltaz[i] / (dt * dt);

	}
	//correccion de las variables del foco
	s += 0.079167 * dt * dt * deltas;
	ds += 0.365 * dt * deltas;
	d3s += 1.5 * deltas / dt;
	d4s += deltas / (dt * dt);

	//factor 0.5 de la cinetica y hay que multiplicar por s para obtener las velocidades reales
	E_cin = E_cin*0.5*s*s;
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
			velocity_output.write((char*)&s, sizeof(double));//no se si esta es la mejor forma pero para tener mayor velocidad guardo v y s por separado para construir la velocidad correcta en el analisis posterior
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

void Sistema_NVT::run_gear(int n_pasos, float dt, int save_rate, std::string results_filename, std::string energy_filename, std::string bin_filename, std::string process_name)
{
	std::ofstream energy_output(energy_filename, std::ios::out);
	std::ofstream bin_output(bin_filename, std::ios::binary);
	progress_bar bar(process_name, n_pasos);
	int array_size = N * sizeof(double);

	//definimos la variable delta que hace falta para el gear
	double* deltax = new double[this->N];
	double* deltay = new double[this->N];
	double* deltaz = new double[this->N];

	//medias que voy a calcular
	double mEc=0;//para comprobar la temperatura
	double mphi=0;//potencial
	double mdphi=0, md2phi=0;//primera y segunda derivada
	double mphi2=0, mdphi2=0;//potencial y primera derivada al cuadrado
	double mphidphi=0;//potencial por su primera derivada

	for (int j = 0; j < n_pasos; j++)
	{
		gear4(dt, deltax, deltay, deltaz);
		//calculamos las energias medias y la media de las derivadas
		mEc=mEc+E_cin;
		mphi=mphi+E_p;
		mdphi=mdphi+dphi;
		mphi2=mphi2+E_p*E_p;
		mdphi2=mdphi2+dphi*dphi;
		md2phi=md2phi+d2phi;
		mphidphi=mphidphi+E_p*dphi;

		if (j % save_rate == 0)
		{
			energy_output << std::fixed << std::setw(10) << std::setprecision(6) << E_p << " " << E_cin << " " << dphi << " " << d2phi << std::endl;
			bin_output.write((char*)rx, array_size);
			bin_output.write((char*)ry, array_size);
			bin_output.write((char*)rz, array_size);
			bin_output.write((char*)vx, array_size);
			bin_output.write((char*)vy, array_size);
			bin_output.write((char*)vz, array_size);
			bin_output.write((char*)&s, sizeof(double));//no se si esta es la mejor forma pero para tener mayor velocidad guardo v y s por separado para construir la velocidad correcta en el analisis posterior
		}
		bar.update(j);
	}
	mEc=(mEc/n_pasos);
	mphi=(mphi/n_pasos);
	mdphi=(mdphi/n_pasos);
	mphi2=(mphi2/n_pasos);
	mdphi2=(mdphi2/n_pasos);
	md2phi=(md2phi/n_pasos);
	mphidphi=(mphidphi/n_pasos);

	//por ultimo usamos esto para calcular los parametros del sistema en unidades reducidas
	double temp = 2 * mEc / (3 * N - 3);//no tengo claro si esto es valido aqui pero para comprobar voy a meterlo
	double p=(N*T/V)-mdphi;
	double Cv=((3*N-3)/2)+(mphi2-mphi*mphi)/(T*T);
	double gamma=(V/Cv)*((N/V)+(mphi*mdphi-mphidphi)/(T*T));
	double inv_kT=(N*T/V)+V*md2phi-(V*(mdphi2-(mdphi*mdphi))/T);

	energy_output.close();
	bin_output.close();

	//escribimos los resultadosa un archivo de forma que sea facil de abrir con numpy e interpretar a simple vista
	std::ofstream results_file(results_filename, std::ios::app);
	results_file << temp << " " << p << " " << Cv << " " << gamma << " " << 1/inv_kT << std::endl;
	results_file.close();

	//borramos las variables auxiliares delta;
	delete[] deltax;
	delete[] deltay;
	delete[] deltaz;

}
