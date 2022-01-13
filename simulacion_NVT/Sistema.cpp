#include "Sistema.h"

Sistema::Sistema()//constructor vacio
{
  
}

Sistema::Sistema(int num_k, int N, float rho, float E_t)
{
	//inicializamos las constantes del sistema
  this->N = N;
	this->rho = rho;
  V = N / rho;
  L = std::cbrt(V);
  double a = L / num_k;
  this->E_t = E_t;
  rc = L / 2;
  rc2 = rc * rc;
  //inicializamos la energia cinetica y potencial a 0 de forma provisional
  E_p = 0;
  E_cin = 0;
  dphi = 0;
  d2phi = 0;
	//inicializamos los pointers de r, v y a a arrays de tama�o N
	rx = new double[N];
	ry = new double[N];
	rz = new double[N];

	vx = new double[N];
	vy = new double[N];
	vz = new double[N];

	ax = new double[N];
	ay = new double[N];
	az = new double[N];

  colocar_red(num_k, a);//colocamos ya la situacion inicial puesto que si se quiere leer una externa ya hay un constructor para eso

}

Sistema::Sistema(std::string ascii_filename)
{
    std::string bin_filename;

    //primero leemos el archivo para inicializar las variables explicitadas
    std::ifstream ascii_in(ascii_filename, std::ios::in);
    ascii_in >> this->N >> this->rho >> this->L >> this->rc;
    ascii_in >> this->E_t >> this->E_p >> this->E_cin;
    ascii_in >> bin_filename;
    ascii_in.close();

    //a partir de estas variables inicializamos el resto
    this->rc2 = rc * rc;
    this->V = L * L * L;

    //inicializamos los pointers de r, v y a a arrays de tama�o N
    rx = new double[N];
    ry = new double[N];
    rz = new double[N];

    vx = new double[N];
    vy = new double[N];
    vz = new double[N];

    ax = new double[N];
    ay = new double[N];
    az = new double[N];

    dphi = 0;
    d2phi = 0;

    std::ifstream in_bin(bin_filename, std::ios::binary);
    in_bin.read((char*)rx, N * sizeof(double));
    in_bin.read((char*)ry, N * sizeof(double));
    in_bin.read((char*)rz, N * sizeof(double));
    in_bin.read((char*)vx, N * sizeof(double));
    in_bin.read((char*)vy, N * sizeof(double));
    in_bin.read((char*)vz, N * sizeof(double));
    in_bin.read((char*)ax, N * sizeof(double));
    in_bin.read((char*)ay, N * sizeof(double));
    in_bin.read((char*)az, N * sizeof(double));
    in_bin.close();
}

Sistema::~Sistema()
{
	delete[] rx;
	delete[] ry;
	delete[] rz;

	delete[] vx;
	delete[] vy;
	delete[] vz;

	delete[] ax;
	delete[] ay;
	delete[] az;
}

void Sistema::colocar_red(int num_k, double a)
{
    int part_asignadas = 0;

    for (int i = 0; i < num_k; i++)
    {
        for (int j = 0; j < num_k; j++)
        {
            for (int k = 0; k < num_k; k++)
            {
                //colocamos las pariculas de cada celda primitiva
                //primera particula
                rx[0 + part_asignadas] = (0 + 2 * i) * a / 2;
                ry[0 + part_asignadas] = (0 + 2 * j) * a / 2;
                rz[0 + part_asignadas] = (0 + 2 * k) * a / 2;
                //segunda particula
                rx[1 + part_asignadas] = (1 + 2 * i) * a / 2;
                ry[1 + part_asignadas] = (1 + 2 * j) * a / 2;
                rz[1 + part_asignadas] = (0 + 2 * k) * a / 2;
                //tercera particula
                rx[2 + part_asignadas] = (1 + 2 * i) * a / 2;
                ry[2 + part_asignadas] = (0 + 2 * j) * a / 2;
                rz[2 + part_asignadas] = (1 + 2 * k) * a / 2;
                //cuarta particula
                rx[3 + part_asignadas] = (0 + 2 * i) * a / 2;
                ry[3 + part_asignadas] = (1 + 2 * j) * a / 2;
                rz[3 + part_asignadas] = (1 + 2 * k) * a / 2;

                part_asignadas = part_asignadas + 4;
            }
        }
    }
    //creo el generador de num aleatorios
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-1, 1);

    //descolocamos las particulas aleatoriamente
    for (int i = 0; i < N; i++)
    {
        rx[i] = rx[i] + distribution(generator) * L / 100;
        ry[i] = ry[i] + distribution(generator) * L / 100;
        rz[i] = rz[i] + distribution(generator) * L / 100;

        //aprovecho este bucle para asignar las velocidades de forma aleatoria y sumar la energia cinetica
        vx[i] = distribution(generator);
        vy[i] = distribution(generator);
        vz[i] = distribution(generator);
    }

    //calculo del potencial
    lennard_jones(N, rc, rc2, L, V, &E_p, &dphi, &d2phi, rx, ry, rz, ax, ay, az);

    //Asignacion de energia cinetica
    double E_cin_obj = E_t - E_p;
    double px = 0, py = 0, pz = 0;
    E_cin = 0;
    for (int i = 0; i < N; i++)
    {
        //calculamos el momento total en cada direccion
        px += vx[i];
        py += vy[i];
        pz += vz[i];
    }
    //calculamos el momento por particula
    px = px / double(N);
    py = py / double(N);
    pz = pz / double(N);

    //ajustamos las velocidades para que el momento total del sistema sea 0 y calculamos la energia cinetica
    for (int i = 0; i < N; i++)
    {
        vx[i] -= px;
        vy[i] -= py;
        vz[i] -= pz;
        E_cin += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
    }
    px = 0; py = 0; pz = 0;
    for (int i = 0; i < N; i++)
    {
        //calculamos el momento total en cada direccion
        px += vx[i];
        py += vy[i];
        pz += vz[i];
    }
    std::cout << px << " " << py << " " << pz << std::endl;

    E_cin = E_cin * 0.5;
    double factor = std::sqrt(E_cin_obj / E_cin);
    //reescalamos las velocidades para que la energia cinetica total sea la que queremos y recalculamos la energia cinetica
    E_cin = 0;
    for (int i = 0; i < N; i++)
    {
        vx[i] = vx[i] * factor;
        vy[i] = vy[i] * factor;
        vz[i] = vz[i] * factor;
        E_cin += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
    }
    E_cin = 0.5 * E_cin;
}
void Sistema::imprimir_energias()
{
    std::cout << "Energias del sistema:---------------------" << std::endl;
    std::cout << "Energia total deseada:       " << E_t << std::endl;
    std::cout << "Energia potencial obtenida: " << E_p << std::endl;
    std::cout << "Energia cinetica obtenida:    " << E_cin << std::endl;
    std::cout << "Energia total obtenida:      " << E_p + E_cin << std::endl;
}

void Sistema::exportar_datos(std::string ascii_filename, std::string bin_filename)
{
    //escribimos el archivo ascii
    std::ofstream ascii_file(ascii_filename, std::ios::out);
    ascii_file << N << " " << rho << " " << L << " " << rc << std::endl;
    ascii_file << std::fixed << std::setw(10) << std::setprecision(6) << E_t << " " << E_p << " " << E_cin << std::endl;
    ascii_file << bin_filename;

    ascii_file.close();

}

void Sistema::exportar_datos(std::string ascii_filename, std::string bin_filename, float dt, int n_pasos)
{
    //escribimos el archivo ascii
    std::ofstream ascii_file(ascii_filename, std::ios::out);
    ascii_file << N << " " << rho << " " << L << " " << rc << n_pasos << dt << std::endl;
    ascii_file << std::fixed << std::setw(10) << std::setprecision(6) << E_t << " " << E_p << " " << E_cin << std::endl;
    ascii_file << bin_filename;

    ascii_file.close();

}

void Sistema::exportar_configuracion(std::string bin_filename)
{
    //escribimos el archivo binario
    int tamano_arrays = N * sizeof(double);
    std::ofstream bin_file(bin_filename, std::ios::binary);
    bin_file.write((char*)rx, tamano_arrays);
    bin_file.write((char*)ry, tamano_arrays);
    bin_file.write((char*)rz, tamano_arrays);
    bin_file.write((char*)vx, tamano_arrays);
    bin_file.write((char*)vy, tamano_arrays);
    bin_file.write((char*)vz, tamano_arrays);
    bin_file.write((char*)ax, tamano_arrays);
    bin_file.write((char*)ay, tamano_arrays);
    bin_file.write((char*)az, tamano_arrays);

    bin_file.close();
}

void Sistema::stab_verlet(int n_pasos, float dt, int save_rate, std::string energy_filename, std::string vel_filename, std::string process_name)
{
    std::ofstream energy_output(energy_filename, std::ios::out);
    std::ofstream velocity_output(vel_filename, std::ios::binary);
    progress_bar bar(process_name, n_pasos);
    int array_size = N * sizeof(double);

    for (int j = 0; j < n_pasos; j++)
    {
        velocity_verlet(dt);

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
}

void Sistema::run_verlet(int n_pasos, float dt, int save_rate, std::string results_file, std::string energy_filename, std::string bin_filename, std::string process_name)
{
    std::ofstream energy_output(energy_filename, std::ios::out);
    std::ofstream bin_output(bin_filename, std::ios::binary);
    progress_bar bar(process_name, n_pasos);
    int array_size = N * sizeof(double);
    double mean_Ec = E_cin;
    double mean_inv_Ec = 1 / E_cin;
    double mean_dphi = dphi;
    double mean_d2phi = d2phi;
    double mean_dphi_invEc = dphi / E_cin;
    double mean_dphi2_invEc = dphi*dphi / E_cin;

    for (int j = 0; j < n_pasos; j++)
    {
        velocity_verlet(dt);
        //calculamos las energias medias y la media de las derivadas
        mean_Ec += E_cin;
        mean_inv_Ec += 1 / E_cin;
        mean_dphi += dphi;
        mean_d2phi += d2phi;
        mean_dphi_invEc += dphi / E_cin;
        mean_dphi2_invEc += dphi*dphi / E_cin;

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
    mean_Ec = mean_Ec / n_pasos;
    mean_inv_Ec = mean_inv_Ec / n_pasos;
    mean_dphi = mean_dphi / n_pasos;
    mean_d2phi = mean_d2phi / n_pasos;
    mean_dphi_invEc = mean_dphi_invEc / n_pasos;
    mean_dphi2_invEc = mean_dphi2_invEc / n_pasos;

    //por ultimo usamos esto para calcular los parametros del sistema en unidades reducidas
    double T = 2 * mean_Ec / (3 * N - 3);
    double p = N * T / V - mean_dphi;
    double Cv = 1 / (1 + (2 / (3 * N - 3) - 1) * mean_Ec * mean_inv_Ec);
    double alpha_E = 1 / (V * (1 - 2 / (3 * N - 3)) * mean_Ec * mean_dphi_invEc - mean_dphi);
    double gamma = N / Cv + V * (0.5 * (3 * N - 3) - 1) * (mean_dphi * mean_inv_Ec - mean_dphi_invEc);
    double inv_ks = (N * T / V) * (1 + 2 * gamma - N / Cv) + V * mean_d2phi - V * ((3 * N - 3) * 0.5 - 1) * (mean_dphi2_invEc - 2 * mean_dphi * mean_dphi_invEc + mean_dphi * mean_dphi * mean_inv_Ec);

    //A partir de estos resultados calculamos el resto:
    double inv_kt = inv_ks - T * Cv * gamma * gamma / V;
    double Cp = Cv * inv_ks / inv_kt;
    double alpha_p = Cv * gamma / (inv_kt * V);
    double alpha_s = -1 / (gamma * T);

    //calcular alpha_E por el segundo metodo y comprobar que dan cosas similares

    energy_output.close();
    bin_output.close();


    std::ofstream results_output(results_file, std::ios::out | std::ios::app);
    //escribimos los resultadosa un archivo de forma que sea facil de abrir con numpy e interpretar a simple vista
    results_output << T << " " << p << " " << Cv << " " << Cp << " " << alpha_E << " " << alpha_p << " " << alpha_s << " " << gamma << " " << 1 / inv_ks << " " << 1 / inv_kt << std::endl;
    results_output.close();
}

void Sistema::velocity_verlet(float dt)
{
    E_cin = 0;//ponemos la energia cinetica a 0 para volver a calcularla con las aceleraciones

        //actualizamos las posiciones en funcion de las velocidades y aceleraciones actuales
    for (int i = 0; i < N; i++)
    {
        rx[i] += vx[i] * dt + 0.5 * ax[i] * dt * dt;
        ry[i] += vy[i] * dt + 0.5 * ay[i] * dt * dt;
        rz[i] += vz[i] * dt + 0.5 * az[i] * dt * dt;

        //actualizamos la velocidad con la aceleracion nueva
        vx[i] += 0.5 * ax[i] * dt;
        vy[i] += 0.5 * ay[i] * dt;
        vz[i] += 0.5 * az[i] * dt;
    }

    //actualizamos energia potencial y aceleraciones llamando al potencial para las nuevas posiciones
    lennard_jones(N, rc, rc2, L, V, &E_p, &dphi, &d2phi, rx, ry, rz, ax, ay, az);

    //actualizamos las velocidades con las nuevas aceleraciones y calculamos energia cinetica
    for (int i = 0; i < N; i++)
    {
        vx[i] += 0.5 * ax[i] * dt;
        vy[i] += 0.5 * ay[i] * dt;
        vz[i] += 0.5 * az[i] * dt;

        E_cin += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
    }
    E_cin = 0.5 * E_cin;
}

void Sistema::corregir_energia()
{
    double factor = std::sqrt((E_t-E_p) / E_cin);
    //reescalamos las velocidades para que la energia cinetica total sea la que queremos y recalculamos la energia cinetica
    E_cin = 0;
    for (int i = 0; i < N; i++)
    {
        vx[i] = vx[i] * factor;
        vy[i] = vy[i] * factor;
        vz[i] = vz[i] * factor;
        E_cin += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
    }
    E_cin = 0.5 * E_cin;
}
