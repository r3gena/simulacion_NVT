#include "potenciales.h"

//llamamos por pointer a los potenciales para poder modificarlos en el origen, por defecto los array se llaman por pointer asi que no hay que especificarlo
void lennard_jones(const int N, const float rc, const float rc2, const float L, const float V, double* E_p, double* dphi, double* d2phi, const double rx[], const double ry[], const double rz[], double ax[], double ay[], double az[])
{
	*E_p = 0;
	*dphi = 0;
	*d2phi = 0;
	
	double rijx, rijy, rijz, rij2, rij6, inv_rij6, inv_rij12, fmod, ep_corr, dphi_corr, d2phi_corr;

	//primero nos aseguramos de que las aceleraciones estan a 0
	for (int i = 0; i < N; i++)
	{
		ax[i] = 0;
		ay[i] = 0;
		az[i] = 0;
	}

	//ahora calculamos el potencial en unnidades reducidas: uij=4*(1/rij^12-1/rij^6)
	for (int i = 0; i < N - 1; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			//primero calculamos la distancia entre particulas para ver si estan dentro del radio de corte
			rijx = rx[i] - rx[j];
			rijy = ry[i] - ry[j];
			rijz = rz[i] - rz[j];

			//corregimos las distancias para que la particula interactue con la copia mas cercana
			rijx = rijx - L * std::round(rijx / L);
			rijy = rijy - L * std::round(rijy / L);
			rijz = rijz - L * std::round(rijz / L);

			rij2 = rijx * rijx + rijy * rijy + rijz * rijz;

			//ahora si las particulas estan en rango se calcula la interaccion
			if (rij2 <= rc2)
			{
				rij6 = rij2 * rij2 * rij2;
				inv_rij6 = 1 / rij6;
				inv_rij12 = inv_rij6 * inv_rij6;
				fmod = (2 * inv_rij12 - inv_rij6) /rij2;

				ax[i] = ax[i] + fmod * rijx;
				ay[i] = ay[i] + fmod * rijy;
				az[i] = az[i] + fmod * rijz;
				ax[j] = ax[j] - fmod * rijx;
				ay[j] = ay[j] - fmod * rijy;
				az[j] = az[j] - fmod * rijz;

				//empezamos a calcular las derivadas del potencial con V
				*E_p = *E_p + (inv_rij12 - inv_rij6);
				*dphi = *dphi + (-2 * inv_rij12 + inv_rij6);
				*d2phi = *d2phi + (26 * inv_rij12 - 7 * inv_rij6);
			}
		}
	}

	//reescalamos las aceleraciones por un factor 24, lo hacemos aqui para reducir el num de operaciones
	for (int i = 0; i < N; i++)
	{
		ax[i] = ax[i] * 24;
		ay[i] = ay[i] * 24;
		az[i] = az[i] * 24;
	}
	//TODO: correcciones a la energia
	double f = 3.141592653589793238462 * N * N / (V * rc * rc * rc);
	ep_corr = 8 * f * (1. / (3 * pow(rc, 6)) - 1) / 3.;
	dphi_corr = 16 * f * (-2. / (3 * pow(rc, 6)) + 1);
	d2phi_corr = 16 * f * (26. / (3 * pow(rc, 6)) - 7);

	//añadimos las correcciones  y las operaciones finales para el calculo de las derivadas de phi
	*E_p = *E_p * 4 + ep_corr;
	*dphi = *dphi * 24 + dphi_corr;
	*d2phi = *d2phi * 24 + d2phi_corr;
	*d2phi = (*d2phi - *dphi * 2) / (9 * V * V);
	*dphi = *dphi / (3 * V);
}