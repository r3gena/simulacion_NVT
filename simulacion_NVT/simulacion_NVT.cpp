#include "Sistema_NVT.h"
#include <iostream>
#include <string>

int Total_pasos=100000;
int num_procesos=1;
int pasos_proceso=Total_pasos/num_procesos;
int save_rate=10;

float dt=0.0001;

std::string energy_name;
std::string binary_name;
std::string process_name;
std::string results_name="../results.txt";

int main()
{
  Sistema_NVT sistema_nvt("../estabilizado.txt");//tiene pinta que esto esta llamando primero a la clase padre

  for(int i=0; i<num_procesos; i++)
  {
    energy_name="/home/pedro/Documentos/Simulacion/run/Energy/energy_";
		binary_name="/home/pedro/Documentos/Simulacion/run/Binary/bin_";
    process_name="Process ";

    energy_name.append(std::to_string(i+1));
		binary_name.append(std::to_string(i+1));
    process_name.append(std::to_string(i+1));

    energy_name.append(".txt");
		binary_name.append(".bin");

    sistema_nvt.run_gear(pasos_proceso, dt, save_rate, results_name, energy_name, binary_name, process_name);
  }
  return 0;
}
