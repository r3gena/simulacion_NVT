#pragma once
#include "Sistema.h"
#include <ostream>
class Sistema_NVT :
    public Sistema
{
public:
    double MT;//masa del foco termico
    double T;//temperatura te�rica que tiene que dar, dato del sistema y la media tiene que salir similar a esto

    Sistema_NVT(std::string filename);

    ~Sistema_NVT();

    void stab_gear(int npasos, float dt, int save_rate, std::string energyfile, std::string velfile, std::string processname);
    void run_gear(int npasos, float dt, int save_rate, std::string resultsfile, std::string energyfile, std::string binfile, std::string processname);
    void exportar_auxiliares(std::string filename);//funcion que exporta el binario de los datos auxiliares
private:
    //derivadas de las posiciones
    double* dax;
    double* day;
    double* daz;
    double* d2ax;
    double* d2ay;
    double* d2az;
    //variable del foco y sus derivadas
    double s;//inicializar a 1 y las derivadas a 0
    double ds;
    double d2s;
    double d3s;
    double d4s;

    //definicion del algoritmo de gear de cuarto orden:
    void gear4(float dt, double deltax[], double deltay[], double deltaz[]);

};
