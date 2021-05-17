#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <ctime>
#include <random>
#include <list>
#include "vecFunc.cpp"
#include "wav_generator.cpp"


class String{
private:
    double T = 102.6; //Tension
    double time = 0.000001; //в референсной статье -- это k -- шаг по времени
    double rho = 2500;
    double r = 0.44e-3; //radius of string
    double E = 8.6e9; //Young module
    double S = M_PI * pow(r,2); //cross area
    double c = pow(T/(rho*S), 0.5); //wave speed
    double Coef_1 = T*time*time/(rho*S);
    double I = M_PI*pow(r,4)/4;
    double L = 1;
    double Coef_2 = E*time*time*I/(rho*S);
	double h = c*time*2; //Размер клетки (расстояние между соседними юнитами)
    long int Number_of_particles = (int)(L/h);
    
    std::vector<std::vector<double> > Data;
    
    std::vector<float> sound;
    
    double coord_i = (int)(0.7*L/h); //Координата считывания звука


	void generate_String(){
        int n = Number_of_particles;
		for (int i = 0; i < Number_of_particles; i++){
		   // линейный профиль
           /*
            if (i < 0.7*Number_of_particles){
                Data[0][i] = i/(7.0*n); //Предыдущее состояние сетки
                Data[1][i] = i/(7.0*n); //Текущее состояние
                Data[2][i] = 0; //Следующее состояние сетки
                
            }else{
                Data[0][i] = 1.0/3.0-i/(3.0*n); //Предыдущее состояние сетки
                Data[1][i] = 1.0/3.0-i/(3.0*n); //Текущее состояние
                Data[2][i] = 0; //Следующее состояние сетки
                    
            }
            */
            /* Работающий профиль 4 степени*/
            Data[0][i] = 1e-13*i*i*(i - Number_of_particles)*(i - Number_of_particles); //Предыдущее состояние сетки
            Data[1][i] = 1e-13*i*i*(i - Number_of_particles)*(i - Number_of_particles); //Текущее состояние
            Data[2][i] = 0; //Следующее состояние сетки

             // let's try 2nd degree parabolic profile
             /*
            Data[0][i] = 1e-4*i*(i - Number_of_particles); //Предыдущее состояние сетки
            Data[1][i] = 1e-4*i*(i - Number_of_particles); //Текущее состояние
            Data[2][i] = 0; //Следующее состояние сетки
              */
		}
	}

    
    std::vector<double> Operator_Dxx(int k){
        std::vector<double> b(Number_of_particles);
        for (int i = 1; i < Number_of_particles - 1; i++)
                b[i] = - 2 * Data[k][i] + Data[k][i - 1] + Data[k][i + 1];
    
        b[0] = -2 * Data[k][0] + Data[k][1];
        b[Number_of_particles - 1] = -2 * Data[k][Number_of_particles - 1] + Data[k][Number_of_particles - 2];
        
        for (int i = 0; i < Number_of_particles; i++)
            b[i] = b[i] / (h*h);
        
        return b;
    }
    
    std::vector<double> Operator_Dxxxx(int k){ //хахахаха
        std::vector<double> b(Number_of_particles);
        for (int i = 2; i < Number_of_particles - 2; i++)
            b[i] = Data[k][i - 2] - 4 * Data[k][i - 1] + 6 * Data[k][i] - 4 * Data[k][i + 1] + Data[k][i + 2];
    
        b[0] = 5 * Data[k][0] - 4 * Data[k][1] + Data[k][2];
        b[1] = - 4 * Data[k][0] + 6 * Data[k][1] - 4 * Data[k][2] + Data[k][3];
        
        b[Number_of_particles - 2] = Data[k][Number_of_particles - 4] - 4 * Data[k][Number_of_particles - 3] + 6 * Data[k][Number_of_particles - 2] - 4 * Data[k][Number_of_particles - 1];
        b[Number_of_particles - 1] = Data[k][Number_of_particles - 3] - 4 * Data[k][Number_of_particles - 2] + 5 * Data[k][Number_of_particles - 1];
        
        for (int i = 0; i < Number_of_particles; i++)
            b[i] = b[i] / (h*h*h*h);
        
        return b;
    }
    
public:
    
    String(){
        Data = std::vector<std::vector<double> > (3, std::vector<double>(Number_of_particles));
        generate_String();
        std::cout << (int)(L/h) << std::endl;
        std::cout << h << std::endl;
    }
    
	void output(std::ofstream &f){
		f << Number_of_particles << "\n\n";
		for (int j = 0; j < Number_of_particles; j++){
			f << "10\t" << (j - Number_of_particles/2)*h << "\t" << "0" << "\t" << Data[1][j] << "\n";
		}
	}
    
    void move(){
        
        std::vector<double> b(Number_of_particles);
        b = Operator_Dxx(1);
        
        std::vector<double> a(Number_of_particles);
        a = Operator_Dxxxx(1);
        
        Data[2] = Coef_1 * b  + 2 * Data[1] - Data[0]- Coef_2 * a;
        Data[0] = Data[1];
        Data[1] = Data[2];
        
        sound.push_back(Data[1][coord_i]);
        
    }
    
    std::vector<float> Get_some_sound(){
        return sound;
    }
    
};


int main(){
	
    int Num_step = 5000; //Количество шагов рассчета
    std::ofstream fout;
    fout.open("Grid.xyz");
	String c = String();
    
    for (int i = 0; i < Num_step; i++){
        if (i % 10 == 0)
            c.output(fout);
        c.move();
    }
    //std::cout << M_PI << std::endl;
    fout.close();
    
    std::vector<float> data;
    data = c.Get_some_sound();
    writeWAV(data, "Attempt1_Non_Stiff");
    
    return 0;
}
