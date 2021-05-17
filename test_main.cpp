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
    double time = 1 / 44100.0; //в референсной статье -- это k -- шаг по времени
    double rho = 2500;
    double r = 0.44e-3; //radius of string
    double E = 8.6e9; //Young module
    double S = M_PI * pow(r,2); //cross area
    double c = pow(T/(rho*S), 0.5); //wave speed
    double Coef_1 = T*time*time/(rho*S);
    double I = M_PI*pow(r,4)/4;
    double L = 0.32;
    double Coef_2 = E*time*time*I/(rho*S);
	double h = 1.02 * std::sqrt((Coef_1 + std::sqrt(Coef_1 * Coef_1 + 16 * Coef_2))/2); //Размер клетки (расстояние между соседними юнитами)
    long int Number_of_particles = (int)(L/h);
    
    long int m1 = 4; //количество первых слагаемых (гаммы)
    long int m2 = 3; //количество вторых слагаемых (кси)
    
    
    std::vector<double> a1 = {1.00000000e-10, 2.15757554e+03, 2.30004601e+03, 4.76438766e+03};
    std::vector<double> b1 = {3.65939041e-04, 4.46237001e-04, 5.81732138e-05, 1.62899478e-03};
    
    std::vector<double> a2 = {2.45283773e+02, 2.42702880e+02, 2.35867275e+05};
    std::vector<double> b2 = {1.83439737e-07, 8.72746675e-08, 1.22966394e-06};
     
    
    std::vector<std::vector<double>> Data;
    
    std::vector<std::vector<std::vector<double>>> chi;
    std::vector<std::vector<std::vector<double>>> gamma;
    
    std::vector<float> sound;
    
    double coord_i = (int)(0.7*L/h); //Координата считывания звука

    
    double A = 0.0;
    double B = 0.0;
    

	void generate_String(){
        int n = Number_of_particles;
		for (int i = 0; i < Number_of_particles; i++){
            
            Data[0][i] = 0.05 / pow(Number_of_particles / 2, 4)*i*i*(i - Number_of_particles)*(i - Number_of_particles); //Предыдущее состояние сетки
            Data[1][i] = 0.05 / pow(Number_of_particles / 2, 4)*i*i*(i - Number_of_particles)*(i - Number_of_particles); //Текущее состояние
            Data[2][i] = 0; //Следующее состояние сетки
            
		}
	}

    
    std::vector<double> Operator_Dxx(std::vector<double> const &x){
        std::vector<double> b(Number_of_particles);
        for (int i = 1; i < Number_of_particles - 1; i++)
                b[i] = - 2 * x[i] + x[i - 1] + x[i + 1];
    
        b[0] = -2 * x[0] + x[1];
        b[Number_of_particles - 1] = -2 * x[Number_of_particles - 1] + x[Number_of_particles - 2];
        
        for (int i = 0; i < Number_of_particles; i++)
            b[i] = b[i] / (h*h);
        
        
        return b;
    }
    
    std::vector<double> Operator_Dxxxx(std::vector<double> const &x){ //хахахаха
        std::vector<double> b(Number_of_particles);
        for (int i = 2; i < Number_of_particles - 2; i++)
            b[i] = x[i - 2] - 4 * x[i - 1] + 6 * x[i] - 4 * x[i + 1] + x[i + 2];
    
        b[0] = 5 * x[0] - 4 * x[1] + x[2];
        b[1] = - 4 * x[0] + 6 * x[1] - 4 * x[2] + x[3];
        
        b[Number_of_particles - 2] = x[Number_of_particles - 4] - 4 * x[Number_of_particles - 3] + 6 * x[Number_of_particles - 2] - 4 * x[Number_of_particles - 1];
        b[Number_of_particles - 1] = x[Number_of_particles - 3] - 4 * x[Number_of_particles - 2] + 5 * x[Number_of_particles - 1];
        
        for (int i = 0; i < Number_of_particles; i++)
            b[i] = b[i] / (h*h*h*h);
    
        
        return b;
    }

    
    std::vector<double> Run(std::vector<double> d, double A, double B, int n) {
        
        std::vector<double> a (Number_of_particles, A);
        std::vector<double> b (Number_of_particles, B);
        std::vector<double> c (Number_of_particles, A);


        
        n--;
        c[0] /= b[0];
        d[0] /= b[0];

        for (int i = 1; i < n; i++) {
            c[i] /= b[i] - a[i]*c[i-1];
            d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
        }

        d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

        for (int i = n; i-- > 0;) {
            d[i] -= c[i]*d[i+1];
        }
        
        std::vector<double> Result (Number_of_particles);
        
        Result = d;
        
        return Result;
    }
    
    
    
public:
    
    void check(std::vector<double> x, int N, double A, double B){
        
        x = Run(x, A, B, N);
        
        for (int i  = 0; i < N; i++)
            std::cout << x[i] << " ";
    }
    
    
    String(){
        Data = std::vector<std::vector<double>> (3, std::vector<double>(Number_of_particles));
        
        chi = std::vector<std::vector<std::vector<double>>> (m1, std::vector<std::vector<double>>(3, std::vector<double>(Number_of_particles, 0)));
        gamma = std::vector<std::vector<std::vector<double>>> (m1, std::vector<std::vector<double>>(3, std::vector<double>(Number_of_particles, 0)));
        
        B += rho*S/(time * time);
        std::cout << B << '\n';
        for (int i = 0; i < m1; i++){
            B += b1[i] / (time * (2 + a1[i] * time));
        }
        
        
        for (int i = 0; i < m2; i++){
            B += 2 * b2[i]/(time * (2 + a2[i] * time) * h * h);
            A -= b2[i]/(time * (2 + a2[i] * time) * h * h);
        }
        
        generate_String();
        std::cout << A << " " << B << '\n';
    }
    
	void output(std::ofstream &f){
		f << Number_of_particles << "\n\n";
		for (int j = 0; j < Number_of_particles; j++){
			f << "10\t" << (j - Number_of_particles/2)*h << "\t" << "0" << "\t" << Data[1][j] << "\n";
		}
	}
    
    void move(){
        
        std::vector<double> super_vector(Number_of_particles, 0);
        
        super_vector = 2 * rho*S*Data[1]/(time * time) + T * Operator_Dxx(Data[1]) - E * I * Operator_Dxxxx(Data[1]) - rho * S *Data[0] / (time * time);
        //std::cout<< super_vector[0] << " ";
        for (int i = 0; i < m1; i++){
            super_vector += (1 - (2 - a1[i] * time) / (2 + a1[i] * time)) * b1[i] / time * gamma[i][1];
            super_vector += b1[i] / (time * (2 + a1[i] * time)) * Data[0];
        }
        //std::cout<< super_vector[0] << " ";

        for (int i = 0; i < m2; i++){
            super_vector -= (1 - (2 - a2[i] * time) / (2 + a2[i] * time)) * b2[i] / time * Operator_Dxx(chi[i][1]);
            super_vector -= b2[i] / (time * (2 + a2[i] * time)) * Operator_Dxx(Data[0]);
        }

        Data[2] = Run(super_vector, A, B, Number_of_particles);
        //std::cout << Data[2][0] << " ";
        for (int i = 0; i < m1; i++){
            gamma[i][2] = (Data[2] + (2 - a1[i] * time) * gamma[i][1] - Data[0]) / (2 + a1[i] * time);
            
            gamma[i][0] = gamma[i][1];
            gamma[i][1] = gamma[i][2];
        }
        //std::cout<<gamma[0][2][0] << '\n';
        
        for (int i = 0; i < m2; i++){
            chi[i][2] = (Data[2] + (2 - a2[i] * time) * chi[i][1] - Data[0]) / (2 + a2[i] * time);
            
            chi[i][0] = chi[i][1];
            chi[i][1] = chi[i][2];
        }

        Data[0] = Data[1];
        Data[1] = Data[2];

        sound.push_back(Data[1][coord_i]);
        
        
    }
    
    std::vector<float> Get_some_sound(){
	    float max_a = 0.0;
	    for(int i = 0; i < sound.size(); i++)
	        max_a = std::max(sound[i], max_a);
	    for(int i = 0; i < sound.size(); i++)
	        sound[i] = sound[i] / max_a;
        return sound;
    }

    
};


int main(){
	
    int Num_step = 50000; //Количество шагов рассчета
    
    std::ofstream fout;
    fout.open("Grid.xyz");
	String c = String();
    
    
    for (int i = 0; i < Num_step; i++){
        if (i % 10 == 0)
            c.output(fout);
        c.move();
    }

    fout.close();
    
    std::vector<float> data;
    data = c.Get_some_sound();
    writeWAV(data, "Attempt1");
    
    
    std::vector<double> a = {3,4,3};
    
    c.check(a,3,1,2);
    
    return 0;
}
