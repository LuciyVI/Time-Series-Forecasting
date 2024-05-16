#include </usr/include/termocolor/termcolor.hpp>
#include <boost/format.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <random>
#include <cmath>
#include <gnuplot-iostream.h>

#define LR 0.4
double epsil = 0.01;

double linearActivation(double x) {

    return 0.5*cos(0.5*x)-0.5;
}

double net(Eigen::VectorXd weights, Eigen::VectorXd inputs, int len_window){
 
    weights(0)=0.0;

    double net = 0.0;
    
        for(int i=1 ;i<len_window;++i){
            net += weights(i)*inputs(i-1);
        }
    
    return net;

}
   
double delta_wieghts(double error,double input){

    
    return LR*error*input;

}

double calculate_error(double input, double right_value)
{
    return right_value - input;
}

void distribution_weights(Eigen::VectorXd &weights, Eigen::VectorXd input,Eigen::VectorXd right_value,int len_window , int number_points,double &sum_error)
   {    
    


        for(int i = 0; i<weights.size()-1;++i)
        {
            //TODO передать изменённые веса в main
            double error = calculate_error(input(i),right_value(i));
            weights(i+1)+=weights(i)*(delta_wieghts(error,input(i)));
            std::cout<<"w : "<<weights<<std::endl;
            
            sum_error += error*error;
            
        }    
        
   }




void ploting_grph(double y , double x, std::vector<double> input){

    Gnuplot gp;
    gp << "set xrange [y:x]\nset yrange [y:x]\n";
    gp << "plot '-' with lines title 'func = 0.5*cos(0.5*x)-0.5 \n";
    gp.send1d(input);

}



int number_points = 20;


int main() {
    
    double a = -5.0;
    double b = 5.0;
    
    int len_window = 4;

    std::vector<double> wights(len_window);
    
    Eigen::VectorXd input(number_points);
    Eigen::VectorXd new_answer(number_points);
    Eigen::VectorXd X(1000);
    Eigen::VectorXd Y(1000);
    int count = 1;
    for(int i = 0;i<len_window;i++)
    {
        wights[i]=1.0;
    }
    

    std::vector<double>inputs_vector;
    Eigen::VectorXd inputs_eigen(number_points+1);
    std::vector<double>right_vector;
    

    for (int i = 0; i < number_points+1; ++i) {
        inputs_vector.push_back(a + i * (b - a) / number_points);
        inputs_eigen(i) = a + i * (b - a) / number_points;
        right_vector.push_back(linearActivation(inputs_eigen(i)));
    }  

    std::cout << "Using vector with push_back:" << std::endl;
        for (const auto &input : inputs_vector) {
            std::cout << input << " ";

        }
    std::cout << std::endl;

    std::cout << "Using Eigen::VectorXd:" << std::endl;
    std::cout << inputs_eigen << std::endl;

    std::cout << "Right vector with push_back:" << std::endl;
        for (const auto &right_value : right_vector) {
         
            std::cout<<"\n"<< right_value << " ";
            
        }


    for(int i = 0; i<1000;i++){
        X(i)=a-0.5+i*(1+2*b-a)/100;
    }
    for(int i = 0; i<1000;i++){
        Y(i)=linearActivation(X(i));
    }
    

    std::cout<<"W[0] "<<wights[0]<<"\t W[1]"<<wights[1]<<"\t W[2]"<<wights[2]<<"\t W[3]"<<wights[3]<<std::endl;
    double sum_error = 0.0;
    double *ptr  = &sum_error;


    std::cout<<"Sum_error : "<<sum_error<<std::endl; 

    ploting_grph(4,5,right_vector);
        
        std::cout<<"Sum_error : "<<sum_error<<std::endl;        
        // distribution_weights(wights,input,right_vector,len_window,number_points,*ptr);
        
        // std::cout<<"Net:"<<net(wights,input,len_window)<<std::endl;
    
    
    return 0;

}