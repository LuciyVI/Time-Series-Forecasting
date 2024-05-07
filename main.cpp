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
            
            double error = calculate_error(input(i),right_value(i));
            std::cout<<"Error : "<<error<<std::endl;
            weights(i+1)+=weights(i)*(delta_wieghts(error,input(i)));
            
            sum_error += error*error;
            
        }    
        
   }




void ploting_grph(double y , double x, Eigen::VectorXd input){

    Gnuplot gp;
    gp << "set xrange [-5:5]\nset yrange [-5:5]\n";
    gp << "plot '-' with lines title 'y = x^2'\n";
    gp.send1d(input);

}



int number_points = 20;


int main() {
    
    double a = -5.0;
    double b = 5.0;
    
    int len_window = 4;

    Eigen::VectorXd wights(len_window);
    Eigen::VectorXd right_vector(number_points);
    Eigen::VectorXd input(number_points);
    Eigen::VectorXd new_answer(number_points);
    Eigen::VectorXd X(1000);
    Eigen::VectorXd Y(1000);
    int count = 1;
    
    do{
        
        input(count-1)=a+(count-1)*(b-a)/number_points;
        right_vector(count-1)=linearActivation(input(count-1));
        input(input.size()-count)=b+count*(b-a)/number_points;
        new_answer(count-1)=linearActivation(input(input.size()-count));
        count++;
        

    }while(!input.all());
    
    std::cout<<"Right_Vector : "<<right_vector<<std::endl;
    std::cout<<"inputs :  "<<input<<std::endl;
    for(int i = 0; i<1000;i++){
        X(i)=a-0.5+i*(1+2*b-a)/100;
    }
    for(int i = 0; i<1000;i++){
        Y(i)=linearActivation(X(i));
    }
    double sum_error = 0.0;
    double *ptr  = &sum_error;

    wights<< 0,1,1,1;

    distribution_weights(wights,input,right_vector,len_window,number_points,*ptr);


    ploting_grph(4,5,Y);

    std::cout<<"Net:"<<net(wights,input,len_window)<<std::endl;

    return 0;

}