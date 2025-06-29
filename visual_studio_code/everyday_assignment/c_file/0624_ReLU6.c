#include <stdio.h>

double ReLU6(double x);

int main()
{
    double x=0;
    
    scanf("%lf",&x);
    printf("x : %lf\n",x);
    
    
    x = ReLU6(x);
    printf("x after ReLU6 : %lf", x);

    return 0;
} 

 
double ReLU6(double x){
    
    // 0 < x 
    if(x < 0) {
        return 0;}

    //0 <= x <= 6
    else if( (0 <= x) && ( x <=6 ) )  {
        return x;}

    // x > 6 
    else {
        return 6;}
} 