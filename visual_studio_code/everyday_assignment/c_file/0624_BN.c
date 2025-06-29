#include <stdio.h>
#include <math.h>
double BN(double x);    /* Batch Normalization */

int main()
{
    double x=0;
    
    scanf("%lf",&x);
    printf("x : %lf\n",x);
    
    
    x = BN(x);
    printf("x after BN : %lf", x);

    return 0;
} 

 
double BN(double x){    /* Batch Normalization */

    
    double r=1, mean=2, var=0.5, b=2, e=0.0005;
    double y;   //result of BN

    y = r * ( (x - mean) / sqrt (var + e) ) + b;  

    return y;
} 