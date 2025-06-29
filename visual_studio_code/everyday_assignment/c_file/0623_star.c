#include <stdio.h>

int main()
{
    int end=0, num = 0;


    while(end==0){
        
        printf("enter number : ");
        scanf("%d",&num);



        if(num<0)
        {
            printf("error!");
            end = 1;
        }



        else if(num==0)
        {
            end = 1;
        }



        else // num > 0
        {
            for(int i = 0; i<num; i++)
            {
                printf("*");
                
                if(i==(num-1)){
                    printf("\n");
                    num--;
                    i=-1;
                }
                
            } 
        }
    }
        
    
    
    return 0;


} 