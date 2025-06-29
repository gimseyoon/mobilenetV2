#include <stdio.h>
#include <stdlib.h>
#include "bitmap.h"

int main()
{
    FILE *fpBmp;
    FILE *fpBmp_new;
    FILE *fpBmp_blue;
    FILE *fpBmp_reverse;

    BITMAPFILEHEADER fileHeader;
    BITMAPINFOHEADER fileInfo;
    unsigned char *PixelData;
    unsigned char *PixelData_BLUE;
    unsigned char *PixelData_REVERSE;


    fpBmp = fopen("C:/seyoon/visual_studio_code/bitmap_file/lenna.bmp","rb");
    fpBmp_new = fopen("C:/seyoon/visual_studio_code/bitmap_file/lenna_new.bmp","wb");
    fpBmp_blue = fopen("C:/seyoon/visual_studio_code/bitmap_file/lenna_blue.bmp","wb");
    fpBmp_reverse = fopen("C:/seyoon/visual_studio_code/bitmap_file/lenna_reverse.bmp","wb");

    fread(&fileHeader, 14, 1, fpBmp);
    fread(&fileInfo, 40, 1, fpBmp);

    printf("fileInfo.biWidth : %d\n", fileInfo.biWidth);
    printf("fileInfo.biHeight : %d\n", fileInfo.biHeight);
    printf("fileInfo.biBitCount : %d\n", fileInfo.biBitCount);
    printf("fileInfo.biSizeImage : %d\n", fileInfo.biSizeImage);
    


    PixelData = (unsigned char *)malloc(fileInfo.biSizeImage);
    PixelData_BLUE = (unsigned char *)malloc(fileInfo.biSizeImage);
    PixelData_REVERSE = (unsigned char *)malloc(fileInfo.biSizeImage); 


    
    
    ///////////////////////////////////////////////////////
    // lenna_new //////////////////////////////////////////
    ///////////////////////////////////////////////////////
    // fread
    fread(PixelData, sizeof(unsigned char), fileInfo.biSizeImage, fpBmp);
    // print
    
    for (int i = 1535; i > 1435; i--) {
        printf("PixelData[%d] : %d\n", i,PixelData[i]);
    }
    printf("\n");
    

    //fwrite
    fwrite(&fileHeader, 14, 1, fpBmp_new);
    fwrite(&fileInfo, 40, 1, fpBmp_new);
    fwrite(PixelData, sizeof(unsigned char),fileInfo.biSizeImage, fpBmp_new );



    ///////////////////////////////////////////////////////
    // lenna_blue /////////////////////////////////////////
    ///////////////////////////////////////////////////////

    // compute
    for(int i = 0; i < fileInfo.biSizeImage; i++)
    {
        if( !(i%3 ==0) ) PixelData_BLUE[i]= 0;
        else PixelData_BLUE[i]=PixelData[i];
    }
    //print
    /*    
    for (int i = 0; i < 100; i++) {
        printf("PixelData_BLUE[%d] :  %d\n ", i,PixelData_BLUE[i]);
    }
    printf("\n");
    */
    // fwrite
    fwrite(&fileHeader, 14, 1, fpBmp_blue);
    fwrite(&fileInfo, 40, 1, fpBmp_blue);
    fwrite(PixelData_BLUE, sizeof(unsigned char),fileInfo.biSizeImage, fpBmp_blue );



    ///////////////////////////////////////////////////////
    // lenna_reverse //////////////////////////////////////
    ///////////////////////////////////////////////////////
    //compute
    for (int i = 0; i < fileInfo.biHeight ; i++) // 512
    {
        for(int j = 0; j< fileInfo.biWidth; j++) //512
        {
            PixelData_REVERSE[(512*i + j)*3 ] = PixelData[(512*i + 512-(j+1))*3 ];     //blue
            PixelData_REVERSE[(512*i + j)*3+1 ] = PixelData[(512*i + 512-(j+1))*3 +1];    //green
            PixelData_REVERSE[(512*i + j)*3+2 ] = PixelData[(512*i + 512-(j+1))*3 +2];      //red
        }
    }

    //print
    for (int i = 0; i < 100; i++) {
        printf("PixelData_REVERSE[%d] : %d\n", i,PixelData_REVERSE[i]);
    }
    printf("\n");

    //fwrite
    fwrite(&fileHeader, 14, 1, fpBmp_reverse);
    fwrite(&fileInfo, 40, 1, fpBmp_reverse);
    fwrite(PixelData_REVERSE, sizeof(unsigned char),fileInfo.biSizeImage, fpBmp_reverse );


    ///////////////////////////////////////////////////////
    // fclose /////////////////////////////////////////////
    ///////////////////////////////////////////////////////
    fclose(fpBmp);
    fclose(fpBmp_new);
    fclose(fpBmp_blue);
    fclose(fpBmp_reverse);
    
    return 0;
}
