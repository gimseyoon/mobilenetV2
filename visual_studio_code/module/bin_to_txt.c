#include <stdio.h>
#include <stdlib.h>

int main() {
    FILE *binFile, *txtFile;
    const char *binFilePath = "C:/Digital_Circuit_Design_Project_ZIP/ETRI/mobilenetV2/parameter/data/bin/mobilenetv2_result_layer08.bin";
    const char *txtFilePath = "C:/Digital_Circuit_Design_Project_ZIP/ETRI/mobilenetV2/parameter/data/txt/mobilenetv2_result_layer08.txt";
    
    float value;
    
    // open binary file (read file)
    binFile = fopen(binFilePath, "rb");
    if (binFile == NULL) {
        printf("Error opening binary file.\n");
        return 1;
    }
    
    // open txt file (write file)
    txtFile = fopen(txtFilePath, "w");
    if (txtFile == NULL) {
        printf("Error creating text file.\n");
        fclose(binFile);
        return 1;
    }
    
    // Repeat to the end of the file
    while (fread(&value, sizeof(float), 1, binFile) == 1) {
        fprintf(txtFile, "%f\n", value);
    }
    
    printf("[ Transformation complete ]\n");
    
    fclose(binFile);
    fclose(txtFile);
    
    return 0;
}
