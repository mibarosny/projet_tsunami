# include "tsunami.h"

int main(void)
{   

        char *meshName = "../data/PacificFine.txt";
        char *resultBaseName = "../output/tsunamiFine";
        
        tsunamiCompute(2.,10000,100,meshName,resultBaseName);
        tsunamiAnimate(2.,10000,100,meshName,resultBaseName);
         
        exit(0);     
}