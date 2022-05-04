# include "tsunami.h"

int main(void)
{   

        char *meshName = "../data/PacificFine.txt";
        char *resultBaseName = "../output/tsunamiFine";
        
        tsunamiCompute(0.1,400,100,meshName,resultBaseName);
        tsunamiAnimate(0.1,400,100,meshName,resultBaseName);
         
        exit(0);     
}