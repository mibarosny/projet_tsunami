# include "tsunami.h"

typedef struct {
    int nElem;
    int nNode;
    int nEdge;
    int size;
    double *nodes;
    int *elems;
    int *edges;
    double *E;
    double *U;
    double *V;
    double *FE;
    double *FU;
    double *FV;
} femTsunamiProblem;

void tsunamiComputeInitialize(femTsunamiProblem * problem) {
    double * nodes = problem->nodes; 
    int * elems = problem->elems; 
    for (int el = 0; el < problem->nElem; el++)
        for (int i = 0; i < 3; i++) {
            problem->E[el*3+i] = tsunamiInitialConditionOkada(nodes[3*elems[el*3+i]], nodes[3*elems[el*3+i]+1]);
            problem->U[el*3+i] = 0.;
            problem->V[el*3+i] = 0.;
        }
}

void femShallowTriangleMap(int index, int map[3]) {
    for (int i = 0; i < 3; i++) 
        map[i] = index*3 + i;
}

int ff = 0;
void femShallowEdgeMap(femTsunamiProblem *myProblem, int index, int map[2][2])
{
    int i,j,k;   
    for (j=0; j < 2; ++j) {
        int node = myProblem->edges[index*4 + j];
        for (k=0; k < 2; k++) {
            int elem = myProblem->edges[index*4 + k + 2];
            map[k][j] = (myProblem->nElem)*3;
            if (elem >= 0) {
                for (i=0; i < 3; i++) {
                    if (myProblem->elems[elem*3 + i] == node) {
                        map[k][j] = elem*3 + i;  }}}}}
}

void femDiscretePhi(int n, double xsi, double eta, double * phi) {
    if (n == 2) {
        phi[0] = 1 - xsi - eta;  
        phi[1] = xsi;
        phi[2] = eta;
    } else {
        phi[0] = (1.0 - xsi)/2.0;
        phi[1] = (1.0 + xsi)/2.0;  
    }
}

double interpolateClassic(double *phi, double *U, int *map, int n)
{
    double u = 0.0; int i;
    for (i=0; i <n; i++)
        u += phi[i]*U[map[i]];
    return u;
}

double interpolate(double *phi, double *U, int *map, int e, int p, int n)
{
    double u = 0.0; int i;
    for (i=0; i < n; i++)
        u += phi[i]*U[e*map[i]+p];
    return u;
}

void femShallowAddIntegralsElements(femTsunamiProblem * myProblem) {
    double *BE = myProblem->FE;
    double *BU = myProblem->FU;
    double *BV = myProblem->FV;
    double *E = myProblem->E;
    double *U = myProblem->U;
    double *V = myProblem->V;    
       
    double  xLoc[3],yLoc[3],phi[3],dphidx[3],dphidy[3];
    double  xsi,eta,weight,jac;
    double  y,e,u,v,x,h;
    int     i,j,k,elem,mapElem[3];

    int n = 3;
    for (elem=0; elem < myProblem->nElem; elem++) {
        femShallowTriangleMap(elem, mapElem);
        int *mapCoord = &(myProblem->elems[elem*n]);
        for (j=0; j < n; ++j) {
        	  xLoc[j] = myProblem->nodes[3*mapCoord[j]];
        	  yLoc[j] = myProblem->nodes[3*mapCoord[j]+1]; }
        jac = (xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]);
        dphidx[0] = (yLoc[1] - yLoc[2])/jac;
        dphidx[1] = (yLoc[2] - yLoc[0])/jac;
        dphidx[2] = (yLoc[0] - yLoc[1])/jac;
        dphidy[0] = (xLoc[2] - xLoc[1])/jac;
        dphidy[1] = (xLoc[0] - xLoc[2])/jac;
        dphidy[2] = (xLoc[1] - xLoc[0])/jac;        
        for (k=0; k < 3; k++) {
            xsi = gaussTriangleXsi[k];
            eta = gaussTriangleEta[k];
            weight = gaussTriangleWeight[k];     
            femDiscretePhi(2,xsi,eta,phi);               
            x = interpolate(phi, myProblem->nodes, mapCoord, 3, 0, n); 
            y = interpolate(phi, myProblem->nodes, mapCoord, 3, 1, n); 
            h = interpolate(phi, myProblem->nodes, mapCoord, 3, 2, n); 
            e = interpolate(phi,E,mapElem,1,0,n);
            u = interpolate(phi,U,mapElem,1,0,n);
            v = interpolate(phi,V,mapElem,1,0,n);
            //if (elem == 0 && k == 0)
            //printf("%f %f %f %f %f %f \n", y, x, h, e, u, v);
            
            double z3d = (4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
            double lat = asin(z3d/R)*180/PI;
            double coriolis = 2.*Omega*sin(lat);     
            for (i=0; i < n; i++) {
                BE[mapElem[i]] += ( (dphidx[i]*h*u + dphidy[i]*h*v)*((4.*R*R+x*x+y*y)/(4.*R*R))
                                + phi[i]*h*((x*u+y*v)/(R*R)) )*jac*weight;
                BU[mapElem[i]] += ( phi[i]*(coriolis*v - Gamma*u) + dphidx[i]*g*e*((4.*R*R+x*x+y*y)/(4.*R*R))
                                + phi[i]*(g*x*e)/(2*R*R))*jac*weight;
                BV[mapElem[i]] += ( phi[i]*(-coriolis*u - Gamma*v) + dphidy[i]*g*e*((4.*R*R+x*x+y*y)/(4.*R*R))
                                + phi[i]*(g*y*e)/(2*R*R))*jac*weight;
            }
        }
    }        
}


void femShallowAddIntegralsEdges(femTsunamiProblem * myProblem) {
    double *BE = myProblem->FE;
    double *BU = myProblem->FU;
    double *BV = myProblem->FV;
    double *E = myProblem->E;
    double *U = myProblem->U;
    double *V = myProblem->V;
    
    double  xEdge[2],yEdge[2],phiEdge[2];
    double  xsi,weight,jac;
    double  eL,eR,uL,uR,vL,vR,unL,unR;
    double  qe,qu,qv,x,y,h,factor;
    int     i,j,k,edge,mapEdge[2][2];

    for (edge=0; edge < myProblem->nEdge; edge++) {
        femShallowEdgeMap(myProblem,edge,mapEdge);
        for (j=0; j < 2; ++j) {
        	  int node = myProblem->edges[edge*4 + j];
        	  xEdge[j] = myProblem->nodes[node*3];
        	  yEdge[j] = myProblem->nodes[node*3 + 1]; }
        int *mapCoord = &(myProblem->edges[edge*4]);
        int boundary = (mapEdge[1][0] == myProblem->size-1);
        
        double dxdxsi = (xEdge[1] - xEdge[0]);
        double dydxsi = (yEdge[1] - yEdge[0]);
        double norm = sqrt(dxdxsi*dxdxsi + dydxsi*dydxsi);
        double nx =  dydxsi/norm;
        double ny = -dxdxsi/norm;
        jac = norm / 2.0;
        for (k=0; k < 2; k++) {
            xsi = gaussEdgeXsi[k];
            weight = gaussEdgeWeight[k];     
            femDiscretePhi(1,xsi,xsi,phiEdge);           
            
            eL = interpolate(phiEdge,E,mapEdge[0],1,0,2);
            eR = boundary ? eL : interpolate(phiEdge,E,mapEdge[1],1,0,2);
            uL = interpolate(phiEdge,U,mapEdge[0],1,0,2);
            uR = interpolate(phiEdge,U,mapEdge[1],1,0,2);
            vL = interpolate(phiEdge,V,mapEdge[0],1,0,2);
            vR = interpolate(phiEdge,V,mapEdge[1],1,0,2);
            y = interpolate(phiEdge,myProblem->nodes,mapCoord,3,1,2);
            x = interpolate(phiEdge,myProblem->nodes,mapCoord,3,0,2);
            h = interpolate(phiEdge,myProblem->nodes,mapCoord,3,2,2);
            unL = uL*nx+ vL*ny;
            unR = boundary ? -unL : uR*nx + vR*ny;
            factor = (4.*R*R+x*x+y*y)/(4.*R*R);

            qe =  0.5*h*   ( (unL+unR) + sqrt(g/h)*( eL-eR ) )*factor;
            qu =  0.5*g*nx*( ( eL+eR ) + sqrt(h/g)*(unL-unR) )*factor;
            qv =  0.5*g*ny*( ( eL+eR ) + sqrt(h/g)*(unL-unR) )*factor;        
            for (i=0; i < 2; i++) {
                BE[mapEdge[0][i]] -= qe*phiEdge[i]*jac*weight; 
                BU[mapEdge[0][i]] -= qu*phiEdge[i]*jac*weight; 
                BV[mapEdge[0][i]] -= qv*phiEdge[i]*jac*weight; 
                BE[mapEdge[1][i]] += qe*phiEdge[i]*jac*weight;
                BU[mapEdge[1][i]] += qu*phiEdge[i]*jac*weight;
                BV[mapEdge[1][i]] += qv*phiEdge[i]*jac*weight; }}}

}

void femShallowMultiplyInverseMatrix(femTsunamiProblem * myProblem) {
    double *BE = myProblem->FE;
    double *BU = myProblem->FU;
    double *BV = myProblem->FV;
         
    double invA[3][3] = {{18.0,-6.0,-6.0},{-6.0,18.0,-6.0},{-6.0,-6.0,18.0}};
    double BEloc[3],BUloc[3],BVloc[3];
 
    double xLoc[3],yLoc[3],jac;
    int    elem,i,j,mapElem[3];
    
    for (elem=0; elem < myProblem->nElem; elem++) {
        femShallowTriangleMap(elem, mapElem);
        int *mapCoord = &(myProblem->elems[elem*3]);
        for (j=0; j < 3; ++j) {
        	  xLoc[j] = myProblem->nodes[3*mapCoord[j]];
        	  yLoc[j] = myProblem->nodes[3*mapCoord[j]+1]; }
        jac = (xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]);
        for (i=0; i < 3; i++) {
            BEloc[i] = BE[mapElem[i]];
            BUloc[i] = BU[mapElem[i]];
            BVloc[i] = BV[mapElem[i]];
            BE[mapElem[i]] = 0.0; 
            BU[mapElem[i]] = 0.0; 
            BV[mapElem[i]] = 0.0; }
        for (i=0; i < 3; i++) { 
            for (j=0; j < 3; j++) {
                BE[mapElem[i]] += invA[i][j] * BEloc[j] / jac; 
                BU[mapElem[i]] += invA[i][j] * BUloc[j] / jac; 
                BV[mapElem[i]] += invA[i][j] * BVloc[j] / jac; }}}
}

void tsunamiUpdateRungeKutta(femTsunamiProblem * myProblem, double dt) {
    int i, j, size = myProblem->size;
    double  *U = myProblem->U,
            *V = myProblem->V,
            *E = myProblem->E;
    double *Uold, *Eold, *Vold, *Upredictor, *Epredictor, *Vpredictor, *KE, *KU, *KV;
    Uold        = malloc(sizeof(double)*size*3);
    Eold        = Uold + size;
    Vold        = Eold + size;
    Upredictor  = malloc(sizeof(double)*size*3);
    Epredictor  = Upredictor + size;
    Vpredictor  = Epredictor + size;
    KE          = malloc(sizeof(double)*size*3);
    KU = KE + size;
    KV = KU + size;
    for (i=0; i < size; i++) {
        KU[i] = 0.0;
        KE[i] = 0.0;
        KV[i] = 0.0;
        Uold[i] = U[i];  
        Vold[i] = V[i];  
        Eold[i] = E[i];  
    }

    const int    nStage   = 2;
    const double beta[2]  = {0.0,   1.0  }; 
    const double gamma[2] = {1./2., 1./2.};
    
    for (j = 0; j < nStage; j++) {
        for (i=0; i < size; i++) {
                U[i] = Uold[i] + dt * beta[j] * KU[i];
                V[i] = Uold[i] + dt * beta[j] * KU[i];
                E[i] = Uold[i] + dt * beta[j] * KU[i];
            }
            
            for (i=0; i < size; i++) {
        	    U[i] += dt * gamma[j] * KU[i]; 
            }
        }
    
    free(Uold);
    free(KE);
}


void tsunamiCompute(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName)
{ 
    femTsunamiProblem *problem = (femTsunamiProblem *) malloc(sizeof(femTsunamiProblem));
    int i, j, trash, *elems, *edges;
    double dtrash, *nodes;     

    /** lecture du fichier */
    
    FILE* file = fopen(meshFileName,"r");
    fscanf(file, "Number of nodes %d \n",&(problem->nNode));   
    /* FOR RAPPORT En fait, l’appel à malloc est coûteux en temps et il vaut mieux déclarer une fois un tableau de
    nNode*3, qui va ensuite être virtuellement découpé en plusieurs, que faire l’inverse. Ici, on va
    avoir besoin de nNOde éléments pour les vecteurs X, Y et bath, on préfère
    donc faire qu’un seul appel à malloc.
    
    Bla bla de cache aussi  */
    nodes = malloc(sizeof(double)*problem->nNode*3);
    for (i = 0; i < problem->nNode; i++) 
        fscanf(file,"%d : %le %le %le\n",&trash,&nodes[i*3],&nodes[i*3+1],&nodes[i*3+2]); 
    problem->nodes = nodes;
        
    fscanf(file, "Number of triangles %d \n",&(problem->nElem)); 
    elems = malloc(sizeof(int)*3*(problem->nElem));
    problem->size = 3*problem->nElem + 1;
    for (i = 0; i < problem->nElem; i++) 
        fscanf(file,"%d : %d %d %d \n",&trash,&elems[i*3],&elems[i*3+1],&elems[i*3+2]);   
    problem->elems = elems;

    fscanf(file, "Number of edges %d \n",&(problem->nEdge)); 
    edges = malloc(sizeof(int)*4*(problem->nEdge));
    for (i = 0; i < problem->nEdge; i++) 
        fscanf(file,"%d : %d %d : %d %d \n",&trash,&edges[i*4],&edges[i*4+1],&edges[i*4+2],&edges[i*4+3]); 
    problem->edges = edges;  

    fclose(file); 

    problem->E = malloc(sizeof(double)*(problem->size)*3);    // dégré de liberté
    problem->U = problem->E + problem->size;
    problem->V = problem->U + problem->size;   
    double  *E = problem->E, 
            *U = problem->U, 
            *V = problem->V;  

    problem->FE = malloc(sizeof(double)*(problem->size)*3);    // dégré de liberté
    problem->FU = problem->FE + problem->size;
    problem->FV = problem->FU + problem->size; 

    tsunamiComputeInitialize(problem);
    tsunamiWriteFile(baseResultName,0,U,V,E,problem->nElem,3);

    FILE * uf = fopen("uf.csv", "w+");
    FILE * ef = fopen("ef.csv", "w+");
    FILE * vf = fopen("vf.csv", "w+");


    for (int it = 0; it < nmax; it++) {
        for (i=0; i < problem->size; i++) {
            problem->FE[i] = 0.0;
            problem->FU[i] = 0.0;
            problem->FV[i] = 0.0; 
        }

        femShallowAddIntegralsElements(problem);
        femShallowAddIntegralsEdges(problem);
        femShallowMultiplyInverseMatrix(problem);
        
        for (i=0; i < problem->size; i++) {
            E[i] += dt * problem->FE[i];
            U[i] += dt * problem->FU[i];
            V[i] += dt * problem->FV[i]; 
        } 
        
        if ((it+1)%sub == 0)
            tsunamiWriteFile(baseResultName,(it+1),U,V,E,problem->nElem,3); 
    }
    /*
        if (it == 0)
        for (i = 0; i < problem->nElem; i++) {
            for (j = 0; j < 3; j++) {
                fprintf(uf, "%13e,", problem->FU[i*3+j]);
                fprintf(ef, "%13e,", problem->FE[i*3+j]);
                fprintf(vf, "%13e,", problem->FV[i*3+j]);
            }
            fprintf(uf, "\n");
            fprintf(ef, "\n");
            fprintf(vf, "\n");
        }
        */

    fclose(uf);
    fclose(ef);
    fclose(vf);
    
    free(E);
    free(problem->FE);
    free(nodes);
    free(elems);
    free(edges);
    free(problem);
 
}
