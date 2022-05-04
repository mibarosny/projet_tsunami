# include "tsunami.h"
// FOR RAPPORT -> a ajouter dans le rapport

/**
 * FOR RAPPORT, on explique la structure particulière de interpolate
 */
double interpolate(double *phi, double *U, int *map, int e, int p, int n)
{
    double u = 0.0; int i;
    for (i=0; i < n; i++)
        u += phi[i]*U[e*map[i]+p];
    return u;
}

void femShallowTriangleMap(int index, int map[3])
{
    int j;
    for (j=0; j < 3; ++j) 
        map[j] = index*3 + j; 
}

void femShallowEdgeMap(int * edges, int * elems, int nElem, int index, int map[2][2])
{
    int i,j,k;   
    for (j=0; j < 2; ++j) {
        int node = edges[index*4+j]; 
        for (k=0; k < 2; k++) {
            int el = edges[index*4+2+k]; 
            map[k][j] = nElem*3;
            if (el >= 0) {
                for (i=0; i < 3; i++) {
                    if (elems[el*3 + i] == node) {
                        map[k][j] = el*3 + i;  }}}}}
}

void tsunamiInitializeCompute(int nElem, double * nodes, int *elems, double *E, double *U, double *V) {
    for (int el = 0; el < nElem; el++)
        for (int i = 0; i < 3; i++) {
            E[el*3+i] = tsunamiInitialConditionOkada(nodes[3*elems[el*3+i]], nodes[3*elems[el*3+i]+1]);
            U[el*3+i] = 0.;
            V[el*3+i] = 0.;
        }
}

void femShallowAddIntegralsElements(int nElem, int* elems, double *nodes, double *E, double *U, double *V, double *FE, double *FU, double *FV) {
    double  xLoc[3],yLoc[3],phi[3],dphidx[3],dphidy[3];
    double  xsi,eta,weight,jac;
    double  x,y,h,u,e,v;
    int map[3], el, i, j;

    for (el = 0; el < nElem; el++) {
        femShallowTriangleMap(el, map);
        int *mapCoord = elems+el*3;
        for (j=0; j < 3; ++j) {
            xLoc[j] = elems[mapCoord[j]];
            yLoc[j] = elems[mapCoord[j]+1];                
        }
        jac = (xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]);
        dphidx[0] = (yLoc[1] - yLoc[2])/jac;
        dphidx[1] = (yLoc[2] - yLoc[0])/jac;
        dphidx[2] = (yLoc[0] - yLoc[1])/jac;
        dphidy[0] = (xLoc[2] - xLoc[1])/jac;
        dphidy[1] = (xLoc[0] - xLoc[2])/jac;
        dphidy[2] = (xLoc[1] - xLoc[0])/jac;

        for (j = 0; j < 3; j++) {
            xsi = gaussTriangleXsi[j];
            eta = gaussTriangleEta[j];
            weight = gaussTriangleWeight[j];
            // Fem discrete
            phi[0] = 1.-xsi-eta;
            phi[1] = xsi;
            phi[2] = eta;
            y = interpolate(phi, nodes, mapCoord, 3, 1, 3);
            x = interpolate(phi, nodes, mapCoord, 3, 0, 3);
            h = interpolate(phi, nodes, mapCoord, 3, 2, 3);
            e = interpolate(phi, E, map, 1, 0, 3);
            u = interpolate(phi, U, map, 1, 0, 3);
            v = interpolate(phi, V, map, 1, 0, 3);
            //if (el == 0 && j == 0)
            //printf("%f %f %f %f %f %f \n", y, x, h, e, u, v);
            double z3d = R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
            double lat = asin(z3d/R)*180/PI;
            double coriolis = 2*Omega*sin(lat);                
            for (i = 0; i < 3; i++) {
                FE[map[i]] += ( h*(dphidx[i]*u + dphidy[i]*v) * ((4*R*R+x*x+y*y)/(4*R*R)) 
                                + phi[i]*h*((x*u+y*v)/(R*R)) )*jac*weight;
                FU[map[i]] += ( phi[i]*(coriolis*v - Gamma*u) + dphidx[i]*g*e*((4*R*R+x*x+y*y)/(4*R*R))
                                + phi[i]*(g*x*e)/(2*R*R))*jac*weight;
                FV[map[i]] += ( phi[i]*(-coriolis*u - Gamma*v) + dphidy[i]*g*e*((4*R*R+x*x+y*y)/(4*R*R))
                                + phi[i]*(g*y*e)/(2*R*R))*jac*weight;
            }

        }
    }
}

void femShallowAddIntegralsEdges(int nElem, int nEdge, int*edges, int size, int* elems, double *nodes, double *E, double *U, double *V, double *FE, double *FU, double *FV) {
    double  xLoc[3],yLoc[3],phi[3],dphidx[3],dphidy[3];
    double  xsi,eta,weight,jac;
    double  x,y,h,u,e,v;
    int map[3], el, i, j;
    double  xEdge[2],yEdge[2],phiEdge[2];
    double  eL,eR,uL,uR,vL,vR,unL,unR;
    double  qe,qu,qv;
    int     edge,mapEdge[2][2];

    for (edge = 0; edge < nEdge; edge++) {
        femShallowEdgeMap(edges, elems, nElem, edge, mapEdge);
        int *mapCoord = edges + edge*4;
        for (j = 0; j < 2; ++j) {
            int node = edges[edge*4+j];
            xEdge[j] = nodes[node*3];
            yEdge[j] = nodes[node*3+1];
        }

        int boundary = (mapEdge[1][0] == size-1);
        double dxdxsi = (xEdge[1] - xEdge[0]);
        double dydxsi = (yEdge[1] - yEdge[0]);
        double norm = sqrt(dxdxsi*dxdxsi + dydxsi*dydxsi);
        double nx =  dydxsi/norm;
        double ny = -dxdxsi/norm;
        jac = norm / 2.0;
        for (j = 0; j < 2; j++) {
            xsi = gaussEdgeXsi[j];
            weight = gaussEdgeWeight[j];
            // femdiscrete
            phiEdge[0] = (1.0 - xsi)/2.0;
            phiEdge[1] = (1.0 + xsi)/2.0;

            eL = interpolate(phiEdge,E,mapEdge[0],1,0,2);
            eR = boundary ? eL : interpolate(phiEdge,E,mapEdge[1],1,0,2);
            uL = interpolate(phiEdge,U,mapEdge[0],1,0,2);
            uR = interpolate(phiEdge,U,mapEdge[1],1,0,2);
            vL = interpolate(phiEdge,V,mapEdge[0],1,0,2);
            vR = interpolate(phiEdge,V,mapEdge[1],1,0,2);
            y = interpolate(phiEdge,nodes,mapCoord,3,1,2);
            x = interpolate(phiEdge,nodes,mapCoord,3,0,2);

            printf("%f %f \n", x, y);
            unL = uL*nx+ vL*ny;
            unR = boundary ? -unL : uR*nx + vR*ny;
            qe =  0.5*h*   ( (unL+unR) + sqrt(g/h)*( eL-eR ) )*((4*R*R+x*x+y*y)/(4*R*R));
            qu =  0.5*g*nx*( ( eL+eR ) + sqrt(h/g)*(unL-unR) )*((4*R*R+x*x+y*y)/(4*R*R));
            qv =  0.5*g*ny*( ( eL+eR ) + sqrt(h/g)*(unL-unR) )*((4*R*R+x*x+y*y)/(4*R*R));   
            for (i = 0; i < 2; i++) {
                FE[mapEdge[0][i]] -= qe*phiEdge[i]*jac*weight; 
                FU[mapEdge[0][i]] -= qu*phiEdge[i]*jac*weight; 
                FV[mapEdge[0][i]] -= qv*phiEdge[i]*jac*weight; 
                FE[mapEdge[1][i]] += qe*phiEdge[i]*jac*weight;
                FU[mapEdge[1][i]] += qu*phiEdge[i]*jac*weight;
                FV[mapEdge[1][i]] += qv*phiEdge[i]*jac*weight;
            }
        }
    }

}

void tsunamiComputeRightHandSide(int nElem, int nEdge, int*edges, int size, int* elems, double *nodes, double *E, double *U, double *V, double *FE, double *FU, double *FV) {
    femShallowAddIntegralsElements(nElem, elems, nodes, E, U, V, FE, FU, FV);
    femShallowAddIntegralsEdges(nElem, nEdge, edges, size, elems, nodes, E, U, V, FE, FU, FV);
}

void femShallowMultiplyInverseMatrix(int nElem, int* elems, double *nodes, double *E, double *U, double *V, double *FE, double *FU, double *FV) {
    double  xLoc[3],yLoc[3],phi[3],dphidx[3],dphidy[3];
    double  xsi,eta,weight,jac;
    double  x,y,h,u,e,v;
    int map[3], el, i, j;
    double invA[3][3] = {{18.0,-6.0,-6.0},{-6.0,18.0,-6.0},{-6.0,-6.0,18.0}};
    double BEloc[3],BUloc[3],BVloc[3];

    for (el=0; el < nElem; el++) {
        femShallowTriangleMap(el, map);
        int *mapCoord = elems+el*3;
        for (j=0; j < 3; ++j) {
            xLoc[j] = elems[mapCoord[j]];
            yLoc[j] = elems[mapCoord[j]+1];                
        }
        jac = (xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]);
        for (i=0; i < 3; i++) {
            BEloc[i] = FE[map[i]];
            BUloc[i] = FU[map[i]];
            BVloc[i] = FV[map[i]];
            FE[map[i]] = 0.0; 
            FU[map[i]] = 0.0; 
            FV[map[i]] = 0.0; }
        for (i=0; i < 3; i++) { 
            for (j=0; j < 3; j++) {
                FE[map[i]] += invA[i][j] * BEloc[j] / jac; 
                FU[map[i]] += invA[i][j] * BUloc[j] / jac; 
                FV[map[i]] += invA[i][j] * BVloc[j] / jac;
            }
        }
    }
}


void tsunamiUpdateRungeKutta(double dt, int nElem, int nEdge, int*edges, int size, int* elems, double *nodes, double *E, double *U, double *V, double *FE, double *FU, double *FV) {

    int i, j;
    double *Uold, *Eold, *Vold, *Upred, *Epred, *Vpred;
    Uold        = malloc(sizeof(double) * size * 6);
    Eold        = Uold + size;
    Vold        = Eold + size;
    Upred       = Vold + size;
    Epred       = Upred + size;
    Vpred       = Epred + size;
    double *KE, *KU, *KV;
    KE = malloc(sizeof(double) * size*3);
    KU = KE + size;
    KV = KU + size;
    for (i=0; i < size; i++) {
        KU[i] = 0.0;
        KE[i] = 0.0;
        KV[i] = 0.0;
        Uold[i] = U[i];  
        Eold[i] = E[i];  
        Vold[i] = V[i];  
    }

    const int    nStage   = 2;
    const double beta[2]  = {0.0,  1.0  }; 
    const double gamma[2] = {1./1., 1./2.};
    
    for (j = 0; j < nStage; j++) {
        for (i=0; i < size; i++) {
                Upred[i] = Uold[i] + dt * beta[j] * KU[i];
                Vpred[i] = Vold[i] + dt * beta[j] * KV[i];
                Epred[i] = Eold[i] + dt * beta[j] * KE[i];
            }
            //legendreComputeRightHandSide(nElem, nEdge, edges, size, elems, nodes, KE, KU, KV, FE, FU, FV);
            for (i=0; i < size; i++) {
        	    U[i] += dt * gamma[j] * KU[i]; 
        	    V[i] += dt * gamma[j] * KV[i]; 
        	    E[i] += dt * gamma[j] * KE[i]; 
            }
        }
    
    free(Uold);
    free(KE);
}

void tsunamiCompute(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName)
{ 
    int i,j,el,nNode,nElem,nEdge,trash, size;
    double dtrash;
    
    double BathMax = 9368;

    double *nodes;               // coordonnés X, Y et bathymétrie 
    int * elems;                 // noeud 1, 2 et 3
    int * edges;                 // noeud 1 et 2, elems gauche, droite

    /** lecture du fichier */
    
    FILE* file = fopen(meshFileName,"r");
    fscanf(file, "Number of nodes %d \n",&nNode);   
    /* FOR RAPPORT En fait, l’appel à malloc est coûteux en temps et il vaut mieux déclarer une fois un tableau de
    nNode*3, qui va ensuite être virtuellement découpé en plusieurs, que faire l’inverse. Ici, on va
    avoir besoin de nNOde éléments pour les vecteurs X, Y et bath, on préfère
    donc faire qu’un seul appel à malloc.
    
    Bla bla de cache aussi  */
    nodes = malloc(sizeof(double)*nNode*3);
    for (i = 0; i < nNode; i++) 
        fscanf(file,"%d : %le %le %le\n",&trash,&nodes[i*3],&nodes[i*3+1],&nodes[i*3+2]); 
        
    fscanf(file, "Number of triangles %d \n",&nElem); 
    elems = malloc(sizeof(int)*3*nElem);
    size = 3*nElem + 1;
    for (i = 0; i < nElem; i++) 
        fscanf(file,"%d : %d %d %d \n",&trash,&elems[i*3],&elems[i*3+1],&elems[i*3+2]);   

    fscanf(file, "Number of edges %d \n",&nEdge); 
    edges = malloc(sizeof(int)*4*nEdge);
    for (i = 0; i < nEdge; i++) 
        fscanf(file,"%d : %d %d %d %d \n",&trash,&edges[i*4],&edges[i*4+1],&edges[i*4+2],&edges[i*4+3]);   

    fclose(file); 

    double *E, *U, *V;  // élévation, vitesse horizontal et vertical 
    E = malloc(sizeof(double)*(size)*3);    // dégré de liberté
    U = E + size;
    V = U + size;    

    double *FU, *FE, *FV;
    FE = malloc(sizeof(double)*(size)*3);
    FU = FE + size;
    FV = FU + size;

    // Conditions initiales
    tsunamiInitializeCompute(nElem, nodes, elems, E, U, V);

    for (int it = 0; it < nmax; it++) {

        for (i=0; i < size; i++) {
            FE[i] = 0.0;
            FU[i] = 0.0;
            FV[i] = 0.0; 
        }
        femShallowAddIntegralsElements(nElem, elems, nodes, E, U, V, FE, FU, FV);
        femShallowAddIntegralsEdges(nElem, nEdge, edges, size, elems, nodes, E, U, V, FE, FU, FV);
        femShallowMultiplyInverseMatrix(nElem, elems, nodes, E, U, V, FE, FU, FV);
        for (i=0; i < size; i++) {
            E[i] += dt * FE[i];
            U[i] += dt * FU[i];
            V[i] += dt * FV[i]; 
        } 
        
        if ((it+1)%sub == 0)
            tsunamiWriteFile(baseResultName,(it+1),U,V,E,nElem,3); 
    
    }
    
    free(E);
    free(elems);
    free(nodes);
    free(edges);
 
}
