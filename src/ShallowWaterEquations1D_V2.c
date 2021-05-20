#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <stdbool.h>
#define PI 3.14
#define max(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})

#define NUM_VARS 2
#define NX 100
#define N_GHOST 2
#define MIN_X (0.0)
#define MAX_X (1.0)

#define CFL 0.8
#define STOPPING_TIME 3

#define NUM_RK_STEPS 1
#define MAX_TIME_ITER 1000000
#define NUMERICAL_METHOD 1
#define INITIAL_CONDITIONS 1
#define BOUNDARY_CONDITIONS 2



 void InitialConditions(double** u_array, double** u_Exact, double** u_init, double* x_centroid_array, double* dx_array, int rkStep);
 double timestepCalculator(double** u_array, double* dx_array,  int rkStep);
 void reconstructVariables(double** u_array, double** uEast, double** uWest,  double* dx_array, int rkStep);
 void updateGhostCells(double** u_array, int rkStep);
 void fluxCalculator(double** u_array, double** uEast, double** uWest,  double** total_flux, double* dx_array);
void TimeIntegration(double** u_array, double** u_Exact, double* dx_array,  double** total_flux, double dt, int rkStep, FILE* IC2_fp);
 void copySolution(double** u_array, double** u_Exact, int rkStep, FILE* IC2_fp1);
void writeSolutionToFile(double** u_array );
void ErrorCalculator(double** u_array, double** u_init);

int main()
{
    //open files to hold initial condition data
    FILE* IC_fp = NULL;
    IC_fp = fopen("./out/Square_Initial_conditions.txt", "w");

    FILE* IC2_fp = NULL;
    IC2_fp = fopen("./out/Sine_Initial_Conditions.txt", "w");

    FILE* IC2_fp1 = NULL;
    IC2_fp1 = fopen("./out/solution.txt", "w");
    FILE* IC3_fp = NULL;
    IC3_fp = fopen("./out/Aolution.txt", "w");

    //error handle
    if (IC_fp == NULL || IC2_fp == NULL)
    {
        printf("\n");
        printf("FILE OPEN FAILED");
        printf("\n");
    }

    double DX = (MAX_X - MIN_X) / NX;
    int ARRAY_SIZE_X = NX + 2 * N_GHOST;


    double** u_array = (double**)calloc(ARRAY_SIZE_X, sizeof(double*));
    double** u_Exact = (double**)calloc(ARRAY_SIZE_X, sizeof(double*));
    double** u_init = (double**)calloc(ARRAY_SIZE_X, sizeof(double*));

    double** uEast = (double**)calloc(ARRAY_SIZE_X, sizeof(double*));
    double** uWest = (double**)calloc(ARRAY_SIZE_X, sizeof(double*));

    double* dx_array = (double*)calloc(ARRAY_SIZE_X, sizeof(double));
    double* x_centroid_array = (double*)calloc(ARRAY_SIZE_X, sizeof(double));

    double** total_flux = (double**)calloc(ARRAY_SIZE_X, sizeof(double*));

    for (int i = 0; i < ARRAY_SIZE_X; i++)
    {
        u_array[i] = (double*)calloc(NUM_VARS, sizeof(double));
        u_Exact[i] = (double*)calloc(NUM_VARS, sizeof(double));
        u_init[i] = (double*)calloc(NUM_VARS, sizeof(double));

        uEast[i] = (double*)calloc(NUM_VARS, sizeof(double));
        uWest[i] = (double*)calloc(NUM_VARS, sizeof(double));
        total_flux[i] = (double*)calloc(NUM_VARS, sizeof(double));


    }

    //populate DX ann cx
    for (int i = 0; i < ARRAY_SIZE_X; i++ )
    {

        dx_array[i] = DX;

        x_centroid_array[i] = MIN_X + DX * (i + 0.5 - N_GHOST);



    }

    double INITIAL_VELOCITY = 1.0;
    double time = 0.0;
    int rkStep = 0;
    bool lastTimeStep = false;
    
    double error = 0.0;
    InitialConditions(u_array, u_Exact, u_init, x_centroid_array,  dx_array, rkStep);

    for (int t = 0; t < MAX_TIME_ITER; t++)
     {
         //printf("\n\n\n");
         double dt =  timestepCalculator(u_array, dx_array, rkStep);

        for (int rkStep = 0; rkStep < NUM_RK_STEPS; rkStep++)
        {
            updateGhostCells(u_array, rkStep);
            reconstructVariables(u_array, uEast, uWest,  dx_array, rkStep);
            fluxCalculator(u_array, uEast, uWest, total_flux, dx_array);
            TimeIntegration(u_array,u_Exact, dx_array,total_flux, dt, rkStep, IC2_fp);
            copySolution(u_array, u_Exact, rkStep, IC2_fp1);
        }
        
       

        if(t % 3 == 0)
        {
            //writeSolutionToFile(u_array);

        }
        //ErrorCalculator(u_array, u_init);
        for (int i = N_GHOST; i < NX + N_GHOST; i++) {
        
            error += fabs(u_array[i][0] - u_init[i][0]);

            // errorL1 += error;
            // errorL2 += error * error;
            // errorLinf += fmax(error, errorLinf);
        }

        if (time + dt > STOPPING_TIME)
        {
            dt = STOPPING_TIME - time;
            lastTimeStep = true;
        }


        time += dt;
        if (lastTimeStep)
        {
            break;
        }

    }

    error /= NX;

     printf("\n\n%f\n", error);
     //ErrorCalculator(u_array, u_init);
    
    


    return EXIT_SUCCESS;
}




void InitialConditions(double** u_array, double** u_Exact, double** u_init, double* x_centroid_array,  double* dx_array, int rkStep)
{

    switch (INITIAL_CONDITIONS) {

        case 1:

            for (int i = 0; i < NX + 2 * N_GHOST; i++)
            {

                double rSqr = x_centroid_array[i] * x_centroid_array[i];
                double Gaussian = 1.0;
                double variance = 0.2 * 0.2;
                double initvars[3] = {Gaussian * exp(-rSqr / 2.0 /variance), 0.0, 0.0};

                for (int var = 0; var < NUM_VARS; var++)
                {
                    u_array[i][var] = initvars[var];
                    u_Exact[i][var] = initvars[var];
                    u_init[i][var] = initvars[var];
                        //printf(" %14f", u_Exact[i][j][0]);
                        //printf("%f\n", Gaussian * exp(-rSqr / 2.0 /variance));
                        //fprintf(IC_fp, "%f\n", u_array[i][j][rkStep][var]);
                }

                    //fprintf(IC_fp, "%f\n", u_array[i][j][0]);
                    
            }
            break;
        
        case 2:

            for (int i = N_GHOST; i < NX + N_GHOST; i++)
            {
                if (x_centroid_array[i] < 0.5)
                {
                    u_array[i][0] = 1.0;
                    u_array[i][1] = 0.0;

                    u_Exact[i][0] = 1.0;
                    u_Exact[i][1] = 0.0;

                    u_init[i][0] = 1.0;
                    u_init[i][1] = 0.0;

                }
                else{
                    u_array[i][0] = 0.5;
                    u_array[i][1] = 0.0;

                    u_Exact[i][0] = 0.5;
                    u_Exact[i][1] = 0.0;

                    u_init[i][0] = 0.5;
                    u_init[i][1] = 0.0;
                }
                //printf("%f\n", h_array[i][rkStep]);
                //printf("%f %f\n", u_array[i].cx, u_array[i].U[rkStep][0]);
            }
            break;
    
    }

}

double timestepCalculator(double** u_array, double* dx_array, int rkStep)
{
    double MIN_DT = (double)(INFINITY);
    for (int i = N_GHOST; i < NX + 2 * N_GHOST - N_GHOST; i++)
    {

        double h;
        double u;
        double v;
        double dtx;
        double dty;

        h = u_array[i][0];
        u = u_array[i][1] / h;

        double gravity = 9.81;
        double c = sqrt(gravity * h);
        dtx = dx_array[i] / (fabs(u) + c);

        double DT_AT_CELL = 1.0 / (1.0 /dtx);
        if (DT_AT_CELL < MIN_DT)
        {
            MIN_DT = DT_AT_CELL;
        }


    }

    //printf("%f\n", MIN_DT * CFL);
    return MIN_DT * CFL;
}

void updateGhostCells(double** u_array, int rkStep)
{

        switch (BOUNDARY_CONDITIONS) {

            case 1: 

                for (int ghostcell = 0; ghostcell < N_GHOST; ghostcell++)
                {
                    for (int var = 0; var < NUM_VARS; var++)
                    {
                        //we want to apply periodic boundary conditions. Whereby at the end of wach
                        //timestep we copy the last two u_array to the first ghost u_array, and the first two ghost u_array
                        //with the last two ghost u_array
                        u_array[ghostcell][var] = u_array[NX + ghostcell][var];
                        u_array[NX + N_GHOST + ghostcell][var] = u_array[N_GHOST + ghostcell][var];
                        //printf("%f\n", u_array[ghostcell][rkStep]);
                    }

                }
                break;
            
            case 2:

                for (int i = 0; i < N_GHOST; i++)
                {
                    int var = 0;
                    u_array[i][var] = u_array[2 * N_GHOST - i - 1][var];
                    u_array[NX + N_GHOST + i][var] = u_array[NX + N_GHOST - i - 1][var];

                    var = 1;
                    u_array[i][var] = -u_array[2 * N_GHOST - i - 1][var];
                    u_array[NX + N_GHOST + i][var] = -u_array[NX + N_GHOST - i - 1][var];
                }

              
        }
        


    return;
}

void reconstructVariables(double** u_array, double** uEast, double** uWest, double* dx_array, int rkStep)
{
    //int SLOPE_LIMITERS = 2; Once we have applied the BCs we want to reconstruct our variables
        switch(NUMERICAL_METHOD)
        {
            case 1:
                //we want to index from one cell before the internal u_array to one cell after the size of the mesh
                //so that we can account for the fluxes at i - 1/2 and i + 1/2 for the first and last u_array
                for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                {

                    for (int var = 0; var < NUM_VARS; var++)
                    {
                        //first order reconstruction we copy the centroid values from u_array into our i+1/2 and i-1/2 arrays
                        //which are represented by uEast and uWest. uEast and uWest reprresents the interface values of u_array
                        uWest[i][var] = u_array[i][var];
                        uEast[i][var] = u_array[i][var];


                    }
                      //printf("%f\n", uNorth[i][j][2]);


                }
                break;
            
            case 2:
                //we want to index from one cell before the internal u_array to one cell after the size of the mesh
                //so that we can account for the fluxes at i - 1/2 and i + 1/2 for the first and last u_array
                for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                {

                    for (int var = 0; var < NUM_VARS; var++)
                    {
                        double r = (u_array[i][var] - u_array[i - 1][var] + 1e-6) / (u_array[i + 1][var] - u_array[i][var] + 1e-6);
                        r = fmax(0, r);
                        double phi = (r * r + r) / (r * r + 1.0);
                        double du_dx = (u_array[i + 1][var] - u_array[i][var]) / 2.0 / dx_array[i];

                        uWest[i][var] = u_array[i][var] - phi *  du_dx * dx_array[i] / 2.0;
                        uEast[i][var] = u_array[i][var] + phi * du_dx * dx_array[i] / 2.0;


                    }
                      //printf("%f\n", uNorth[i][j][2]);


                }
                break;


         }
    return;

}

void fluxCalculator(double** u_array, double** uEast, double** uWest,  double** total_flux, double* dx_array)
{
    for (int i = 0; i < NX + 2 * N_GHOST; i++)
    {

        for (int var = 0; var < NUM_VARS; var++)
        {
            total_flux[i][var] = 0.0;
             //printf("%f", total_flux[i][j][var]);
        }

    }

    double dx = dx_array[0];


    //we want to look through all of the internal u_array



    //double flux;
    int i;
    for (i = N_GHOST; i < NX + N_GHOST + 1; i++) {
        double* UL = uEast[i - 1];
        double* UR = uWest[i];

        double gravity = 9.81;
        double hL = UL[0];
        double uL = UL[1] / hL;
                    //double vL = leftValue[2] / hL;
        double cL = sqrt(gravity * hL);
        double maxSpeedL = fabs(uL) + cL;

        double hR = UR[0];
        double uR = UR[1] / hR;
                    //double vL = leftValue[2] / hL;
        double cR = sqrt(gravity * hR);
        double maxSpeedR = fabs(uR) + cR;

        double maxWaveSpeed = fmax(maxSpeedL, maxSpeedR);

        double FL[2] = {
                 hL * uL,
                 hL * uL * uL + 0.5 * gravity * hL * hL};
        double FR[2] = {
                 hR * uR,
                 hR * uR * uR + 0.5 * gravity * hR * hR};

        int var;
        for (var = 0; var < NUM_VARS; var++) {

            double flux = 0.5 * (FL[var] + FR[var]) - 0.5 * maxWaveSpeed * (UR[var] - UL[var]);
            total_flux[i - 1][var] += flux;
            total_flux[i][var] -= flux;
        }
    }



   // printf("%f\n", leftValue[1]);




    return;
}



void TimeIntegration(double** u_array, double** u_Exact, double* dx_array,  double** total_flux, double dt, int rkStep, FILE* IC2_fp)
{
    //printf("%f\n", area);

    // for (int i = N_GHOST; i < NX + N_GHOST; i++)
    // {

    //     for (int var = 0; var < NUM_VARS; var++)
    //     {
    //         u_Exact[i][var] = u_array[i][var] + (-dt / dx_array[0]) * total_flux[i][var];
    //                 //printf("hello");
    //         //printf(" %14f", u_Exact[i][var]);
    //                 //fprintf(IC2_fp, "%f\n", u_Exact[i][j][rkStep][var]);

    //                 //u_array[i][j][rkStep + 1] = 0.5 * (u_array[i][j - 1][rkStep] + u_array[i][j + 1][rkStep]) - 0.5 * (dt / area) * total_flux[i][j];

    //     }
    // }


    double dx = dx_array[0];


    switch(NUM_RK_STEPS)
    {
        case 1:
            for (int i = N_GHOST; i < NX + N_GHOST; i++)
            {
                for (int var = 0; var < NUM_VARS; var++)
                {
                    u_Exact[i][var] = u_array[i][var] + (-dt / dx_array[0]) * total_flux[i][var];
                }

            }
            break;

        case 2:
            for (int i = N_GHOST; i < NX + N_GHOST; i++)
            {
                for (int var = 0; var < NUM_VARS; var++)
                {


                    switch (rkStep)
                    {
                        case 0:
                            u_Exact[i][var] = u_array[i][var] + (-dt / dx_array[0]) * total_flux[i][var];
                            break;

                        case 1:
                            u_Exact[i][var] = 0.5 * (u_array[i][var] + u_array[i][var] + (-dt / dx_array[0]) * total_flux[i][var]);
                            //u_array[i].U[rkStep + 1][var] = 0.5 * (u_array[i].U[rkStep - 1][var] + u_array[i].U[rkStep][var] + (-dt / dx)  * u_array[i].totalFlux[var]);
                            //u_array[i][rkStep + 1] = 0.5 * (u_array[i][rkStep - 1] + u_array[i][rkStep] + dt / dx_array[i] * total_flux[i]);
                            break;
                        default:
                            printf("Timestep not available");
                    }
                }
            }
            break;
        case 3:
            for (int i = N_GHOST; i < NX + N_GHOST; i++)
            {
                for (int var = 0; var < NUM_VARS; var++)
                {


                    switch(rkStep)
                    {
                        case 0:
                            u_Exact[i][var] = 0.5 * (u_array[i][var] + u_array[i][var] + (-dt / dx_array[0]) * total_flux[i][var]);
                            
                            break;

                        case 1:
                            u_Exact[i][var] = 3.0 / 4.0 * u_array[i][var] + 1.0 / 4.0 * u_array[i][var] + 1.0 / 4.0 * (-dt / dx_array[0]) * total_flux[i][var];
                            // u_array[i].U[rkStep + 1][var] = 3.0 / 4.0 * u_array[i].U[rkStep - 1][var] + 1.0/ 4.0 * u_array[i].U[rkStep][var]
                            //                             + 1.0 / 4.0 * (-dt / dx)  * u_array[i].totalFlux[var];
                            //u_array[i][rkStep + 1] = 3.0 / 4.0 * u_array[i][rkStep - 1] + 1.0 / 4.0 * u_array[i][rkStep] + 1.0 / 4.0 * dt / dx_array[i] * total_flux[i];
                            break;

                        case 2:
                            u_Exact[i][var] = 1.0 / 3.0 * u_array[i][var] + 2.0 / 3.0 * u_array[i][var] + 2.0 / 3.0 * (-dt / dx_array[0]) * total_flux[i][var];
                            // u_array[i].U[rkStep + 1][var] = 1.0 / 3.0 * u_array[i].U[rkStep - 2][var] + 2.0 / 3.0 * u_array[i].U[rkStep][var]
                            //                             + 2.0 / 3.0 * (-dt / dx)  * u_array[i].totalFlux[var];
                             //u_array[i][rkStep + 1] = 1.0 / 3.0 * u_array[i][rkStep - 2] + 2.0 / 3.0 * u_array[i][rkStep] + 2.0 / 3.0 * dt / dx_array[i] * total_flux[i];
                            break;

                    }
                }
            }
        break;
    }


    return;
}

void copySolution(double** u_array, double** u_Exact, int rkStep, FILE* IC2_fp1)
{
    for (int i = N_GHOST; i < NX + N_GHOST; i++)
    {
        for (int var = 0; var < NUM_VARS; var++)
        {
            u_array[i][var] = u_Exact[i][var];




          //printf(" %14f", u_Exact[i][j][0]);



        }
         //printf(" %14f", u_array[i][0]);
    }

}

void writeSolutionToFile(double** u_array) {
    printf("\nWriting to file... ");
    //printf("\nWriting to file... ");
    static int fileIndex = 0;
    char fileName[100];


    sprintf(fileName, "sol_%d.txt", fileIndex);
    FILE* file = fopen(fileName, "w");


    // Write cell data
    const int rkStep = 0;
    int i;

    for (i = N_GHOST; i < NX + N_GHOST; i++)
    {

        for (int var = 0; var < NUM_VARS; var++)
        {



        }
        fprintf(file, "%f\n", u_array[i][0]);


    }

    fclose(file);

    fileIndex++;
    //printf("done.\n");
    printf("done.\n");
}

void ErrorCalculator(double** u_array, double** u_init) {

    double errorL1 = 0.0;
    double errorL2 = 0.0;
    double errorLinf = 0.0;

    for (int i = N_GHOST; i < NX + N_GHOST; i++) {
        
        double error = fabs(u_array[i][0] - u_init[i][0]);

         errorL1 += error;
         errorL2 += error * error;
         errorLinf += fmax(error, errorLinf);
    }

    errorL1 /=  NX;
    // errorL2 = sqrt(errorL2 / NX);

    printf("\n\n%f \n", errorL1);
}


     