#include <stdio.h>
#include <stdlib.h>
// #include "config.h"
// #include "data.h"
#include "math.h"
#include <stdbool.h>


#define NUM_VARS 2
#define NUM_CELLS 100
#define NUM_GHOST_CELLS 2
#define MIN_X (0.0)
#define MAX_X (1.0)

#define CFL 0.6
#define STOPPING_TIME 2

#define NUM_RK_STEPS 1
#define MAX_TIME_ITER 1000000

struct Cell {
        double U[NUM_RK_STEPS + 1][NUM_VARS]; // Average Conservative variables
        double UWest[NUM_VARS]; // Reconstructed values
        double UEast[NUM_VARS]; // Reconstructed values
        double totalFlux[NUM_VARS];

        double dx; // width of cell
        double cx; // centroid of cell
    };
typedef struct Cell Cell;

void generateMesh(Cell cells[]);
void setupInitialConditions(Cell cells[]);
double calculateTimeStep(Cell cells[]);
void updateGhostCells(Cell cells[], int rkStep);
void reconstructSolution(Cell cells[], int rkStep);
void calculateFlux(Cell cells[], int rkStep);
void TimeIntegration(Cell cells[], int rkStep, double dt);
void CopySolution(Cell cells[], FILE* IC_fp);
void writeSolutionToFile(Cell cells[]);



int main(int argc, char** argv) {

    FILE* IC_fp = NULL;
    IC_fp = fopen("Sol.txt", "w");


    Cell cells[NUM_CELLS + 2 * NUM_GHOST_CELLS];
    generateMesh(cells);


    setupInitialConditions(cells);


    bool lastTimeStep = false;
    double time = 0.0;
    int rkStep = 0;
    for (int timeIter = 0; timeIter < MAX_TIME_ITER; timeIter++) {


        double dt = calculateTimeStep(cells);
        if (time + dt > STOPPING_TIME) {

            dt = STOPPING_TIME - time;
            lastTimeStep = true;
        }


        int rkStep;
        for (rkStep = 0; rkStep < NUM_RK_STEPS; rkStep++) {

            updateGhostCells(cells, rkStep);

            reconstructSolution(cells, rkStep);

            calculateFlux(cells, rkStep);

            TimeIntegration(cells, rkStep, dt);
        }


        CopySolution(cells, IC_fp);


        time = time + dt;

        if (lastTimeStep) {
            break;
        }
        if ((timeIter + 1) % 15 == 0){
            writeSolutionToFile(cells);
        }
     }

    //writeSolutionToFile(cells);
    return (EXIT_SUCCESS);
}



void generateMesh(Cell cells[]) {
    printf("\nSetting up the cell geometry data...");
    int i;
    double dx = (MAX_X - MIN_X) / NUM_CELLS;
    double half_dx = dx / 2.0;
    for (i = -NUM_GHOST_CELLS; i < NUM_CELLS + NUM_GHOST_CELLS; i++) {
        struct Cell* cell = (cells + i + NUM_GHOST_CELLS);
        cell -> cx = MIN_X + i * dx + half_dx;
        cell -> dx = dx;
    }
    printf(" done.");
}

void setupInitialConditions(Cell cells[]) {
    printf("\nInitializing solution in the domain...");

   int rkStep = 0;
    for (int i = NUM_GHOST_CELLS; i < NUM_CELLS + NUM_GHOST_CELLS; i++)
    {
        if (cells[i].cx < 0.5)
        {
            cells[i].U[rkStep][0] = 1.0;
            cells[i].U[rkStep][1] = 0.0;

        }
        else{
            cells[i].U[rkStep][0] = 0.5;
            cells[i].U[rkStep][1] = 0.0;
        }
        //printf("%f\n", h_array[i][rkStep]);
        //printf("%f %f\n", cells[i].cx, cells[i].U[rkStep][0]);
    }
    printf(" done.");
}


double calculateTimeStep(Cell cells[]) {


    double Min_DT = (double)(INFINITY);
    for (int i = NUM_GHOST_CELLS; i < NUM_CELLS + 2 - NUM_GHOST_CELLS; i++)
    {
        int rkStep = 0;
        double h;
        double u;
        double dtx;

        h = cells[i].U[rkStep][0];
        u = cells[i].U[rkStep][1] / h;

        double gravity = 9.81;
        double c = sqrt(gravity * h);
        dtx = cells[i].dx / (fabs(u) + c);
        double DT_AT_CELL = 1.0 / (1.0 / dtx);
        if(DT_AT_CELL < Min_DT ) {
            Min_DT = DT_AT_CELL;
        }
    }



    //printf("%f\n", MIN_DT * CFL);
    return Min_DT * CFL;
}




void updateGhostCells(Cell cells[], int rkStep) {
    // TODO: simple extrapolation as of now, needs to be changed
     for (int i = 0; i < NUM_GHOST_CELLS; i++)
    {
        int var = 0;
        cells[i].U[rkStep][var] = cells[2 * NUM_GHOST_CELLS - i - 1].U[rkStep][var];
        cells[NUM_CELLS + NUM_GHOST_CELLS + i].U[rkStep][var] = cells[NUM_CELLS + NUM_GHOST_CELLS - i - 1].U[rkStep][var];

        var = 1;
        cells[i].U[rkStep][var] = -cells[2 * NUM_GHOST_CELLS - i - 1].U[rkStep][var];
        cells[NUM_CELLS + NUM_GHOST_CELLS + i].U[rkStep][var] = -cells[NUM_CELLS + NUM_GHOST_CELLS - i - 1].U[rkStep][var];
    }
}

void reconstructSolution(Cell cells[], int rkStep) {
    // TODO: Exercise for students - Increase the order of the reconstruction
    int i;
    for (i = NUM_GHOST_CELLS - 1; i < NUM_CELLS + NUM_GHOST_CELLS + 1; i++) {
        // Simple first order reconstruction
        int var;
        for (var = 0; var < NUM_VARS; var++) {
            cells[i].UWest[var] = cells[i].U[rkStep][var];
            cells[i].UEast[var] = cells[i].U[rkStep][var];
        }
    }
}




void calculateFlux(Cell cells[], int rkStep) {
    int i;
    for (i = 0; i < NUM_CELLS + 2 * NUM_GHOST_CELLS; i++) {
        int var;
        for (var = 0; var < NUM_VARS; var++) {
            cells[i].totalFlux[var] = 0.0;
        }
    }



    double flux;
    
    for (i = NUM_GHOST_CELLS; i < NUM_CELLS + NUM_GHOST_CELLS + 1; i++) {
        double* UL = cells[i - 1].UEast;
        double* UR = cells[i].UWest;

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

        double FL[NUM_VARS] = {
                 hL * uL,
                 hL * uL * uL + 0.5 * gravity * hL * hL};
        double FR[NUM_VARS] = {
                 hR * uR,
                 hR * uR * uR + 0.5 * gravity * hR * hR};

        int var;
        for (var = 0; var < NUM_VARS; var++) {

            flux = 0.5 * (FL[var] + FR[var]) - 0.5 * maxWaveSpeed * (UR[var] - UL[var]);
            cells[i - 1].totalFlux[var] += flux;
            cells[i].totalFlux[var] -= flux;
        }
    }
}

void TimeIntegration(Cell cells[], int rkStep, double dt) {


    double dx = cells[0].dx;


    switch(NUM_RK_STEPS)
    {
        case 1:
            for (int i = NUM_GHOST_CELLS; i < NUM_CELLS + NUM_GHOST_CELLS; i++)
            {
                for (int var = 0; var < NUM_VARS; var++)
                {
                    cells[i].U[rkStep + 1][var] = cells[i].U[rkStep][var] + (-dt / dx)  * cells[i].totalFlux[var];
                }

            }
            break;

        case 2:
            for (int i = NUM_GHOST_CELLS; i < NUM_CELLS + NUM_GHOST_CELLS; i++)
            {
                for (int var = 0; var < NUM_VARS; var++)
                {


                    switch (rkStep)
                    {
                        case 0:
                            cells[i].U[rkStep + 1][var] = cells[i].U[rkStep][var] + (-dt / dx)  * cells[i].totalFlux[var];
                            break;

                        case 1:
                            cells[i].U[rkStep + 1][var] = 0.5 * (cells[i].U[rkStep - 1][var] + cells[i].U[rkStep][var] + (-dt / dx)  * cells[i].totalFlux[var]);
                            //u_array[i][rkStep + 1] = 0.5 * (u_array[i][rkStep - 1] + u_array[i][rkStep] + dt / dx_array[i] * total_flux[i]);
                            break;
                        default:
                            printf("Timestep not available");
                    }
                }
            }
            break;
        case 3:
            for (int i = NUM_GHOST_CELLS; i < NUM_CELLS + NUM_GHOST_CELLS; i++)
            {
                for (int var = 0; var < NUM_VARS; var++)
                {


                    switch(rkStep)
                    {
                        case 0:
                            cells[i].U[rkStep + 1][var] = cells[i].U[rkStep][var] + (-dt / dx)  * cells[i].totalFlux[var];
                            break;

                        case 1:
                            cells[i].U[rkStep + 1][var] = 3.0 / 4.0 * cells[i].U[rkStep - 1][var] + 1.0/ 4.0 * cells[i].U[rkStep][var]
                                                        + 1.0 / 4.0 * (-dt / dx)  * cells[i].totalFlux[var];
                            //u_array[i][rkStep + 1] = 3.0 / 4.0 * u_array[i][rkStep - 1] + 1.0 / 4.0 * u_array[i][rkStep] + 1.0 / 4.0 * dt / dx_array[i] * total_flux[i];
                            break;

                        case 2:
                            cells[i].U[rkStep + 1][var] = 1.0 / 3.0 * cells[i].U[rkStep - 2][var] + 2.0 / 3.0 * cells[i].U[rkStep][var]
                                                        + 2.0 / 3.0 * (-dt / dx)  * cells[i].totalFlux[var];
                             //u_array[i][rkStep + 1] = 1.0 / 3.0 * u_array[i][rkStep - 2] + 2.0 / 3.0 * u_array[i][rkStep] + 2.0 / 3.0 * dt / dx_array[i] * total_flux[i];
                            break;

                    }
                }
            }
        break;
    }
}

void CopySolution(Cell cells[], FILE* IC_fp) {
    int i;
    for (i = NUM_GHOST_CELLS; i < NUM_CELLS +  NUM_GHOST_CELLS; i++) {
        int var;
        for (var = 0; var < NUM_VARS; var++) {
            cells[i].U[0][var] = cells[i].U[NUM_RK_STEPS][var];
        }
        fprintf(IC_fp, "%f\n", cells[i].U[0][1]);
        //printf(" %14f", cells[i].U[0][0]);

    }


}

void writeSolutionToFile(Cell cells[]) {
    printf("\nWriting to file... ");
    //printf("\nWriting to file... ");
    static int fileIndex = 0;
    char fileName[100];


    sprintf(fileName, "sol_%d.txt", fileIndex);
    FILE* file = fopen(fileName, "w");


    // Write cell data
    const int rkStep = 0;
    int i;

    for (i = NUM_GHOST_CELLS; i < NUM_CELLS + NUM_GHOST_CELLS; i++)
    {

        for (int var = 0; var < NUM_VARS; var++)
        {



        }
        fprintf(file, "%f\n", cells[i].U[0][0]);


    }

    fclose(file);

    fileIndex++;
    //printf("done.\n");
    printf("done.\n");
}
















