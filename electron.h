//
//  electron.h
//  HARM2D
//
//  Created by Sean Ressler on 3/5/15.
//  Copyright (c) 2015 Home. All rights reserved.
//

/*
Conduction(): Handles electron heat conduction
eVol(): Updates electron internal energy due to volume changes
felcalc(): Calculates electron fraction
fixupuel(): Applies floors to electron internal energy across the grid
fixupuel1zone(): Applies floors to electron internal energy in a single zone
Heating(): Updates electron entropy due to heating terms
init_entropy(): Initializes electron entropy based on initial conditions
*/
void Conduction(double pi[][N2M][N3M][NPR], double prh[][N2M][N3M][NPR], double pr[][N2M][N3M][NPR], int i, int j, int k, double Dt, int was_floor_activated);
void eVol(double pi[][N2M][N3M][NPR],double prh[][N2M][N3M][NPR], double pr[][N2M][N3M][NPR], double Dt, int i, int j, int k, int was_floor_activated,int is_after_fixup);
double felcalc(double pi[][N2M][N3M][NPR], int i, int j, int k, double Tel);
void fixupuel(double (* pv)[N2M][N3M][NPR]);
void fixupuel1zone( int i, int j, int k, double prim[NPR] ) ;
void Heating(double pi[][N2M][N3M][NPR], double prh[][N2M][N3M][NPR], double pr[][N2M][N3M][NPR], int i, int j, int k, double Dt, int was_floor_activated, int is_after_fixup);
void init_entropy(void);




