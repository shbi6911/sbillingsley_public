/*This program uses functions named rotate and write_csv to take x,y coordinates of a ground station,
and model the station rotating with the earth throughout the course of a day.  It then outputs a
csv file of x,y,z coordinates for each minute of the day.*/
#include <iostream>
#include "utilities.h"

int main (){
    double cbtemp[] = {-4460.49,2682.08};   //initial coordinates
    double theta = ((360.0/86400.0) * 60);  //angle of rotation per minute
    double cb[4320];                        //preallocating results array
    cb[0] = cbtemp[0];                      //assigning initial values
    cb[1] = cbtemp[1];
    cb[2] = -3674.26;
    for (int i=0; i < 4320; i+=3){         //iterate over every triplet of coords in 1D array
        rotate(cbtemp,2,theta);           //rotate coords, and store in cb
        cb[i+3] = cbtemp[0];
        cb[i+4] = cbtemp[1];
        cb[i+5] = -3674.26;
    }
    write_csv(cb,4320,"CBPosition.csv");  //writing csv file from cb array
    
return 0;    
}