//These are member functions for class Coord, defined in file coord.h
#include "coord.h"

/*setGroundPoint is a mutator function
    Input: x, y and z values which are set to Coord values ground_xpoint, ground_ypoint and ground_zpoint
        this represents the position of a ground station at a given time point
    Output: none*/
void Coord::setGroundPoint (double x,double y, double z){
    ground_xpoint = x;
    ground_ypoint = y;
    ground_zpoint = z;
    
}
/*setSpacePoint is a mutator function
    Input: x, y and z values which are set to Coord values space_xpoint, space_ypoint and space_zpoint
        this represents the position of a satellite at a given time point
    Output: none*/
void Coord::setSpacePoint (double x,double y, double z){
    space_xpoint = x;
    space_ypoint = y;
    space_zpoint = z;
    
}
/*setFlag calculates the masking angle phi for a given ground station and satellite position
        it then stores the masking angle phi in degrees and sets visible_flag to TRUE if phi is greater than 10 degrees
    Input: two sets of x, y, and z values, one for ground, one for satellite
    Output: none (visible_flag set to TRUE or FALSE)*/
void Coord::setFlag (){
    double dx = (space_xpoint - ground_xpoint);
    double dy = (space_ypoint - ground_ypoint);
    double dz = (space_zpoint - ground_zpoint);
    double dotprod = ((dx*ground_xpoint) + (dy*ground_ypoint) + (dz*ground_zpoint));
    double gnorm = sqrt((ground_xpoint*ground_xpoint) + (ground_ypoint*ground_ypoint) + (ground_zpoint*ground_zpoint));
    double dnorm = sqrt((dx*dx) + (dy*dy) + (dz*dz));
    double phirad = ((M_PI/2) - acos(dotprod/(gnorm*dnorm)));
    phi = phirad * (180/M_PI);
    visible_flag = (phi > 10);
}
//displayInfo outputs selected values of class Coord to the screen
void Coord::displayInfo(){
    cout << "Ground Station: " << endl;
    cout << "x: " << ground_xpoint << endl;
    cout << "y: " << ground_ypoint << endl;
    cout << "z: " << ground_zpoint << endl;
    cout << "Satellite: " << endl;
    cout << "x: " << space_xpoint << endl;
    cout << "y: " << space_ypoint << endl;
    cout << "z: " << space_zpoint << endl;
    cout << "Masking Angle " << phi << " degrees" << endl;
    cout << "Visibility: " << visible_flag << endl;
}
//export_flag is a simple get function that returns the value of visible_flag
    //Input: none
    //Output: boolean value of visible_flag
bool Coord::export_flag(){
    return visible_flag;
}