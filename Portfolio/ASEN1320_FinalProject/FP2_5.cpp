/*This main function will create two arrays of class Coord, defined in coord.h and coord.cpp.
It will fill the arrays with data from csv files Sat1Position, Sat2Position, and CBPostion,
and calculate masking angles and set visibility flags for these points.  It will output the
visibility flag data to two csv files Sat1Visibility and Sat2Visibility.*/
#include "coord.h"

//This function will read three coordinates from a csv file
void read_csv (ifstream&, double&, double&, double&);

int main (){
    
    Coord sat1_array[1441];                     //declare arrays of class Coord
    Coord sat2_array[1441];
    
    ifstream sat1Position;                      //open files for input
    sat1Position.open ("Sat1Position.csv");
    ifstream sat2Position;
    sat2Position.open ("Sat2Position.csv");
    ifstream CBPosition;
    CBPosition.open ("CBPosition.csv");
    
    double s1x, s1y, s1z; 
    double s2x, s2y, s2z;
    double gx, gy, gz;            //temp variables

    //this loop will define values of the two arrays of class Coord
    for (int i = 0; i < 1441; i++){
        read_csv(sat1Position, s1x, s1y, s1z);               //read from three csv files
        read_csv(sat2Position, s2x, s2y, s2z);
        read_csv(CBPosition, gx, gy, gz);
        
        sat1_array[i].setSpacePoint(s1x,s1y,s1z);           //set values for Sat1
        sat1_array[i].setGroundPoint(gx,gy,gz);
        sat1_array[i].setFlag();
        
        sat2_array[i].setSpacePoint(s2x,s2y,s2z);           //set values for Sat2
        sat2_array[i].setGroundPoint(gx,gy,gz);
        sat2_array[i].setFlag();
    }
    sat1Position.close();                            //close input files
    sat2Position.close();
    CBPosition.close();
    
    ofstream sat1Visibility;                        //open files for output
    sat1Visibility.open ("Sat1Visibility.csv");
    ofstream sat2Visibility;
    sat2Visibility.open ("Sat2Visibility.csv");
    
    bool temp1; bool temp2;                         //temp variables
    
    //this loop will output visible_flag data for each array to a csv file
    for (int i=0; i < 1441; i++){
        temp1 = sat1_array[i].export_flag();            //extract values
        temp2 = sat2_array[i].export_flag();
        sat1Visibility << temp1 << endl;                //write to csv
        sat2Visibility << temp2 << endl;
        //sat1_array[i].displayInfo();
        //sat2_array[i].displayInfo();
    }
    sat1Visibility.close();                         //close output files
    sat2Visibility.close();
return 0;
}



/*read_csv is a function to read x,y,z coordinate data from a csv file
    Inputs: an ifstream object and three double values, all passed by reference
    Outputs: none (three double values set to x,y,z coordinates from csv file)*/
void read_csv (ifstream &csv_file, double &x, double &y, double &z){
    string line1, line2, line3;
    getline(csv_file, line1, ',');
    getline(csv_file, line2, ',');
    getline(csv_file, line3);
    x = stod (line1);
    y = stod (line2);
    z = stod (line3);
    }