#include "coord.h"

void read_csv (ifstream&, double&, double&, double&);

int main (){
    ifstream sat1Position;                      //open files for input
    sat1Position.open ("Sat1Position.csv");
    ifstream sat2Position;
    sat2Position.open ("Sat2Position.csv");
    ifstream CBPosition;
    CBPosition.open ("CBPosition.csv");
    
    double s1x, s1y, s1z;
    //sat1Position >> s1x >> s1y >> s1z;
    //cout << "x: " << s1x << " y: " << s1y << " z: " << s1z << endl;
    
    //string line1, line2, line3;
    double x,y,z;
    
    read_csv(sat1Position, x, y, z);
    cout << x << " , " << y << " , " << z << endl;
    read_csv(sat2Position, x, y, z);
    cout << x << " , " << y << " , " << z << endl;
    read_csv(CBPosition, x, y, z);
    cout << x << " , " << y << " , " << z << endl;
    
    // getline(sat1Position, line1, ',');
    // getline(sat1Position, line2, ',');
    // getline(sat1Position, line3);
    // x = stod (line1);
    // y = stod (line2);
    // z = stod (line3);
    // cout << line1 << " , " << line2 << " , " << line3 << endl;
    // cout << x << " , " << y << " , " << z << endl;
    
    // getline(sat1Position, line1, ',');
    // getline(sat1Position, line2, ',');
    // getline(sat1Position, line3);
    // x = stod (line1);
    // y = stod (line2);
    // z = stod (line3);
    // cout << line1 << " , " << line2 << " , " << line3 << endl;
    // cout << x << " , " << y << " , " << z << endl;
    
    sat1Position.close();                       //close files
    sat2Position.close();
    CBPosition.close();
    
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