#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

class Coord {
  public:
    void setGroundPoint(double, double, double );
    void setSpacePoint(double, double, double);
    void setFlag();
    void displayInfo();
    bool export_flag();
  private:
    double ground_xpoint;
    double ground_ypoint;
    double ground_zpoint;
    double space_xpoint;
    double space_ypoint;
    double space_zpoint;
    double phi;
    bool visible_flag;
};