/*ASEN 1400 Gateway to Space
 * This code runs on the arduino nano to gather altitude data and pop a balloon for return to ground
 * Author: Matt Rhode
 */

/////////////////////////////////---------- Declarations ----------///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include <Adafruit_Sensor.h>                                  //includes Adafruit Sensor library
#include <Adafruit_BMP280.h>                                  //includes Adafruit BME280 pressure sensor library
#include <Wire.h>                                             //Includes Wire library, necessary for I2C communication
#define SEALEVELPRESSURE_HPA (1013.25)                        //Defines Sea Level Pressure for the BME280 sensor

Adafruit_BMP280 bmp;                                          // bme is identifying the sensor and preceeds functions in the bme280 library

float init_height;                                            //creates decimal-enabled variable for the initial height when powered
float max_height;                                             //creates decimal-enabled variable for the height at which the user would like to pop the balloon
float current_height;                                         //creates a decimal-enabled variable for measuring the instantaneous height of the balloon through flight

/////////////////////////////////---------- Setup Block ----------/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void setup()                                                  //setup block

{
  Serial.begin(9600);                                         //Begin communications with the serial monitor for data acquisition
  unsigned status;                                            //create an unsigned variable (non-negative max 32 bits) called "status"
  status = bmp.begin(0x76);                                   //assigns the variable "status" the address for the BME280 pressure sensor
  pinMode (4, OUTPUT);                                        //Sets digital pin 4 to output mode for activating burn wire  
  init_height = (bmp.readAltitude(SEALEVELPRESSURE_HPA));     //reads the pressure sensor at power-up and sets the height at the ground
}

/////////////////////////////////---------- Loop Block ----------//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void loop()                                                  //loop block (repeats)

{ 
  digitalWrite (4, LOW);                                      //force digital pin 4 logic to low (no power)
    current_height =(bmp.readAltitude(SEALEVELPRESSURE_HPA)); //reads the current height of the balloon (meters) and set's it to variable current height
    max_height = 15;                                           //set height (meters) at which you'd like to pop the balloon
      
    if (current_height-init_height > max_height) {            //IF statement checks to see if the current height
          digitalWrite (4, HIGH);                             //is higher than desired max height.  If so, it 
          delay(1000);                                        //activates the burn wire and waits 1 second
          digitalWrite (4, LOW);}                             //stops power to burn wire
          
    //delay(1000);                                              // for testing only, comment out for flight
    Serial.print(millis());                                   //Writes out a time stamp to the serial monitor
    Serial.print(";");                                        //Writes out semi-colon to serial monitor so data can be delimited
    Serial.println(current_height-init_height);               //Writes out altitude reading to serial monitor, uses read pressure, compares to sea level (defined above), then returns altitude
}
