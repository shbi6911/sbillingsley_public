//////////////////////////////////////////////////
// Radiation-Watch.org
// URL http://www.radiation-watch.org/
//////////////////////////////////////////////////

///　Digital I/O PIN Settings　///
int signPin1 = 2; //Radiation Pulse (Yellow)
int signPin2 = 9;

int noisePin1 = 5; //Vibration Noise Pulse (White)
int noisePin2 = 11;
//VCC 5V (Red)
//GND (Blue)
////////////////////////////////

const double alpha = 53.032; // cpm = uSv x alpha

int index = 0; //Number of loops
char msg[256] = ""; //Message buffer for serial output
char debug1[128] = "";
char debug2[128] = "";

int signCount1 = 0; //Counter for Radiation Pulse
int signCount2 = 0;
int noiseCount1 = 0; //Counter for Noise Pulse
int noiseCount2 = 0;

int sON1 = 0; //Lock flag for Radiation Pulse
int sON2 = 0;
int nON1 = 0; //Lock flag for Noise Pulse
int nON2 = 0;

double cpm1 = 0; //Count rate [cpm] of current
double cpmHistory1[100]; //History of count rates
int cpmIndex1 = 0; //Position of current count rate on cpmHistory[]
int cpmIndexPrev1 = 0; //Flag to prevent duplicative counting

double cpm2 = 0; //Count rate [cpm] of current
double cpmHistory2[100]; //History of count rates
int cpmIndex2 = 0; //Position of current count rate on cpmHistory[]
int cpmIndexPrev2 = 0; //Flag to prevent duplicative counting

//Timing Settings for Loop Interval
int prevTime = 0;
int currTime = 0;

int totalSec = 0; //Elapsed time of measurement [sec]
int totalHour = 0; //Elapsed time of measurement [hour]

//Time settings for CPM calcuaration
int cpmTimeMSec1 = 0;
int cpmTimeSec1 = 0;
int cpmTimeMin1 = 0;

int cpmTimeMSec2 = 0;
int cpmTimeSec2 = 0;
int cpmTimeMin2 = 0;

//String buffers of float values for serial output
char cpmBuff[20];
char uSvBuff[20];
char uSvdBuff[20];

void setup()
{
  //Serial setup
  //9600bps
  Serial.begin(9600);

  digitalWrite(13, HIGH); //turn on LED

  //PIN setting for Radiation Pulse
  pinMode(signPin1, INPUT);
  digitalWrite(signPin1, HIGH);
  pinMode(signPin2, INPUT);
  digitalWrite(signPin2, HIGH);

  //PIN setting for Noise Pulse
  pinMode(noisePin1, INPUT);
  digitalWrite(noisePin1, HIGH);
  pinMode(noisePin2, INPUT);
  digitalWrite(noisePin2, HIGH);

  //CSV-formatting for serial output (substitute , for _)
  Serial.println("ctr_hour[h]_sec[s]_count_cpm_uSv/h_uSv/hError");

  //Initialize cpmHistory[]
  for (int i = 0; i < 100; i++ )
  {
    cpmHistory1[i] = 0;
    cpmHistory2[i] = 0;
  }

  //Get start time of a loop
  prevTime = millis();
}

void loop()
{
  // Raw data of Radiation Pulse: Not-detected -> High, Detected -> Low
  int sign1 = digitalRead(signPin1);  //(GEIGER 1)
  int sign2 = digitalRead(signPin2);  //(GEIGER 2)

  // Raw data of Noise Pulse: Not-detected -> Low, Detected -> High
  int noise1 = digitalRead(noisePin1);  //(GEIGER 1)
  int noise2 = digitalRead(noisePin2);  //(GEIGER 2)

  //Radiation Pulse normally keeps low for about 100[usec]  //(GEIGER 1)
  if (sign1 == 0 && sON1 == 0)
  { //Deactivate Radiation Pulse counting for a while
    sON1 = 1;
    signCount1++;
  } else if (sign1 == 1 && sON1 == 1) {
    sON1 = 0;
  }
  if (sign2 == 0 && sON2 == 0)
  { //Deactivate Radiation Pulse counting for a while
    sON2 = 1;
    signCount2++;
  } else if (sign2 == 1 && sON2 == 1) {
    sON2 = 0;
  }

  //Noise Pulse normally keeps high for about 100[usec]   //(GEIGER 2)
  if (noise1 == 1 && nON1 == 0)
  { //Deactivate Noise Pulse counting for a while
    nON1 = 1;
    noiseCount1++;
  } else if (noise1 == 0 && nON1 == 1) {
    nON1 = 0;
  }
  if (noise2 == 1 && nON2 == 0)
  { //Deactivate Noise Pulse counting for a while
    nON2 = 1;
    noiseCount2++;
  } else if (noise2 == 0 && nON2 == 1) {
    nON2 = 0;
  }

  //Output readings to serial port, after 10000 loops
  if (index == 10000) //About 160-170 msec in Arduino Nano(ATmega328)
  {
    //Get current time
    currTime = millis();

    //No noise detected in 10000 loops (GEIGER 1)
    if (noiseCount1 == 0)
    {
      //Shift an array for counting log for each 6 sec.
      if ( totalSec % 6 == 0 && cpmIndexPrev1 != totalSec)
      {
        cpmIndexPrev1 = totalSec;
        cpmIndex1++;

        if (cpmIndex1 >= 100)
        {
          cpmIndex1 = 0;
        }

        if (cpmHistory1[cpmIndex1] > 0)
        {
          cpm1 -= cpmHistory1[cpmIndex1];
        }
        cpmHistory1[cpmIndex1] = 0;
      }

      //Store count log
      cpmHistory1[cpmIndex1] += signCount1;
      //Add number of counts
      cpm1 += signCount1;

      //Get ready time for 10000 loops
      cpmTimeMSec1 += abs(currTime - prevTime);
      //Transform from msec. to sec. (to prevent overflow)
      if (cpmTimeMSec1 >= 1000)
      {
        cpmTimeMSec1 -= 1000;
        //Add measurement time to calcurate cpm readings (max=20min.)
        if ( cpmTimeSec1 >= 10 * 60 )
        {
          cpmTimeSec1 = 10 * 60;
        } else {
          cpmTimeSec1++;
        }

        //Total measurement time
        totalSec++;
        //Transform from sec. to hour. (to prevent overflow)
        if (totalSec >= 3600)
        {
          totalSec -= 3600;
          totalHour++;
        }
      }

      //Elapsed time of measurement (max=20min.)
      double min = cpmTimeSec1 / 60.0;
      if (min != 0)
      {
        //Calculate cpm, uSv/h and error of uSv/h
        dtostrf(cpm1 / min, -1, 3, cpmBuff);
        dtostrf(cpm1 / min / alpha, -1, 3, uSvBuff);
        dtostrf(sqrt(cpm1) / min / alpha, -1, 3, uSvdBuff);
      } else {
        //Devision by zero
        dtostrf(0, -1, 3, cpmBuff);
        dtostrf(0, -1, 3, uSvBuff);
        dtostrf(0, -1, 3, uSvdBuff);
      }

      //Create message for serial port
      sprintf(msg, "1,%d,%d.%03d,%d,%s,%s,%s",
              totalHour, totalSec,
              cpmTimeMSec1,
              signCount1,
              cpmBuff,
              uSvBuff,
              uSvdBuff
             );
      sprintf (debug1, "sign1, %d, noise1, %d",signCount1, noiseCount1);
      //Send message to serial port
      Serial.println(msg);
      Serial.println(debug1);
    }

    //No noise detected in 10000 loops (GEIGER 2)
    if (noiseCount2 == 0)
    {
      //Shift an array for counting log for each 6 sec.
      if ( totalSec % 6 == 0 && cpmIndexPrev2 != totalSec)
      {
        cpmIndexPrev2 = totalSec;
        cpmIndex2++;

        if (cpmIndex2 >= 100)
        {
          cpmIndex2 = 0;
        }

        if (cpmHistory2[cpmIndex2] > 0)
        {
          cpm2 -= cpmHistory2[cpmIndex2];
        }
        cpmHistory2[cpmIndex2] = 0;
      }

      //Store count log
      cpmHistory2[cpmIndex2] += signCount2;
      //Add number of counts
      cpm2 += signCount2;

      //Get ready time for 10000 loops
      cpmTimeMSec2 += abs(currTime - prevTime);
      //Transform from msec. to sec. (to prevent overflow)
      if (cpmTimeMSec2 >= 1000)
      {
        cpmTimeMSec2 -= 1000;
        //Add measurement time to calcurate cpm readings (max=20min.)
        if ( cpmTimeSec2 >= 10 * 60 )
        {
          cpmTimeSec2 = 10 * 60;
        } else {
          cpmTimeSec2++;
        }
      }

      //Elapsed time of measurement (max=20min.)
      double min = cpmTimeSec2 / 60.0;
      if (min != 0)
      {
        //Calculate cpm, uSv/h and error of uSv/h
        dtostrf(cpm2 / min, -1, 3, cpmBuff);
        dtostrf(cpm2 / min / alpha, -1, 3, uSvBuff);
        dtostrf(sqrt(cpm2) / min / alpha, -1, 3, uSvdBuff);
      } else {
        //Devision by zero
        dtostrf(0, -1, 3, cpmBuff);
        dtostrf(0, -1, 3, uSvBuff);
        dtostrf(0, -1, 3, uSvdBuff);
      }

      //Create message for serial port
      sprintf(msg, "2,%d,%d.%03d,%d,%s,%s,%s",
              totalHour, totalSec,
              cpmTimeMSec2,
              signCount2,
              cpmBuff,
              uSvBuff,
              uSvdBuff
             );
      sprintf (debug2, "sign2, %d, noise2, %d",signCount2, noiseCount2);
      //Send message to serial port
      Serial.println(msg);
      Serial.println(debug2);
    }

    //Initialization for next 10000 loops
    prevTime = currTime;
    signCount1 = 0;
    signCount2 = 0;
    noiseCount1 = 0;
    noiseCount2 = 0;
    index = 0;
  }
  index++;
}
