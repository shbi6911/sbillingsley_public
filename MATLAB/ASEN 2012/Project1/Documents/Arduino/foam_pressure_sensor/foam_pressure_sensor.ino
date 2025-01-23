float volt;                                 //define variable of type float (allows decimals) for read voltage                                 

void setup() {                              // put your setup code here, to run once:
    Serial.begin(9600);                     //start serial communications
}

void loop() {                               // put your main code here, to run repeatedly:
    volt = analogRead(A5)*5.0/1024.0-1;     // reads in voltage from pin A5, digitizes, then convert back to voltage
    Serial.println(volt);                   // print voltage value to serial monitor/plotter           
    delay(100);                             // delay so you can read values
}
