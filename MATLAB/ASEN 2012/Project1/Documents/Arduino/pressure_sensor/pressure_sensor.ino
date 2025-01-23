void setup() {
  // put your setup code here, to run once:
Serial.begin(9600);
int sensor;
float sensorvolt;
float sensorUnits;


}

void loop() {
  // put your main code here, to run repeatedly:
sensor = analogRead(A3);
sensorvolt = sensor *(5/1023);
sensorUnits = (sensorvolt-0.5)*(15.0/4.0);
Serial.print(sensor);
Serial.print("\t voltage");
Serial.print(sensorvolt);
Serial.print("\t psi");
Serial.print(sensorUnits);
}
