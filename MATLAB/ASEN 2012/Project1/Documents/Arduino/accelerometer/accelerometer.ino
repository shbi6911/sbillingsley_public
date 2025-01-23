void setup() {
  // put your setup code here, to run once:
int sensorx;
int sensorz;
float sensorVoltx;
float sensorVoltz;
float sensorUnitx;
float sensorUnitz;
Serial.begin(9600);
}

void loop() {
  // put your main code here, to run repeatedly:
sensorx = analogRead(A4);
sensorz = analogRead(A5);
sensorVoltx = sensorx*(5/1023);
sensorVoltz = sensorz*(5/1023);
sensorUnitx = (sensorVoltx - (3.3/2))/(0.330);
sensorUnitz = (sensorVoltz - (3.3/2))/(0.330);
Serial.print("X ");
Serial.print(sensorUnitx);
Serial.print("Z ")
Serial.print(sensorUnitz);

}
