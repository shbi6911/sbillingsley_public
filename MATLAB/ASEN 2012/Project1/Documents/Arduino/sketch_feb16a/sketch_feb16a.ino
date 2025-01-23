void setup() {
  // put your setup code here, to run once:
Serial.begin(9600);
pinMode(5,OUTPUT);   //sets pin 5 to output mode
pinMode(6,OUTPUT);   //sets pin 6 to output mode
pinMode(7,OUTPUT);   //sets pin 7 to output mode
pinMode(9,OUTPUT);   //sets pin 9 to output mode
}

void loop() {
  // put your main code here, to run repeatedly:
Serial.println("hello");
digitalWrite(5,LOW);
digitalWrite(6,LOW);
digitalWrite(7,LOW);
digitalWrite(9,LOW);
delay(1000);              //delays one second (1000 ms)
digitalWrite(5,HIGH);    //HIGH voltage is 5V, LOW voltage is 0V
delay(500);
digitalWrite(6,HIGH);    //HIGH voltage is 5V, LOW voltage is 0V
delay(500);
digitalWrite(7,HIGH);    //HIGH voltage is 5V, LOW voltage is 0V
delay(500);
digitalWrite(9,HIGH);    //HIGH voltage is 5V, LOW voltage is 0V
delay(500);
}
