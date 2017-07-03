 // temp logger based around an LM335 TO92 IC
// nb gain is 10mV/K
// -40C < Toperating < 100C
// More: bit.ly/1SAiQ9K

int Uf=0;  // voltage from LM335
int lastReading=0;  // last reading
long elapsed = 0;
int LED=13;
int LEDgreen=12;
int LEDorangeA=8;
int LEDorangeB=7;
int LEDred=4;


void setup()
{
  Serial.begin(9600);
  pinMode(LED,OUTPUT);
  pinMode(LEDgreen,OUTPUT);
  pinMode(LEDorangeA,OUTPUT);
  pinMode(LEDorangeB,OUTPUT);
  pinMode(LEDred,OUTPUT);
}

void loop()
{
  ledBlink();
  ledBlink();
  ledBlink();
  ledBlink();
  elapsed = millis();
  // Read sensor output
  Uf = analogRead(0);
  // Convert bitwise output to mV
  // (A0 width 1024 bits, range 5000mV)
  float mV = (Uf/1024.0)*5000.0;
  // Convert mV to Kelvin (sens is 10mV/K)
  float temp_K = mV / 10.0;
  // Convert K to C
  float temp_C = temp_K - 273.15;
  
  // if writing to a terminal and moving cursor back is needed
  //Serial.write(27);
  //Serial.print("[2J");
  //Serial.write(27);
  //Serial.print("[H");
  //Serial.write(13);
  //Serial.print(Uf);
  for(int i=0;i<10;i++){
    ledBlink();
    delay(10);
  }
  // if writing to a file and each value on a new line
  Serial.print(elapsed);
  Serial.print("\t");
  Serial.print(mV);
  Serial.print("\t");
  Serial.print(temp_K);
  Serial.print("\t");
  Serial.print(temp_C);
  Serial.print("\t");
  Serial.println(Uf);
  //delay(10000);
  // turn LEDs off and on depending on temp
  if(temp_C < -10.0){
    digitalWrite(LEDgreen,HIGH);
  }else{
    digitalWrite(LEDgreen,LOW);
  }
  
  if(Uf < lastReading){
    digitalWrite(LEDorangeA,HIGH);
  }else{
    digitalWrite(LEDorangeA,LOW);
  }

  if(Uf > lastReading){
    digitalWrite(LEDorangeB,HIGH);
  }else{
    digitalWrite(LEDorangeB,LOW);
  }

  if(temp_C > 6.0){
    digitalWrite(LEDred,HIGH);
  }else{
    digitalWrite(LEDred,LOW);
  }
  lastReading = Uf;
}

void ledBlink(){
  digitalWrite(LED,HIGH);
  delay(50);
  digitalWrite(LED,LOW);
  delay(50);
}


