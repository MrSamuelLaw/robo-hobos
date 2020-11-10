#include "./readPinFast.h"

// Interrupt pins in Mega 2560 are : 2 (interrupt 0), 3 (interrupt 1), 18 (interrupt 5),
// 19 (interrupt 4), 20 (interrupt 3), and 21 (interrupt 2)

// The #define construct is a compiler pre-processor directive
// The value is substituted for the name, wherever the name occurs in the code.
#define enc1Interrupt 2    //Interrupt pin number
#define enc1PinA      21   //Arduino pin number corresponding to interrupt pin
#define enc1PinB      31   //Digital pin number // Do not use pins 0, 1, 3, 8, 9, 11, 12, 13
#define MotorRpwm     3
#define MotorRbrk     9
#define MotorRdir     12

// A variable should be declared volatile whenever its value may change
// at any time-without any action being taken by the code found nearby.
volatile long enc1Ticks = 0;  // Variable used to store encoder ticks

void setup() {
  Serial.begin(9600); // Sets baud rate for data transmission to the PC.
  pinMode(MotorRpwm,OUTPUT);
  pinMode(MotorRbrk,OUTPUT);
  pinMode(MotorRdir,OUTPUT);

  // Set up Encoder-1
  pinMode(enc1PinA, INPUT);      // sets Encoder-1 pin A as input
  pinMode(enc1PinB, INPUT);      // sets Encoder-1 pin B as input
  attachInterrupt(enc1Interrupt, enc1ReadA, RISING); // Executes the function 'enc1ReadA()' at rising edge of signal A from Encoder-1
}

void loop() {
  // Change the calculation in following line to convert encoder ticks into angle in degrees or radians
  double enc1Angle = enc1Ticks / 5.688;

  // Print encoder positions
  Serial.print("Enc-1: "); Serial.println(enc1Angle);

  digitalWrite(12,HIGH);
  digitalWrite(9,LOW);
  analogWrite(MotorRpwm,150);
}

void enc1ReadA()
{
  // Write your code here.
  // Remember, you are supposed to read port-B value using 'readPinFast(pinNumber)' function.
  // If this value is 0, then add 1 to 'enc1Ticks' else subtract 1 from 'enc1Ticks'.



  if (readPinFast(enc1PinB) == 0) {
    --enc1Ticks;
  }
  else{
    ++enc1Ticks;
  }
}
