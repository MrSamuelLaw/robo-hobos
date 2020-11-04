#include <Arduino.h>
#include <Encoder.h>

// ============ define constants =============
double  ki = (10.52 - 0.108605472*4.5)/45;
double kt = (0.980665/2.6666666660);
float R = 4.5;
int kp = 5;
int kv = 1;

// ============= define link properties =============
float bodyA_mass{1};    // kg
float bodyA_izz{1};     // kg*m^2
float bodyA_length{1};  // m

// ============= define right motor variables =============
Encoder right_encoder(20, 21);
long right_old_pos{0};
long right_new_pos{0};
double right_ang{0};
double right_omega{0};

constexpr long CPR_right{1000};
int motRpwm{3};
int motRbrk{9};
int motRdir{12};

// define variables for torque calculations
double right_desired{3};  // radians
double right_Vemf{0};
double right_Vs{0};
double right_Tau{0};
double right_U{0};

// ============== setup clocks ===============
unsigned long old_time{0};
unsigned long new_time{0};
unsigned long dt{0};

// ============== utility funcs ==============
double ticks_to_rad(int ticks){
  return (double)(-2*M_PI*ticks/CPR_right);
}

// ================== setup ==================
void setup() {
  // start the coms
  Serial.begin(9600);

  // define outputs
  pinMode(motRpwm, OUTPUT);
  pinMode(motRbrk, OUTPUT);
  pinMode(motRdir, OUTPUT);

  // establish the connection
  digitalWrite(motRdir, HIGH);
  digitalWrite(motRbrk, LOW);
  analogWrite(motRpwm, 50);

}

// ================== loop ===================
void loop() {
  // read the encoder values
  right_new_pos = right_encoder.read();

  // check if any new values
  if (right_new_pos != right_old_pos) {
    // reset the clock
    new_time = millis();

    // calculate dt
    dt = new_time - old_time;

    // calculate angles in radians
    float rad_new = ticks_to_rad(right_new_pos);
    float rad_old = ticks_to_rad(right_old_pos);

    // calculate omega
    right_omega = (rad_new - rad_old)/dt;

    // calculate back voltage
    right_Vemf = ki*right_omega;

    // calculate controller output
    right_U = kv*(-right_omega) + kv*(right_desired - rad_new);

    // calculate the torqu3
    right_Tau = ((bodyA_length*bodyA_length*bodyA_mass) + bodyA_izz)*right_U;


    // calculate the desired voltage assuming Tau known
    right_Vs = (right_Tau*R/kt) + right_Vemf;

    // constrain the voltage
    right_Vs = constrain(right_Vs, 0, 12);
    right_Vs = map(right_Vs, 0, 12, 0, 255);

    // output it to the motor
    analogWrite(motRpwm, right_Vs);


    // update the values
    right_old_pos = right_new_pos;
    old_time = new_time;
  }
}

