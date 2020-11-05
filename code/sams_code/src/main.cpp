#include <Arduino.h>
#include <Encoder.h>


//   __  __       _                _____ _
//  |  \/  |     | |              / ____| |
//  | \  / | ___ | |_ ___  _ __  | |    | | __ _ ___ ___
//  | |\/| |/ _ \| __/ _ \| '__| | |    | |/ _` / __/ __|
//  | |  | | (_) | || (_) | |    | |____| | (_| \__ \__ \
//  |_|  |_|\___/ \__\___/|_|     \_____|_|\__,_|___/___/
class Motor{
  public:
  /*
  This is a class meant to implement a PD or operational space controller for a motor
  connected to a robot arm.
  1. For the initial setup, set the encoder pins such that as the
     motor spindle turns counter clockwise, the encoder reads a positive angle.

  2. Once the encoder is reading a positive angle for a counter clockwise rotation,
     run a test with the COUNTER_CLOCKWISE and CLOCKWISE values set as the are.
     If the motor starts to spin the wrong way, swap the COUNTER_CLOCKWISE and
     CLOCKWISE assignments.
  */

  // ============================= constructor =============================

  Motor(int pwm_pin, int brk_pin, int dir_pin, int enc_pin_A, int enc_pin_B)
    :m_encoder{Encoder(enc_pin_A, enc_pin_B)}
  {
    // set motor pins numbers
    m_pwm_pin = pwm_pin;
    m_brk_pin = brk_pin;
    m_dir_pin = dir_pin;

    // set the pin modes
    pinMode(m_pwm_pin, OUTPUT);
    pinMode(m_brk_pin, OUTPUT);
    pinMode(m_dir_pin, OUTPUT);

    // set pin values
    analogWrite(m_pwm_pin, 0);
    digitalWrite(m_brk_pin, LOW);
    // if it turns clockwise by mistake, flip the HIGH and LOW assignments
    COUNTER_CLOCKWISE = HIGH;
    CLOCKWISE = LOW;
    digitalWrite(m_dir_pin, COUNTER_CLOCKWISE);

  }

  // ============================== variables ==============================
  // pin variables
  int m_pwm_pin{0};
  int m_brk_pin{0};
  int m_dir_pin{0};
  int m_dir_state_pin{0};
  uint8_t CLOCKWISE{LOW};
  uint8_t COUNTER_CLOCKWISE{HIGH};

  // encoder variables
  Encoder m_encoder;

  // controller variables
  double q[2] = {0, 0};  // rad
  double q_desired{0};   // rad
  double qdot = {0};     // rad/s
  double U{0};           // units/s^2
  double Tau{0};         // Nm
  double Vemf{0};        // Volts
  double Vs{0};          // Volts
  long kp{80000};                    // proportional gain
  long kv{sqrt(kp)*2.828427};        // derivative gain

  // controller constants
  float bodyA_mass{0.03190818};    // kg
  float bodyA_izz{3.98e-06};     // kg*m^2
  float bodyA_length{0.02162775};  // m

  // motor/encoder constants
  int CPR{1000};
  double  ki = (10.52 - 0.108605472*4.5)/45;
  double kt = (0.980665/2.6666666660);
  float R = 4.5;

  // ================================ functions ================================

  // 1 for counter clockwise, -1 for clockwise
  // same as the convection for positive angle rotations
  void set_direction(int dir){

    if (dir == 1){
      digitalWrite(m_dir_pin, COUNTER_CLOCKWISE);
    }
    else if (dir == -1){
      digitalWrite(m_dir_pin, CLOCKWISE);
    }
  }

  // returns the current angle of the motor in radians
  double angle(){
    // returns the current angle in radians
    return (double)(2*M_PI*m_encoder.read()/CPR);
  }


  // updates the angles to where q1 is always the most recent angle
  void update_angles(){
    // q[1] is the newest
    q[0] = q[1];
    q[1] = angle();
  }


  // implements the controller
  double calculate_Vs(long millis){
    // get the newest angles
    update_angles();
    Serial.print("q1 = ");
    Serial.println(q[1]);

    // calculate qdot
    qdot = (double)(q[1] - q[0])*(1000/millis);
    // Serial.println(qdot);

    // calculate back emf
    Vemf = ki*qdot;
    // Serial.println(Vemf);

    // calculate controller value for U
    U = (kv*(0-qdot)) + (kp*(q_desired - q[1]));
    // Serial.println(U);

    // calculate Tau from controller value
    Tau = ((bodyA_length*bodyA_length*bodyA_mass) + bodyA_izz)*U;
    // Serial.println(Tau);

    // calculate Vs
    Vs = (Tau*R/kt) + Vemf;
    // Serial.println(Vs);

    // toggle the direction if needed
    set_direction(1);
    if (Vs < 0){
      Vs = -Vs;
      set_direction(-1);
    }

    // constrain the voltage
    Vs = constrain(Vs, 0, 12);
    // Serial.println(Vs);

    // map the voltage
    Vs = map(Vs, 0, 12, 0, 255);
    // Serial.println(Vs);
    return Vs;
  }

};


//    _____ _       _           _
//   / ____| |     | |         | |
//  | |  __| | ___ | |__   __ _| |___
//  | | |_ | |/ _ \| '_ \ / _` | / __|
//  | |__| | | (_) | |_) | (_| | \__ \
//   \_____|_|\___/|_.__/ \__,_|_|___/

// controller variables


// time variables
long new_time{0};
long old_time{0};
long dt{0};
float V{0};

// define motors
Motor right_motor(3, 9, 12, 20, 21);

//                    _       _                _____      _
//      /\           | |     (_)              / ____|    | |
//     /  \   _ __ __| |_   _ _ _ __   ___   | (___   ___| |_ _   _ _ __
//    / /\ \ | '__/ _` | | | | | '_ \ / _ \   \___ \ / _ \ __| | | | '_ \
//   / ____ \| | | (_| | |_| | | | | | (_) |  ____) |  __/ |_| |_| | |_) |
//  /_/    \_\_|  \__,_|\__,_|_|_| |_|\___/  |_____/ \___|\__|\__,_| .__/
//                                                                 | |
//                                                                 |_|

void setup(){
  Serial.begin(9600);
  // set desired angle
  right_motor.q_desired = 12.5;
}

//                    _       _               _
//      /\           | |     (_)             | |
//     /  \   _ __ __| |_   _ _ _ __   ___   | |     ___   ___  _ __
//    / /\ \ | '__/ _` | | | | | '_ \ / _ \  | |    / _ \ / _ \| '_ \
//   / ____ \| | | (_| | |_| | | | | | (_) | | |___| (_) | (_) | |_) |
//  /_/    \_\_|  \__,_|\__,_|_|_| |_|\___/  |______\___/ \___/| .__/
//                                                             | |
//                                                             |_|

void loop(){
  // capture the current angle and time
  new_time = millis();
  dt = (new_time - old_time);

  // calculate Vs
  V = right_motor.calculate_Vs(dt);

  // output it to the motor
  analogWrite(right_motor.m_pwm_pin, V);

  // overwrite the old time
  old_time = new_time;
}


















// // ============ define constants =============
// double  ki = (10.52 - 0.108605472*4.5)/45;
// double kt = (0.980665/2.6666666660);
// float R = 4.5;
// long kv = 25000;
// long kp = sqrt(kv);

// // ============= define link properties =============
// float bodyA_mass{0.03190818};    // kg
// float bodyA_izz{3.98e-06};     // kg*m^2
// float bodyA_length{0.02162775};  // m

// // ============= define right motor variables =============
// Encoder encoder_right(20, 21);
// long old_pos_right{0};
// long new_pos_right{0};
// double ang_right{0};
// double omega_right{0};

// constexpr long CPR_right{1000};
// int motor_pwm_right{3};
// int motor_brk_right{9};
// int motor_dir_right{12};
// bool motor_dir_bool_right{false};

// // define variables for torque calculations
// double desired_right{3};  // radians
// double Vemf_right{0};
// double Vs_right{0};
// double Tau_right{0};
// double U_right{0};

// // ============== setup clocks ===============
// unsigned long old_time{0};
// unsigned long new_time{0};
// unsigned long dt{0};

// // ============== utility funcs ==============
// double ticks_to_rad(int ticks){
//   // converts ticks to radians
//   return (double)(-2*M_PI*ticks/CPR_right);
// }

// void toggle_direction_right(){
//   // flips the direction of the right motor
//   if (motor_dir_bool_right) {
//     motor_dir_bool_right = false;
//     digitalWrite(motor_dir_right, HIGH);
//   }
//   else {
//     motor_dir_bool_right = true;
//     digitalWrite(motor_dir_right, LOW);
//   }

// }

// void toggle_direction_left(){
//   // flips the direction of the left motor

// }



// // ================== setup ==================
// void setup() {
//   // start the coms
//   Serial.begin(9600);

//   // define outputs
//   pinMode(motor_pwm_right, OUTPUT);
//   pinMode(motor_brk_right, OUTPUT);
//   pinMode(motor_dir_right, OUTPUT);

//   // establish the connection
//   digitalWrite(motor_dir_right, LOW);
//   digitalWrite(motor_brk_right, LOW);
//   analogWrite(motor_pwm_right, 50);

// }

// // ================== loop ===================
// void loop() {
//   // read the encoder values
//   new_pos_right = encoder_right.read();

//   // check if any new values
//   if (new_pos_right != old_pos_right) {
//     // reset the clock
//     new_time = millis();

//     // calculate dt
//     dt = new_time - old_time;

//     // calculate angles in radians
//     float rad_new = ticks_to_rad(new_pos_right);
//     float rad_old = ticks_to_rad(old_pos_right);
//     Serial.print("q1 = ");
//     Serial.println(rad_new);

//     // calculate omega
//     omega_right = (rad_new - rad_old)/dt;

//     // calculate back voltage
//     Vemf_right = ki*omega_right;

//     // calculate controller output
//     U_right = kv*(-omega_right) + kv*(desired_right - rad_new);
//     Serial.print("U = ");
//     Serial.println(U_right);

//     // calculate the torqu3
//     Tau_right = ((bodyA_length*bodyA_length*bodyA_mass) + bodyA_izz)*U_right;
//     Serial.print("Tau = ");
//     Serial.println(Tau_right);

//     // calculate the desired voltage assuming Tau known
//     Vs_right = (Tau_right*R/kt) + Vemf_right;
//     Serial.print("Vs = ");
//     Serial.println(Vs_right);

//     // flip the direction if Vs is negative
//     if (Vs_right < 0) {
//       toggle_direction_right();
//       Vs_right = -1*Vs_right;
//     }

//     // constrain the voltage
//     Vs_right = constrain(Vs_right, 0, 12);
//     Vs_right = map(Vs_right, 0, 12, 0, 255);

//     // output it to the motor
//     analogWrite(motor_pwm_right, Vs_right);

//     // update the values
//     old_pos_right = new_pos_right;
//     old_time = new_time;
//   }
// }

