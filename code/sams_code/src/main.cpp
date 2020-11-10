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
  args:
    pwm_pin: Pin number that controls the voltage supplied to the motor.
    brk_pin: Pin number that cuts off the motor (LOW means the motor turns).
    dir_pin: Pin number that controls the direction of the pin.
    COUNTER_CLOCKWISE_: HIGH or LOW, this defines what the dir_pin is set to when
                        the motor is asked to turn counter clockwise.
    CLOCKWISE_: HIGH or LOW, this defines what the dir_pin is set to when the
                motor is asked to go clockwise.
    enc_pin_A: Pin that encoder channel A is plugged into (must be an interrupt pin)
    enc_pin_B: Pin that encoder challen B is plugged into (must be an interrupt pin)
    encoder_dir_: Set to either 1 or -1, used to correct the reading for the angle.
                  If a counter clockwise rotation causes the encoder to read negative,
                  flip the value from -1 to 1 or 1 to -1 depending on what it is
                  currently set to.
  */

  // ============================= constructor =============================

  Motor(int pwm_pin, int brk_pin, int dir_pin,
          uint8_t COUNTER_CLOCKWISE_, uint8_t CLOCKWISE_,
          int enc_pin_A, int enc_pin_B, int encoder_dir_)
    // set default values
    :m_encoder{Encoder(enc_pin_A, enc_pin_B)},  // instantiates the encoder
     encoder_dir{encoder_dir_},                // sets the encoder corrector
     m_pwm_pin{pwm_pin},                       // sets the motor pwm pin
     m_brk_pin{brk_pin},                       // sets the motor brk pin
     m_dir_pin{dir_pin},                       // sets the motor pwm pin
     COUNTER_CLOCKWISE{COUNTER_CLOCKWISE_},     // sets the CCW variable
     CLOCKWISE{CLOCKWISE_}                     // sets the CC variable
    // perform default function calls
  {
    // set the pin modes
    pinMode(m_pwm_pin, OUTPUT);
    pinMode(m_brk_pin, OUTPUT);
    pinMode(m_dir_pin, OUTPUT);

    // set pin values
    analogWrite(m_pwm_pin, 0);
    digitalWrite(m_brk_pin, LOW);
    // if it turns clockwise by mistake, flip the HIGH and LOW assignments
    digitalWrite(m_dir_pin, COUNTER_CLOCKWISE);
  }

  // ============================== variables ==============================
  // safe to playwith controller variables
  float bodyA_mass{0.03190818};      // kg
  float bodyA_izz{3.98e-06};         // kg*m^2
  float bodyA_length{0.02162775};    // m
  long kp{80000};                    // proportional gain
  long kv{sqrt(kp)*2.828427};        // derivative gain

  // pin variables
  int m_pwm_pin{0};
  const int m_brk_pin;
  int m_dir_pin{0};
  int m_dir_state_pin{0};
  uint8_t CLOCKWISE{LOW};
  uint8_t COUNTER_CLOCKWISE{HIGH};
  int encoder_dir{1};

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
    return (double)(encoder_dir*2*M_PI*m_encoder.read()/CPR);
  }

  // kills the motor in the event that the probe is tripped, preventint
  // damage to the robot. Is completely optional, and should be set in the
  // setup loop of the ardino code using a lambda
  static void probe(int probe_pin){
    digitalWrite(probe_pin, HIGH);
  }

  // updates the angles to where q1 is always the most recent angle
  void update_angles(){
    // q[1] is the newest
    q[0] = q[1];
    q[1] = angle();
  }

  // =============================== controller ================================

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

// time variables
long new_time{0};
long old_time{0};
long dt{0};
float V{0};

// define motors
Motor right_motor(3, 9, 12, HIGH, LOW,
                  20, 21, 1);

//                    _       _                _____      _
//      /\           | |     (_)              / ____|    | |
//     /  \   _ __ __| |_   _ _ _ __   ___   | (___   ___| |_ _   _ _ __
//    / /\ \ | '__/ _` | | | | | '_ \ / _ \   \___ \ / _ \ __| | | | '_ \
//   / ____ \| | | (_| | |_| | | | | | (_) |  ____) |  __/ |_| |_| | |_) |
//  /_/    \_\_|  \__,_|\__,_|_|_| |_|\___/  |_____/ \___|\__|\__,_| .__/
//                                                                 | |
//                                                                 |_|

void setup(){
  // start the coms
  Serial.begin(9600);

  // setup the probe that trips if the robot is going to break something
  pinMode(2, INPUT);
  attachInterrupt(
    digitalPinToInterrupt(2),
    [](){Motor::probe(right_motor.m_brk_pin);},
    HIGH
  );

  // set desired angle
  right_motor.q_desired = 12.5;

  // digitalWrite(right_motor.m_brk_pin, HIGH);
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