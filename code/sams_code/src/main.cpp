#include <Arduino.h>
#include <Encoder.h>

///////////////  __  __       _                _____ _                ///////////////
/////////////// |  \/  |     | |              / ____| |               ///////////////
/////////////// | \  / | ___ | |_ ___  _ __  | |    | | __ _ ___ ___  ///////////////
/////////////// | |\/| |/ _ \| __/ _ \| '__| | |    | |/ _` / __/ __| ///////////////
/////////////// | |  | | (_) | || (_) | |    | |____| | (_| \__ \__ \ ///////////////
/////////////// |_|  |_|\___/ \__\___/|_|     \_____|_|\__,_|___/___/ ///////////////
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
                  currently set to. */


  // ================== constructor ==================
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

  // ==================== variables ====================

  // pin variables
  int m_pwm_pin{0};
  const int m_brk_pin;
  int m_dir_pin{0};
  uint8_t CLOCKWISE{LOW};
  uint8_t COUNTER_CLOCKWISE{HIGH};
  int encoder_dir{1};

  // encoder variables
  Encoder m_encoder;
  float q[2] = {0, 0};  // rad

  // motor variables
  int CPR{1000};        // count per revolution
  float ki = 2.29183;   // V*s/rad
  float kt = 0.891514;  // Nm*ohms/V
  float R = 7.5;        // ohms


  // ===================== functions =====================

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
  float angle(){
    return (float)(encoder_dir*2*M_PI*m_encoder.read()/CPR);
  }

  // updates the angles to where q1 is always the most recent angle
  void update_angles(){
    // q[1] is the newest
    q[0] = q[1];
    q[1] = angle();
  }

  // takes the time step in millis
  // returns the current angular velocity in rad/s
  float qdot(float dt){
    return (float)(1000*(q[1] - q[0])/dt);
  }

  // takes the angular velocity in rad/s
  // returns the current emf voltage given dt in millis
  float Vemf(float qdot) {
    return ki*qdot;
  }

  // takes the a desired voltage to send to the motor, t
  // transforms it, then sends the signal to the driver
  void set_voltage(float Vs){

    // toggle the direction if needed
    set_direction(1);
    if (Vs < 0){
      Vs = -Vs;
      set_direction(-1);
    }

    // // constrain the voltage
    Vs = constrain(Vs, 0, 12);

    // map the voltage
    Vs = map(Vs, 0, 12, 0, 255);

    // send the signal tothe pin
    analogWrite(m_pwm_pin, Vs);
  }
};

////////////////////////   _____ _       _           _      //////////////////////
////////////////////////  / ____| |     | |         | |     //////////////////////
//////////////////////// | |  __| | ___ | |__   __ _| |___  //////////////////////
//////////////////////// | | |_ | |/ _ \| '_ \ / _` | / __| //////////////////////
//////////////////////// | |__| | | (_) | |_) | (_| | \__ \ //////////////////////
////////////////////////  \_____|_|\___/|_.__/ \__,_|_|___/ //////////////////////

// timer variables
long new_time{0};
long old_time{0};
long dt{0};


// define motors
Motor right_motor(3, 9, 12, HIGH, LOW,
                  20, 21, 1);

Motor left_motor(11, 8, 13, HIGH, LOW,
                  18, 19, 1);

// gains
float Kp{1};
float Kv{1};
float K{1};

// hardware controls
int button{2};

//////////////////////////   _____      _                //////////////////////////
//////////////////////////  / ____|    | |               //////////////////////////
////////////////////////// | (___   ___| |_ _   _ _ __   //////////////////////////
//////////////////////////  \___ \ / _ \ __| | | | '_ \  //////////////////////////
//////////////////////////  ____) |  __/ |_| |_| | |_) | //////////////////////////
////////////////////////// |_____/ \___|\__|\__,_| .__/  //////////////////////////
//////////////////////////                       | |     //////////////////////////
//////////////////////////                       |_|     //////////////////////////

void setup(){
   // start the coms
  Serial.begin(9600);

  // Setup the emergency stop function
  pinMode(button, INPUT);  // e-stop button is on pin 2
  attachInterrupt(
    digitalPinToInterrupt(2),
    [](){
      digitalWrite(right_motor.m_brk_pin, HIGH);
      digitalWrite(left_motor.m_brk_pin, HIGH);
    },
    FALLING  // trigger when it starts going
  );

  // Wait for user to push the start button
  while (digitalRead(button) != HIGH) {
    delay(250);
  }

  // Define the start position as in line with the N.y axis
  right_motor.m_encoder.write(right_motor.CPR/4);
  left_motor.m_encoder.write(left_motor.CPR/4);
}

/////////////////////////////  _                        /////////////////////////////
///////////////////////////// | |                       /////////////////////////////
///////////////////////////// | |     ___   ___  _ __   /////////////////////////////
///////////////////////////// | |    / _ \ / _ \| '_ \  /////////////////////////////
///////////////////////////// | |___| (_) | (_) | |_) | /////////////////////////////
///////////////////////////// |______\___/ \___/| .__/  /////////////////////////////
/////////////////////////////                   | |     /////////////////////////////
/////////////////////////////                   |_|     /////////////////////////////

void loop(){
  // ================= controller ==================

  // capture the new time
  new_time = millis();

  // calculate dt
  dt = (new_time - old_time);

  // get the newest angles
  right_motor.update_angles();
  left_motor.update_angles();
  Serial.print("q1 = "); Serial.print(right_motor.q[1]);
  Serial.print(", q4 = "); Serial.println(left_motor.q[1]);

  // calculate controller signal
  // ...calculations go here...

  // calculate Vs
  // ...Vs = (Tau*R/kt) + Vemf...;

  // output to the motors
  // right_motor.set_voltage(Vs);
  // left_motor.set_voltage(Vs);

  // // overwrite the old time
  old_time = new_time;
}
