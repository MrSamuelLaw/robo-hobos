#include <Arduino.h>
#include <BasicLinearAlgebra.h>
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
  float read_angle(){
    return (float)(encoder_dir*2*M_PI*m_encoder.read()/CPR);
  }

  // reads the newest angle, shifts the vector,
  // and returns the new angle
  float update_angles(){
    // q[1] is the newest
    q[0] = q[1];
    q[1] = read_angle();
    return q[1];
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

    // constrain the voltage
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
Motor right_motor(3, 9, 12, LOW, HIGH,
                  20, 21, 1);

Motor left_motor(11, 8, 13, LOW, HIGH,
                  18, 19, 1);

// controller variables
float Kp{400};                      // proportional gain
float Kv{sqrt(Kp)*sqrt(8)*1.3};     // derivative gain
float K{2};   // constant gain
float q_right_desired{1.30823735};  // radians
float q_left_desired{1.81830097};   // radians

// hardware controls
int button{2};

// constants
constexpr float LA_1 {0.06095};
constexpr float LC_1 {0.06095};
constexpr float LB_1 {0.08128};
constexpr float LD_1 {0.08128};
constexpr float L {0.1016};
constexpr float LA_COM {0.027};
constexpr float LC_COM {0.027};
constexpr float LB_COM {0.04};
constexpr float LD_COM {0.04};
constexpr float bodyA_mass {0.007};
constexpr float bodyB_mass {0.008};
constexpr float bodyC_mass {0.007};
constexpr float bodyD_mass {0.008};
constexpr float bodyA_izz {6.36e-6};
constexpr float bodyB_izz {1.01e-5};
constexpr float bodyC_izz {6.36e-6};
constexpr float bodyD_izz {1.01e-5};
constexpr float motorR_izz {3.3e-07};
constexpr float motorL_izz {3.3e-07};
constexpr float gr {144};


// variable declarations
float q1{0};
float q2{0};
float q4{0};
float q5{0};
float q1d{0};
float q2d{0};
float q4d{0};
float q5d{0};
float q1dd{0};
float q2dd{0};
float q4dd{0};
float q5dd{0};
float x{0};
float y{0};
float xd{0};
float yd{0};
float xdd{0};
float ydd{0};
float x_desired{0};
float y_desired{0};
float xd_desired{0};
float yd_desired{0};
float gamma{0};
float alpha1{0};
float alpha2{0};
BLA::Matrix<4, 4> J;
BLA::Matrix<4, 4> Jd;
BLA::Matrix<4, 1> Xdd;
BLA::Matrix<4, 1> Qd;
BLA::Matrix<4, 1> Qdd;
float tau1{0};
float tau2{0};
float Vs_right{0};
float Vs_left{0};
float x_des_vec[4] = {-0.0762, -0.0762, -0.0254, -0.0254};
float y_des_vec[4] = {0.07366, 0.09906, 0.09906, 0.07366};
int head{0};
float eps = 0.01;

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
  Serial.println("Press and hold button to start...");
  while (digitalRead(button) != HIGH) {
    delay(100);
  }

  // Define the start position as in line with the N.y axis
  right_motor.m_encoder.write(right_motor.CPR/4);
  left_motor.m_encoder.write(left_motor.CPR/4);

  // Zero fill the various matricies
  J.Fill(0);
  Jd.Fill(0);

  // Set desired x, y final position
  x_desired = x_des_vec[head];
  y_desired = y_des_vec[head];
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

  // calculate angles
  q1 = right_motor.update_angles();
  q4 = left_motor.update_angles();

  gamma = acos((1.0/2.0)*sqrt(pow(L, 2) + 2*L*LA_1*cos(q1) - 2*L*LC_1*cos(q4)
          + pow(LA_1, 2) - 2*LA_1*LC_1*(sin(q1)*sin(q4) + cos(q1)*cos(q4))
          + pow(LC_1, 2))/LB_1);

  alpha1 = acos((L*LA_1*cos(q1) + pow(LA_1, 2) - LA_1*LC_1*(sin(q1)*sin(q4)
           + cos(q1)*cos(q4)))/(sqrt(pow(L, 2) + 2*L*LA_1*cos(q1) - 2*L*LC_1*cos(q4)
           + pow(LA_1, 2) - 2*LA_1*LC_1*(sin(q1)*sin(q4) + cos(q1)*cos(q4))
           + pow(LC_1, 2))*sqrt(pow(LA_1, 2))));

  alpha2 = q1 - q4 - alpha1 + M_PI;

  q2 = M_PI - (alpha1 + gamma);
  q5 = M_PI + (alpha2 + gamma);

  // calculate angular velocities
  q1d = right_motor.qdot(dt);
  q4d = left_motor.qdot(dt);
  q2d = 0.25*(3.0*q1d*sin(-q1 + q4 + q5) - 4.0*q1d*sin(q1 + q2 - q4 - q5)
        - 3.0*q4d*sin(q5))/sin(q1 + q2 - q4 - q5);
  q5d = 0.25*(3.0*q1d*sin(q2) - 3.0*q4d*sin(q1 + q2 - q4)
        - 4.0*q4d*sin(q1 + q2 - q4 - q5))/sin(q1 + q2 - q4 - q5);

  // calculate the task space coordinates
  x = LA_1*cos(q1) + LB_1*(-sin(q1)*sin(q2) + cos(q1)*cos(q2));
  y = LA_1*sin(q1) + LB_1*(sin(q1)*cos(q2) + sin(q2)*cos(q1));

  // calculate the task space velocities
  xd = -LA_1*q1d*sin(q1) + LB_1*(q1d + q2d)*(-sin(q1)*cos(q2) - sin(q2)*cos(q1));
  yd = LA_1*q1d*cos(q1) + LB_1*(q1d + q2d)*(-sin(q1)*sin(q2) + cos(q1)*cos(q2));

  // calculate the controller signal
  xdd = K*(Kp*(-x + x_desired) + Kv*(-xd + xd_desired));
  ydd = K*(Kp*(-y + y_desired) + Kv*(-yd + yd_desired));

  // create the jacobian
  J(0, 0) = -LA_1*sin(q1) + LB_1*(-sin(q1)*cos(q2) - sin(q2)*cos(q1));
  J(1, 0) = LA_1*cos(q1) + LB_1*(-sin(q1)*sin(q2) + cos(q1)*cos(q2));
  J(0, 1) = LB_1*(-sin(q1)*cos(q2) - sin(q2)*cos(q1));
  J(1, 1) = LB_1*(-sin(q1)*sin(q2) + cos(q1)*cos(q2));
  J(2, 2) = -LC_1*sin(q4) + LD_1*(-sin(q4)*cos(q5) - sin(q5)*cos(q4));
  J(3, 2) = LC_1*cos(q4) + LD_1*(-sin(q4)*sin(q5) + cos(q4)*cos(q5));
  J(2, 3) = LD_1*(-sin(q4)*cos(q5) - sin(q5)*cos(q4));
  J(3, 3) = LD_1*(-sin(q4)*sin(q5) + cos(q4)*cos(q5));

  // create the jocobian's derivative
  Jd(0, 0) = -LA_1*q1d*cos(q1) + LB_1*(q1d*sin(q1)*sin(q2) - q1d*cos(q1)*cos(q2)
             + q2d*sin(q1)*sin(q2) - q2d*cos(q1)*cos(q2));
  Jd(1, 0) = -LA_1*q1d*sin(q1) + LB_1*(-q1d*sin(q1)*cos(q2) - q1d*sin(q2)*cos(q1)
             - q2d*sin(q1)*cos(q2) - q2d*sin(q2)*cos(q1));
  Jd(0, 1) = LB_1*(q1d*sin(q1)*sin(q2) - q1d*cos(q1)*cos(q2) + q2d*sin(q1)*sin(q2)
             - q2d*cos(q1)*cos(q2));
  Jd(1, 1) = LB_1*(-q1d*sin(q1)*cos(q2) - q1d*sin(q2)*cos(q1) - q2d*sin(q1)*cos(q2)
             - q2d*sin(q2)*cos(q1));
  Jd(2, 2) = -LC_1*q4d*cos(q4) + LD_1*(q4d*sin(q4)*sin(q5) - q4d*cos(q4)*cos(q5)
             + q5d*sin(q4)*sin(q5) - q5d*cos(q4)*cos(q5));
  Jd(3, 2) = -LC_1*q4d*sin(q4) + LD_1*(-q4d*sin(q4)*cos(q5) - q4d*sin(q5)*cos(q4)
             - q5d*sin(q4)*cos(q5) - q5d*sin(q5)*cos(q4));
  Jd(2, 3) = LD_1*(q4d*sin(q4)*sin(q5) - q4d*cos(q4)*cos(q5) + q5d*sin(q4)*sin(q5)
             - q5d*cos(q4)*cos(q5));
  Jd(3, 3) = LD_1*(-q4d*sin(q4)*cos(q5) - q4d*sin(q5)*cos(q4) - q5d*sin(q4)*cos(q5)
             - q5d*sin(q5)*cos(q4));

  // calculate the desired joint accelerations
  Xdd << xdd, ydd, xdd, ydd;
  Qd << q1d, q2d, q4d, q5d;
  Qdd = J.Inverse()*(Xdd - Jd*Qd);

  // assign them to their own variabels
  q1dd = Qdd(0);
  q2dd = Qdd(1);
  q4dd = Qdd(2);
  q5dd = Qdd(3);

  // calculate the desired torque
  tau1 = -LA_1*LB_COM*bodyB_mass*pow(q1d, 2)*sin(q2) + LA_1*LB_COM*bodyB_mass*pow(q1d + q2d, 2)*sin(q2) - q1dd*(pow(LA_COM, 2)
         *bodyA_mass + bodyA_izz + bodyB_izz + bodyB_mass*(pow(LA_1, 2) + 2*LA_1*LB_COM*cos(q2) + pow(LB_COM, 2)) + pow(gr, 2)
         *motorR_izz) - q2dd*(bodyB_izz + bodyB_mass*(LA_1*LB_COM*cos(q2) + pow(LB_COM, 2))) - (LA_1*cos(q1) + LB_1*(-sin(q1)
         *sin(q2) + cos(q1)*cos(q2)) - (-LA_1*sin(q1) + LB_1*(-sin(q1)*cos(q2) - sin(q2)*cos(q1)))*(-sin(q1)*sin(q2) + cos(q1)
         *cos(q2))/(-sin(q1)*cos(q2) - sin(q2)*cos(q1)))*(-LC_1*LD_COM*bodyD_mass*pow(q4d, 2)*sin(q5) - q4dd*(bodyD_izz +
         bodyD_mass*(LC_1*LD_COM*cos(q5) + pow(LD_COM, 2))) - q5dd*(pow(LD_COM, 2)*bodyD_mass + bodyD_izz))/(LD_1*(-sin(q1)*sin
         (q2) + cos(q1)*cos(q2))*(-sin(q4)*cos(q5) - sin(q5)*cos(q4))/(-sin(q1)*cos(q2) - sin(q2)*cos(q1)) - LD_1*(-sin(q4)*sin
         (q5) + cos(q4)*cos(q5))) - (-LA_1*sin(q1) + LB_1*(-sin(q1)*cos(q2) - sin(q2)*cos(q1)) + LD_1*(-sin(q4)*cos(q5) - sin
         (q5)*cos(q4))*(LA_1*cos(q1) + LB_1*(-sin(q1)*sin(q2) + cos(q1)*cos(q2)) - (-LA_1*sin(q1) + LB_1*(-sin(q1)*cos(q2) -
         sin(q2)*cos(q1)))*(-sin(q1)*sin(q2) + cos(q1)*cos(q2))/(-sin(q1)*cos(q2) - sin(q2)*cos(q1)))/(LD_1*(-sin(q1)*sin(q2)
         + cos(q1)*cos(q2))*(-sin(q4)*cos(q5) - sin(q5)*cos(q4))/(-sin(q1)*cos(q2) - sin(q2)*cos(q1)) - LD_1*(-sin(q4)*sin(q5)
         + cos(q4)*cos(q5))))*(-LA_1*LB_COM*bodyB_mass*pow(q1d, 2)*sin(q2) - q1dd*(bodyB_izz + bodyB_mass*(LA_1*LB_COM*cos(q2)
         + pow(LB_COM, 2))) - q2dd*(pow(LB_COM, 2)*bodyB_mass + bodyB_izz))/(LB_1*(-sin(q1)*cos(q2) - sin(q2)*cos(q1)));

  tau2 = -LC_1*LD_COM*bodyD_mass*pow(q4d, 2)*sin(q5) + LC_1*LD_COM*bodyD_mass*pow(q4d + q5d, 2)*sin(q5) - q4dd*(pow(LC_COM, 2)
          *bodyC_mass + bodyC_izz + bodyD_izz + bodyD_mass*(pow(LC_1, 2) + 2*LC_1*LD_COM*cos(q5) + pow(LD_COM, 2)) + pow(gr, 2)
          *motorL_izz) - q5dd*(bodyD_izz + bodyD_mass*(LC_1*LD_COM*cos(q5) + pow(LD_COM, 2))) - (-LC_1*cos(q4) - LD_1*(-sin(q4)
          *sin(q5) + cos(q4)*cos(q5)) - (LC_1*sin(q4) - LD_1*(-sin(q4)*cos(q5) - sin(q5)*cos(q4)))*(-sin(q1)*sin(q2) + cos(q1)
          *cos(q2))/(-sin(q1)*cos(q2) - sin(q2)*cos(q1)))*(-LC_1*LD_COM*bodyD_mass*pow(q4d, 2)*sin(q5) - q4dd*(bodyD_izz +
          bodyD_mass*(LC_1*LD_COM*cos(q5) + pow(LD_COM, 2))) - q5dd*(pow(LD_COM, 2)*bodyD_mass + bodyD_izz))/(LD_1*(-sin(q1)*
          sin(q2) + cos(q1)*cos(q2))*(-sin(q4)*cos(q5) - sin(q5)*cos(q4))/(-sin(q1)*cos(q2) - sin(q2)*cos(q1)) - LD_1*(-sin(q4)
          *sin(q5) + cos(q4)*cos(q5))) - (LC_1*sin(q4) - LD_1*(-sin(q4)*cos(q5) - sin(q5)*cos(q4)) + LD_1*(-sin(q4)*cos(q5)
          -sin(q5)*cos(q4))*(-LC_1*cos(q4) - LD_1*(-sin(q4)*sin(q5) + cos(q4)*cos(q5)) - (LC_1*sin(q4) - LD_1*(-sin(q4)*
          cos(q5) - sin(q5)*cos(q4)))*(-sin(q1)*sin(q2) + cos(q1)*cos(q2))/(-sin(q1)*cos(q2) - sin(q2)*cos(q1)))/(LD_1*
          (-sin(q1)*sin(q2) + cos(q1)*cos(q2))*(-sin(q4)*cos(q5) - sin(q5)*cos(q4))/(-sin(q1)*cos(q2) - sin(q2)*cos(q1))
          - LD_1*(-sin(q4)*sin(q5) + cos(q4)*cos(q5))))*(-LA_1*LB_COM*bodyB_mass*pow(q1d, 2)*sin(q2) - q1dd*(bodyB_izz +
          bodyB_mass*(LA_1*LB_COM*cos(q2) + pow(LB_COM, 2))) - q2dd*(pow(LB_COM, 2)*bodyB_mass + bodyB_izz))/(LB_1*(-sin(q1)
          *cos(q2) - sin(q2)*cos(q1)));

  // calculate Vs for the right and left motor
  // ...Vs = (Tau*R/kt) + Vemf...;
  Vs_right = -((tau1*right_motor.R/right_motor.kt) + (right_motor.ki*q1d));
  Vs_left = -((tau2*left_motor.R/left_motor.kt) + (left_motor.ki*q4d));

  // output to the motors
  right_motor.set_voltage(Vs_right);
  left_motor.set_voltage(Vs_left);

  // change target if the x or y values are with in range
  if ((fabs(x - x_desired) <= eps) && (fabs(y - y_desired) <= eps)) {
    head += 1;
    if (head == 4) {head = 0;}
    x_desired = x_des_vec[head];
    y_desired = y_des_vec[head];
  }

  // overwrite the old time
  old_time = new_time;

  // print out for debugging
  Serial.println(head);
  // Serial.print("x= "); Serial.print(x);
  // Serial.print(", y= "); Serial.println(y);
}
