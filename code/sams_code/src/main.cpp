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

  // ================ internal variables ================

  // pin variables
  int m_pwm_pin{0};
  const int m_brk_pin;
  int m_dir_pin{0};
  uint8_t CLOCKWISE{LOW};
  uint8_t COUNTER_CLOCKWISE{HIGH};
  int encoder_dir{1};

  // encoder variables
  Encoder m_encoder;
  double q[2] = {0, 0};  // rad

  // motor variables
  int CPR{1000};        // count per revolution
  double ki = 2.29183;   // V*s/rad
  double kt = 0.891514;  // Nm*ohms/V
  double R = 7.5;        // ohms


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
  double read_angle(){
    return (double)(encoder_dir*2*M_PI*m_encoder.read()/CPR);
  }

  // reads the newest angle, shifts the vector,
  // and returns the new angle
  double update_angles(){
    // q[1] is the newest
    q[0] = q[1];
    q[1] = read_angle();
    return q[1];
  }

  // takes the time step in millis
  // returns the current angular velocity in rad/s
  double qdot(double dt){
    return (double)(1000*(q[1] - q[0])/dt);
  }

  // takes the angular velocity in rad/s
  // returns the current emf voltage given dt in millis
  double Vemf(double qdot) {
    return ki*qdot;
  }

  // takes the a desired voltage to send to the motor, t
  // transforms it, then sends the signal to the driver
  void set_voltage(double Vs){
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
                  18, 19, 1);

Motor left_motor(11, 8, 13, LOW, HIGH,
                  20, 21, 1);

// controller variables
double Kp{200};                      // proportional gain
double Kv{sqrt(Kp*8)*0.1};           // derivative gain
double K{2.5};                       // voltage gain
double q_right_desired{1.30823735};  // radians
double q_left_desired{1.81830097};   // radians

// hardware controls
int button{2};

// constants
constexpr double LA_1 {0.06095};
constexpr double LC_1 {0.06095};
constexpr double LB_1 {0.08128};
constexpr double LD_1 {0.08128};
constexpr double L {0.1016};
constexpr double LA_COM {0.027};
constexpr double LC_COM {0.027};
constexpr double LB_COM {0.04};
constexpr double LD_COM {0.04};
constexpr double bodyA_mass {0.007};
constexpr double bodyB_mass {0.008};
constexpr double bodyC_mass {0.007};
constexpr double bodyD_mass {0.008};
constexpr double bodyA_izz {6.36e-6};
constexpr double bodyB_izz {1.01e-5};
constexpr double bodyC_izz {6.36e-6};
constexpr double bodyD_izz {1.01e-5};
constexpr double motorR_izz {3.3e-07};
constexpr double motorL_izz {3.3e-07};
constexpr double gr {144};


// variable declarations
double q1{0};
double q2{0};
double q4{0};
double q5{0};
double q1d{0};
double q2d{0};
double q4d{0};
double q5d{0};
double q1dd{0};
double q2dd{0};
double q4dd{0};
double q5dd{0};
double x{0};
double y{0};
double xd{0};
double yd{0};
double xdd{0};
double ydd{0};
double x_desired{0};
double y_desired{0};
double xd_desired{0};
double yd_desired{0};
double gamma{0};
double alpha1{0};
double alpha2{0};
BLA::Matrix<4, 4> J;
BLA::Matrix<4, 4> Jd;
BLA::Matrix<4, 1> Xdd;
BLA::Matrix<4, 1> Qd;
BLA::Matrix<4, 1> Qdd;
double tau1{0};
double tau2{0};
double Vs_right{0};
double Vs_left{0};
// ============================ target coordinates ============================
int head{0};
constexpr int vec_size{75};
double eps = 0.005;

double x_des_vec[vec_size] = {
-0.024130000000000002, -0.025109694891845433, -0.027940772047069844, -0.03231152000762576, -0.037742050406485934,
-0.043638982741677504, -0.04936304043722211, -0.05430183707487175, -0.05793953911243276, -0.05991544862125914,
-0.06006481576422884, -0.05843722123776428, -0.055290431739842306, -0.051060433906428296, -0.04631107392098367,
-0.04166906019567737, -0.037751758933673205, -0.03509603594308051, -0.034096278695686474, -0.034958683944632196,
-0.037677038040337056, -0.04203276169998632, -0.04761921890801684, -0.053887517549479486, -0.06020857400112338,
-0.0659443557417252, -0.07052016738026737, -0.07348972614351527, -0.07458559645041507, -0.07374922564018348,
-0.07113715314750665, -0.06710268721665746, -0.06215514566636449, -0.05690132009825523, -0.05197585348931449,
-0.0479684873453353, -0.04535649185902958, -0.04445, -0.04535649185902955, -0.04796848734533526,
-0.051975853489314404, -0.056901320098255125, -0.06215514566636439, -0.06710268721665741, -0.0711371531475066,
-0.07374922564018344, -0.07458559645041507, -0.07348972614351532, -0.07052016738026741, -0.06594435574172532,
-0.060208574001123484, -0.0538875175494796, -0.0476192189080169, -0.04203276169998642, -0.03767703804033712,
-0.03495868394463221, -0.03409627869568646, -0.03509603594308047, -0.03775175893367311, -0.041669060195677274,
-0.04631107392098358, -0.051060433906428226, -0.055290431739842265, -0.05843722123776423, -0.0600648157642288,
-0.05991544862125915, -0.05793953911243282, -0.054301837074871837, -0.04936304043722219, -0.04363898274167759,
-0.03774205040648601, -0.0323115200076259, -0.027940772047069913, -0.02510969489184546, -0.024130000000000002
};

double y_des_vec[vec_size] = {
0.1016, 0.0969614901814798, 0.09294513617911893, 0.09010135103969591, 0.088845263319553,
0.0894086404960349, 0.09181277471109715, 0.09586543310033867, 0.10118222556401539, 0.1072299530957212,
0.11338698841483608, 0.1190137945636598, 0.12352552866172861, 0.12645843950934463, 0.1275224764493718,
0.12663410186561616, 0.12392556098743099, 0.11972955194459331, 0.1145410490401533, 0.10896064195312585,
0.10362586489986711, 0.09913836146798602, 0.09599520618428542, 0.09453222858529324, 0.09488581392825457,
0.09697754341701223, 0.10052342712539639, 0.10506667282789922, 0.11003024472190248, 0.11478320478103199,
0.11871325455397722, 0.12129718654719589, 0.122161192868283, 0.12112413730960561, 0.11821884316764986,
0.11368896048982106, 0.1079617662038626, 0.10160000000000004, 0.09523823379613751, 0.089511039510179,
0.08498115683235019, 0.08207586269039441, 0.08103880713171698, 0.08190281345280408, 0.08448674544602272,
0.08841679521896788, 0.0931697552780974, 0.09813332717210065, 0.10267657287460354, 0.10622245658298769,
0.1083141860717454, 0.10866777141470675, 0.1072047938157146, 0.10406163853201406, 0.099574135100133,
0.09423935804687424, 0.08865895095984677, 0.08347044805540675, 0.07927443901256909, 0.07656589813438387,
0.07567752355062818, 0.07674156049065534, 0.07967447133827134, 0.0841862054363401, 0.08981301158516376,
0.09597004690427871, 0.10201777443598448, 0.10733456689966123, 0.1113872252889028, 0.11379135950396507,
0.114354736680447, 0.11309864896030417, 0.11025486382088114, 0.10623850981852027, 0.10160000000000007
};

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

  gamma = acos((1.0/2.0)*sqrt(pow(L, 2) + 2*L*LA_1*cos(q1) - 2*L*LC_1*cos(q4) + pow(LA_1, 2)
          - 2*LA_1*LC_1*(sin(q1)*sin(q4) + cos(q1)*cos(q4)) + pow(LC_1, 2))/LB_1);

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
  right_motor.set_voltage(K*Vs_right);
  left_motor.set_voltage(K*Vs_left);

  // change target if the x or y values are with in range
  if (sqrt((x_desired - x)*(x_desired - x) + (y_desired - y)*(y_desired - y)) < eps) {
    head += 1;
    if (head == vec_size) {head = 0;}
    x_desired = x_des_vec[head];
    y_desired = y_des_vec[head];
  }

  // overwrite the old time
  old_time = new_time;

  // print out for debugging
  Serial.print(head);
  Serial.print(","); Serial.print(x);
  Serial.print(","); Serial.println(y);
}
