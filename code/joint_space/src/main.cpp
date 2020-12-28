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
double Kp{375};                      // proportional gain
double Kv{sqrt(Kp*8)*0.1};           // derivative gain
double K{2.5};                       // constant gain
double Ks{3.0};                      // voltage gain

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
double q1_desired{0};
double q2_desired{0};
double q4_desired{0};
double q5_desired{0};
double gamma{0};
double alpha1{0};
double alpha2{0};
double tau1{0};
double tau2{0};
double Vs_right{0};
double Vs_left{0};


// ======================================= target angles =======================================
int head{0};
constexpr int vec_size{75};
double eps = 0.05;

double q1_des_vec[vec_size] = {
0.9140829292879634,0.8780748888112396,0.8765237998246486,0.9117214846703798,0.9794736769381787,
1.0700335584853868,1.1717506400091862,1.2748073938384386,1.3726842732247824,1.4614019855374178,
1.5376777428257513,1.5965587618246606,1.6287286582151137,1.6205960174754594,1.5653189721558074,
1.474950622633884,1.3706641873354002,1.2702763816028337,1.1869514948161128,1.130846871179979,
1.109366910770917,1.1258961519004627,1.1782262590223704,1.2585596119913167,1.3560087324008454,
1.4601026569492255,1.562922191839539,1.6590292416728052,1.7437019982978705,1.809550215102155,
1.8409040220700112,1.8144529946506578,1.7279389914405,1.608505746060645,1.4809950464837305,
1.3606259505080478,1.258030864784392,1.1819801558964707,1.139585633103628,1.1348459860162343,
1.1665126984693526,1.2271133951480995,1.304747761168127,1.3868314868201967,1.4630478178423414,
1.5259842950061204,1.5701413043159202,1.5907235109606421,1.5833479958808463,1.5453102048117178,
1.4777745994285287,1.3867023176107216,1.2814056357463746,1.1723716830545148,1.070046556416459,
0.9845181300064579,0.9251990457961629,0.8996741926374362,0.9114154925850264,0.9573635373447067,
1.027741005348061,1.1093644426983262,1.1902031973163665,1.261867272799953,1.319370706311709,
1.3596446079586755,1.380174230066638,1.3784244263174659,1.3523683977893048,1.3020243900186868,
1.2309841849022127,1.1465541625826963,1.058474094501223,0.9774363435804101,0.9140829293105299
};

double q2_des_vec[vec_size] = {
1.512005955975962,1.6001186599615926,1.6620537272463465,1.6881505038905453,1.6723483349469537,
1.6136256504338453,1.5158995195794474,1.3864673333652855,1.2342697541750154,1.0691143228126643,
0.9023309488028849,0.7491341453422273,0.6324363778746434,0.582336886851152,0.616713304702449,
0.7200292337372963,0.8606991204272497,1.0131214535752715,1.1597435093812574,1.2872894937753019,
1.3845157635035434,1.4420197330607591,1.4534586711345134,1.4170429948695018,1.3359586301821191,
1.2172767524730503,1.0702730007884442,0.9054961070503568,0.7356292212204733,0.5795369286930956,
0.4709716973952665,0.4568887271603801,0.546935554573628,0.698913349505417,0.8746627440650226,
1.053157176694752,1.2216132273668217,1.3700696086884605,1.4894756122261363,1.5717406610557758,
1.6110986018493314,1.605956403366329,1.5598986209706234,1.4810378918961493,1.380342336984325,
1.2702244092443638,1.1639880370114803,1.0756898197939577,1.0192140249707855,1.0053755898647592,
1.0378917907493328,1.1118529298027298,1.2165069361986747,1.3393618538013106,1.4684946133590302,
1.592937529580173,1.702319394663482,1.7867583240337706,1.8376262813073259,1.8492555503464612,
1.8206992957245398,1.7561328497028392,1.6635162534825612,1.5526649381087685,1.4339777540660381,
1.3181018471793287,1.2160594513042597,1.1390424396330492,1.0970218806405507,1.095999787214903,
1.1354763456580432,1.2084380743048642,1.3038268439349745,1.4091842659284395,1.5120059559608572
};

double q4_des_vec[vec_size] = {
1.4501431919842653,1.5105053544277935,1.5881678655433786,1.676889879040713,1.7696038696766312,
1.8569219002756099,1.9272314047066292,1.9690393105634487,1.9744642718687353,1.9413683434389897,
1.8728933495972788,1.7755604470406756,1.6577774868569162,1.5301130535105898,1.4094431261526181,
1.3296125523781954,1.3263445145691906,1.3804903143125107,1.4604707155922145,1.5534808434261353,
1.6549648064898665,1.7620175912460359,1.8705295253846754,1.973613773631946,2.0614035387680656,
2.1230178023792714,2.150222905771786,2.1404101821761565,2.0969773893719164,2.0277330999528567,
1.9431500360916591,1.8553969330515525,1.777616877813887,1.7219186347994069,1.6956099341165247,
1.6990004257421927,1.7277833074489426,1.7769656283073063,1.842690713089488,1.9219264877770064,
2.011099627982356,2.1043530583822236,2.192265462007121,2.2624944312141886,2.3032134998362275,
2.3075258432732086,2.2754837028941344,2.2128035920836426,2.128288557197152,2.0318671939961686,
1.933619840734234,1.8432356963513257,1.7692345579015716,1.7177810139717204,1.6917709261069733,
1.6910695933864621,1.7137170312574441,1.7569966466292481,1.8175781143128888,1.890677308681323,
1.9687165571027725,2.0405931083763904,2.0930794985973984,2.1146615253465795,2.0994329281485316,
2.048175555042894,1.9666544074550645,1.8633090256939953,1.7479430358979713,1.6316959329624485,
1.527711521341505,1.4506822424433459,1.4122575096848415,1.4144011288341014,1.4501431919728236
};

double q5_des_vec[vec_size] = {
-0.9202370900305313,-1.0501188161336188,-1.1835896343364918,-1.3056754777694455,-1.403640026672514,
-1.4665855964238714,-1.4863615847895764,-1.4590832187408411,-1.3859327185225954,-1.272401432852303,
-1.1265194945228225,-0.9573660926956098,-0.7748695434933234,-0.5919764536108529,-0.4327157281342907,
-0.3510311047862648,-0.4040789401553962,-0.5524745705023526,-0.7346153549337585,-0.9215900572879615,
-1.09921598536357,-1.2575008470372626,-1.3874800336388668,-1.4808371010756118,-1.5310258596942605,
-1.5351256844754455,-1.4952065776264236,-1.418167560581888,-1.314314466581098,-1.1959422522485055,
-1.076844077790058,-0.9726171646389766,-0.9004183999927079,-0.8757348137776664,-0.9056878245513432,
-0.9849667298859952,-1.0998264196236214,-1.234927214974529,-1.376728844883258,-1.5136486868086885,
-1.6352434649517125,-1.7317523264518682,-1.7946026386286498,-1.8179578256156972,-1.8005253769031266,
-1.746207546653974,-1.6629762828687462,-1.5609292869741185,-1.4508402740847768,-1.3436148421688634,
-1.250221553850163,-1.1813343150521545,-1.1459874138378183,-1.1492757597737946,-1.1904813360798772,
-1.263281752959082,-1.3578201521162214,-1.4628603493272954,-1.5669137844536054,-1.6585512676415652,
-1.7266904029035552,-1.7615689583599459,-1.7565268623731254,-1.709715448336642,-1.6244072479875675,
-1.5076551532596272,-1.368473207098514,-1.2167901473101246,-1.0636171591134775,-0.922290872992877,
-0.8100571022158749,-0.7472049470050838,-0.7487451292533358,-0.8124861996450135,-0.9202370900120341
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

  // set up the initial targets
  q1_desired = q1_des_vec[head];
  q2_desired = q2_des_vec[head];
  q4_desired = q4_des_vec[head];
  q5_desired = q5_des_vec[head];
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


  // ============================== joint space code ==============================

  // calculate the desired accellerations
  q1dd = K*(Kp*(q1_desired - q1)) + (Kv*(0 - q1d));
  q2dd = K*(Kp*(q2_desired - q2)) + (Kv*(0 - q2d));
  q4dd = K*(Kp*(q4_desired - q4)) + (Kv*(0 - q4d));
  q5dd = K*(Kp*(q5_desired - q5)) + (Kv*(0 - q5d));

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
  right_motor.set_voltage(Ks*Vs_right);
  left_motor.set_voltage(Ks*Vs_left);

  // change target if the q1 and q4 values are with in range
  if ((fabs(q1_desired - q1) <= eps) && (fabs(q4_desired - q4) <= eps)) {
    head += 1;
    if (head == vec_size) {head = 0;}
    q1_desired = q1_des_vec[head];
    q4_desired = q4_des_vec[head];
  }

  // overwrite the old time
  old_time = new_time;

  // print out for debugging
  Serial.print(head);
  Serial.print(","); Serial.print(q1);
  Serial.print(","); Serial.println(q4);
}
