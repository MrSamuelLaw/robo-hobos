#include <Arduino.h>
#include <Encoder.h>

// ============= define encoders =============
Encoder right_encoder(2, 3);
Encoder left_encoder(5, 6);
long right_old_pos{0};
long left_old_pos{0};
long right_new_pos{0};
long left_new_pos{0};
double left_deg{0};
double right_deg{0};
constexpr long CPR_left{250};
constexpr long CPR_right{1000};

// ================== setup ==================
void setup() {
  Serial.begin(9600);
}

// ================== loop ===================
void loop() {
  // read the encoder values
  left_new_pos = left_encoder.read();
  right_new_pos = right_encoder.read();

  // check if any new values
  if ((left_new_pos != left_old_pos) || (right_new_pos != right_old_pos)) {
    // update the values
    left_old_pos = left_new_pos;
    right_old_pos = right_new_pos;
    // print out the angles
    left_deg = (double)(360*left_new_pos/CPR_left);
    right_deg = (double)(360*right_new_pos/CPR_right);
    Serial.print("left angle = ");
    Serial.println(left_deg);
    Serial.print("right angle = ");
    Serial.println(right_deg);
  }
}

