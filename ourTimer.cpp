#include <chrono>
#include "ourTimer.h"

//The clock datatypes are really annoying in c++
//This file is basically just a wrapper to make timing easier to use

typedef std::chrono::time_point<std::chrono::steady_clock> clock_type;
clock_type start;

void startTimer(){ //THIS FUNCTION IS NOT THREAD SAFE OBVIOUSLY
  start = std::chrono::steady_clock::now();
}

float timeSince(){
  float ms = std::chrono::duration<float,std::milli>(std::chrono::steady_clock::now()-start).count();
  return ms*0.001; //convert to seconds
}

