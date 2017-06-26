#include <iostream>
#include <chrono>
#include <thread>

extern "C" void cppfunction_()
{

  std::cout<<"Hello from Cpp "<<std::endl;
  std::chrono::time_point<std::chrono::system_clock> tStart = std::chrono::system_clock::now();

  std::this_thread::sleep_for(std::chrono::milliseconds(1));

  std::chrono::time_point<std::chrono::system_clock> tEnd = std::chrono::system_clock::now();
  std::cout<<"duration: "<<std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd-tStart).count()<<std::endl;
}
