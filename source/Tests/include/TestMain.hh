#pragma once

#include "TestUtilities.hh"

int runTest(int argn, char** argc);

int main(int argn, char** argc) {
  auto &logger = streamlog::logstream::global() ;
  logger.setName( "Testing" ) ;
  logger.setSinks( { streamlog::logstream::console() } ) ;
  logger.setLevel("DEBUG6");

  try {
    return runTest(argn, argc);
  } catch (std::out_of_range& e) {
    std::cerr << "Geometry does not agree with energy in the trees:" << e.what() << std::endl;
    return 1;
  } catch (std::runtime_error& e) {
    std::cerr << "Runtime Error: " << e.what() << std::endl;
    return 1;
  } catch (std::invalid_argument& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  return 0;
}
