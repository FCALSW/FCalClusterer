#ifndef Global_hh
#define Global_hh 1

#define _VERY_SMALL_NUMBER 1e-9
#define _CLUSTER_MERGE_DEBUG 1
#define _GLOBAL_COUNTERS_UPDATE_DEBUG 1
#define _BHABHA_SELECTION_CUTS_DEBUG 1
#define _CLUSTER_RESET_STATS_DEBUG 1

#include "LCCluster.hh"

#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <iomanip>

template<class T, class G>
void printMap (std::string callFrom, int line , std::map< T, G> const& /*myMap*/) {
  std::cout << __PRETTY_FUNCTION__  << std::endl;
  std::cout << "Calling from " << callFrom << " " << line  << std::endl;
  // for( auto& entry: myMap ){
  //   std::cout << entry.first << "  " << entry.second  << std::endl;
  // }
  return;
}


template<class T>
void printVector (std::vector<T> const& myVector ) {
  // for( auto& entry: myVector ){
  //   std::cout << std::setw(14) << entry;
  // }
  std::cout << std::endl;
  return;
}


template<class T, class G>
void printMapVector (std::string callFrom, int line , std::map< T, std::vector<  G > > const& myMap) {
  std::cout << "Calling from " << callFrom << " " << line  << std::endl;
  // for( auto const& entry: myMap ){
  //   std::cout << entry.first;
  //   printVector(entry.second);
  // }
  return;
}

template<class T>
void printMapVector (std::string callFrom, int line , std::map< T, LCCluster > const& myMap) {
  std::cout << "Calling from " << callFrom << " " << line  << std::endl;
  // for( auto const& entry: myMap ){
  //   std::cout << std::setw(4) << entry.first << entry.second << std::endl;
  // }
  return;
}

#endif // Global_hh
