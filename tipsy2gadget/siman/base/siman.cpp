#include "siman.hpp"


CScripted siman::config("config",true);

bool siman::fileExists(string filename) {
  ifstream fs(filename.c_str());
  return(fs.is_open());
}
