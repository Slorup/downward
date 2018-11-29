#include "collections.h"

namespace utils {

template <> std::ostream& operator <<<bool> (std::ostream& os, const std::vector<bool>& v){
  os << "[";
  for (size_t element = 0; element < v.size(); element++){
    if(v[element]){
      os << element << ",";
    }
  }
  os << "]";
  return os;
}
}
