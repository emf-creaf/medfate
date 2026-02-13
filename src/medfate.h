#include "RcppArmadillo.h"
#include <vector>
#include <string>
#include <limits>
#include <cmath>

#ifndef MEDFATE_H
#define MEDFATE_H

// ============================================================================
// PURE C++ STRUCTURES (Internal use)
// ============================================================================

class AbstractModelInput {
public:
  std::string input_class;
  std::string version;
  AbstractModelInput(Rcpp::List x);
  
  std::string getInputClass() {return(input_class);}
  void checkInputClass(Rcpp::List x);
  virtual ~AbstractModelInput() = default;
};

struct ABSTRACTMODEL_RESULT {
  virtual ~ABSTRACTMODEL_RESULT() = default;
};


namespace medfate {
    constexpr double NA_DOUBLE = std::numeric_limits<double>::quiet_NaN();
    constexpr int NA_INTEGER = std::numeric_limits<int>::quiet_NaN();

    inline bool is_na(double val) {
        return std::isnan(val);
    }
    
    class MedfateInternalError : public std::exception {
    private:
      std::string message;  
    public:
      MedfateInternalError(const std::string& msg) : message(msg) {}
      virtual const char* what() const noexcept override {
        return message.c_str();
      }
    };      

    
    
}

#endif