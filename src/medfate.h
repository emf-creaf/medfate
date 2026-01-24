#include <vector>
#include <string>
#include <limits>
#include <cmath>

#ifndef MEDFATE_H
#define MEDFATE_H

// ============================================================================
// PURE C++ STRUCTURES (Internal use)
// ============================================================================

namespace medfate {
    constexpr double NA_DOUBLE = std::numeric_limits<double>::quiet_NaN();

    inline bool is_na(double val) {
        return std::isnan(val);
    }

    template<typename T>
    using Matrix = std::vector<std::vector<T>>;
    
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