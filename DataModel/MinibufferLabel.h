// Enum class that represents the type of trigger that led to the
// creation of a particular minibuffer
//#pragma once

#ifndef MINIBUFFERLABEL_H
#define MINIBUFFERLABEL_H

#include <fstream>
#include <sstream>
#include <string>

// TODO: update this to allow for all known minibuffer labels
enum class MinibufferLabel : uint8_t { Unknown = 0u, LED, Soft, Beam, Cosmic,
  Source, Hefty, HeftySource, MRDCR, AmBe };

// Operator for printing a description of the minibuffer label to a
// std::ostream
inline std::ostream& operator<<(std::ostream& out, const MinibufferLabel& mbl)
{
  switch (mbl) {
    case MinibufferLabel::LED:
      out << "LED";
      break;
    case MinibufferLabel::Soft:
      out << "Soft";
      break;
    case MinibufferLabel::Beam:
      out << "Beam";
      break;
    case MinibufferLabel::Cosmic:
      out << "Cosmic";
      break;
    case MinibufferLabel::Source:
      out << "Source";
      break;
    case MinibufferLabel::Hefty:
      out << "Hefty";
      break;
    case MinibufferLabel::HeftySource:
      out << "HeftySource";
      break;
    case MinibufferLabel::AmBe:
      out << "AmBe";
      break;
    case MinibufferLabel::MRDCR:
      out << "MRDCR";
      break;
    default:
      out << "Unknown";
      break;
  }
  return out;
}

inline std::string minibuffer_label_to_string(const MinibufferLabel& mbl){
  std::stringstream temp_stream;
  temp_stream << mbl;
  return temp_stream.str();
}

#endif
