//
// Copyright (c) 2012-2014 Institute of Industrial Science, The University of Tokyo.
// All right reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All right reserved.
//

#ifndef COMMAND_LINE_PARSER_H
#define COMMAND_LINE_PARSER_H

#include <sstream>
#include <string>
#include <vector>
#include <cstring>
#include <cassert>

class CLParser {

public:

  struct Param {
    const char* name;
    const char* value;
  };

  struct bad_value_type {
    const char* name;
    const char* value;
    bad_value_type(const char* name, const char* value)
      : name(name), value(value) {}
  };

public:

  CLParser(int nParams, Param params[]) : nParams(nParams), params(params) {}

  ~CLParser() {};

  bool parse(int nArg, const char* const args[]) {
    int nErrors = 0;
    for (int i = 1; i < nArg; i++) {
      Param* param = findParam(args[i]);
      if (param) {
        if (!parseArg(param, args[i])) nErrors++;
      } else {
        errors.push_back(std::string("*** Unknown parameter: ")
                         + std::string(args[i]));
        nErrors++;
      }
    }
    // if (nErrors == 0) dumpParams();
    return (nErrors == 0) ? true : false;
  }

  template<class T> T getValue(const char* name) const ;

  void printErrors() const {
    for (size_t i = 0; i < errors.size(); i++) {
      std::cerr << errors[i] << std::endl;
    }
  }

  bool checkParams() {
    int nErrors = 0;
    for (int i = 0; i < nParams; i++) {
      if (!params[i].value) {
        errors.push_back(std::string("*** Parameter not specified: ")
                         + std::string(params[i].name));
        nErrors++;
      }
    }
    return (nErrors == 0) ? true : false;
  }


private:

  int nParams;

  Param* params;

  std::vector<std::string> errors;

private:

  Param* findParam(const char* name) const {
    for (int i = 0; i < nParams; i++) {
      if (strncmp(name, params[i].name, strlen(params[i].name)) == 0) {
        return &params[i];
      }
    }
    return 0;
  }

  const char* getValueString(const char* name) const {
    Param* param = findParam(name);
    if (!param) return 0;
    return param->value;
  }

  bool parseArg(Param* param, const char* arg) {
    bool ret = true;
    static const char* TRUE = "true";
    const char* s = &arg[strlen(param->name)];
    if (*s == '\0') {
      param->value = TRUE;
    } else if (*s == '=' && s[1] != '\0') {
      param->value = &s[1];
    } else {
      errors.push_back(std::string("*** Bad parameter assignement: ")
                       + std::string(arg));
      return false;
    }
    return ret;
  }

  void dumpParams() {
    for (int i = 0; i < nParams; i++) {
      printf("[%s][%s]\n", params[i].name, params[i].value);
    }
  }

};


template<class T>
inline T CLParser::getValue(const char* name) const {
  T t;
  const char* s = getValueString(name);
  assert(s);
  std::istringstream ist(s);
  ist.exceptions(std::ios::failbit);
  try {
    ist >> t;
  }
  catch (const std::ios::failure& e) {
    throw bad_value_type(name, s); 
  }
  return t;
}


template<>
inline bool CLParser::getValue(const char* name) const {
  const char* s = getValueString(name);
  assert(s);
  if (strcasecmp(s, "yes") == 0 ||
      strcasecmp(s, "on") == 0 ||
      strcasecmp(s, "true") == 0) {
    return true;
  }
  else if (strcasecmp(s, "no") == 0 ||
           strcasecmp(s, "off") == 0 ||
           strcasecmp(s, "false") == 0) {
    return false;
  }
  else {
    throw bad_value_type(name, s); 
  }
}


template<>
inline std::string CLParser::getValue(const char* name) const {
  const char* s = getValueString(name);
  assert(s);
  return std::string(s);
}

#endif // COMMAND_LINE_PARSER_H
