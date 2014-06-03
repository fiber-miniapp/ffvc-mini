//
// Copyright (c) 2012-2014 Institute of Industrial Science, The University of Tokyo.
// All right reserved.
//
// Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
// All right reserved.
//

#ifndef COMMAND_LINE_H
#define COMMAND_LINE_H

#define USE_MPI

#include <iostream>
#include <cstdlib>
#include "CLParser.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

class CommandLine {

public:

  std::string scale;
  int size;
  int division[3];

  int step;
  double dt;

  std::string comm_mode;

  int p_itr, vp_itr;
  bool practical;

  int output_interval;

public:

  CommandLine() {
    static const char* usage = 
      "usage: ffvc_mini --scale=str --size=int [options]\n"
      "    --scale=str            strong or weak\n"
      "    --size=int             system size\n"
      "options:\n"
      "    --division=str         MPI node division (LxMxN) [1x1x1]\n"
      "    --step=int             time steps [20]\n"
      "    --dt=double            time delta [0.2]\n"
      "    --comm=str             sync or async [async]\n"
      "    --p_itr=int            max number of Poisson iteration [30]\n"
      "    --vp_itr=int           max number of VP iteration [20]\n"
      "    --practical            run in practical calculation mode\n"
      "    --output_interval=int  output interval [0]\n"
      "    --help                 print this message\n";

    static CLParser::Param params[] = {
      { "--scale", 0 },
      { "--size", 0 },
      { "--division", "1x1x1" },
      { "--step", "20" },
      { "--dt", "0.2" },
      { "--comm", "async" },
      { "--p_itr", "30" },
      { "--vp_itr", "20" },
      { "--practical", "false" },
      { "--output_interval", "0" },
      { "--help", "false" },
    };

    int nParams = sizeof(params) / sizeof(CLParser::Param);
    parser = new CLParser(nParams, params);

    usageString = usage;
  };


  ~CommandLine() {
    delete parser;
  };


#ifdef USE_MPI
  void parse(int argc, const char *const argv[], MPI_Comm comm=MPI_COMM_WORLD) {
    int myRank;
    MPI_Comm_rank(comm, &myRank);
#else
  void parse(int argc, const char *const argv[]) {
    int myRank = 0;
#endif

    if (!parser->parse(argc, argv)) {
      if (myRank == 0) parser->printErrors();
      if (myRank == 0) showUsage();
      exitProg(1);
    }

    int nError = 0;

    try {
      bool help = parser->getValue<bool>("--help");
      if (help) {
        if (myRank == 0) showUsage();
        exitProg(0);
      }

      if (!parser->checkParams()) {
        if (myRank == 0) parser->printErrors();
        if (myRank == 0) showUsage();
        exitProg(1);
      }

      scale = parser->getValue<std::string>("--scale");
      if (scale != "strong" && scale != "weak") {
        if (myRank == 0) std::cerr << "*** scale must be strong or weak: " << scale << std::endl;
        nError++;
      }

      size = parser->getValue<int>("--size");
      if (size <= 0) {
        if (myRank == 0) std::cerr << "*** size must be positive: " << size << std::endl;
        nError++;
      }
      if (size%2 != 0) {
        if (myRank == 0) std::cerr << "*** size must be even: " << size << std::endl;
        nError++;
      }

      std::string div = parser->getValue<std::string>("--division");
      if (sscanf(div.c_str(), "%dx%dx%d",
                 &division[0], &division[1], &division[2]) != 3) {
        if (myRank == 0) std::cerr << "*** division must be NxNxN form: " << div << std::endl;
        nError++;
      }
      else if (division[0] <= 0 || division[1] <= 0 || division[2] <= 0) {
        if (myRank == 0) std::cerr << "*** division number must be positive: "  << div << std::endl;
        nError++;
      }

      step = parser->getValue<int>("--step");
      if (step <= 0) {
        if (myRank == 0) std::cerr << "*** step must be positive: " << step << std::endl;
        nError++;
      }

      dt = parser->getValue<double>("--dt");
      if (dt <= 0.0) {
        if (myRank == 0) std::cerr << "*** dt must be positive: " << dt << std::endl;
        nError++;
      }

      comm_mode = parser->getValue<std::string>("--comm");
      if (comm_mode != "sync" && comm_mode != "async") {
        if (myRank == 0) std::cerr << "*** comm must be sync or async: " << comm_mode << std::endl;
        nError++;
      }

      p_itr = parser->getValue<int>("--p_itr");
      if (p_itr <= 0.0) {
        if (myRank == 0) std::cerr << "*** p_itr must be positive: " << p_itr << std::endl;
        nError++;
      }

      vp_itr = parser->getValue<int>("--vp_itr");
      if (vp_itr <= 0.0) {
        if (myRank == 0) std::cerr << "*** vp_itr must be positive: " << vp_itr << std::endl;
        nError++;
      }

      practical = parser->getValue<bool>("--practical");

      output_interval = parser->getValue<int>("--output_interval");

    }
    catch (const CLParser::bad_value_type e) {
      std::cerr << "*** bad value type: "
                << e.name << "=" << e.value << std::endl;
      nError++;
    }

    if (nError > 0) {
      if (myRank == 0) showUsage();
      exitProg(1);
    }
  }

private:

  CLParser* parser;

  const char* usageString;

private:

  void showUsage() {
    std::cerr << std::endl << usageString << std::endl;
  }


  void exitProg(int code = 0) {
#ifdef USE_MPI
    MPI_Finalize();
#endif
    exit(code);
  }

};

#endif // COMMAND_LINE_H
