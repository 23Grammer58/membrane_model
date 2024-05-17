//
// Created by alex on 21.09.2020.
//

#ifndef AORTIC_VALVE_BENCHMARKS_H
#define AORTIC_VALVE_BENCHMARKS_H

#include "BenchCommon.h"

int Benchmark1(int argc, char* argv[]);
int Benchmark2(int argc, char* argv[]);
int OneAxisStretch(int argc, char* argv[]);
int Benchmark4(int argc, char* argv[]);
int BenchKyriacou(int argc, char **argv);
int BenchmarkBendAnnulus(int argc, char* argv[]);
int BenchThreeLeaflet(int argc, char* argv[]);
int processThreeLeaflet(int argc, char* argv[]);
int BenchDataDriven(int argc, char* argv[]);
int testHalfSphereStatDef(int argc, char* argv[]);
int VanLoonBench(int argc, char* argv[]);
int LinVsNonlin(int argc, char* argv[]);

int ContactTest(int argc, char* argv[]);

int ElasticBenchmark(int argc, char* argv[]);



#endif //AORTIC_VALVE_BENCHMARKS_H
