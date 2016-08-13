//
// Created by Hugo Schneider on 8/3/16.
//

#ifndef LONG_PREDICT_CONFIGFILEPARSER_H
#define LONG_PREDICT_CONFIGFILEPARSER_H

#include <stdlib.h>

typedef struct {
    char * modelFile;
    char ** attributes;
    size_t attributeVectorSize;
} Config;


Config * parseConfigFile(char * configFilePath);

void configDestroy(Config * config);


#endif //LONG_PREDICT_CONFIGFILEPARSER_H
