//
// Created by Hugo Schneider on 8/3/16.
//

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <libgen.h>

#include "configFileParser.h"

#define BUFFER_SIZE 2048

size_t trimwhitespace(char *out, size_t len, const char *str)
{
    if(len == 0)
        return 0;

    const char *end;
    size_t out_size;

    // Trim leading space
    while(isspace(*str)) str++;

    if(*str == 0)  // All spaces?
    {
        *out = 0;
        return 1;
    }

    // Trim trailing space
    end = str + strlen(str) - 1;
    while(end > str && isspace(*end)) end--;
    end++;

    // Set output size to minimum of trimmed string length and buffer size minus 1
    out_size = (end - str) < len-1 ? (end - str) : len-1;

    // Copy trimmed string and add null terminator
    memcpy(out, str, out_size);
    out[out_size] = 0;

    return out_size;
}

Config *parseConfigFile(const char *configFilePath) {

    Config *config = 0;
    FILE *configFile;

    char *buffer;
    char *propBuffer;
    char *valueBuffer;
    size_t bufferSize = BUFFER_SIZE;
    ssize_t charactersRead;
    size_t tokenCount;
    size_t tokenSize;
    char * realPath = realpath(configFilePath, NULL);
    char * dirPath = dirname(realPath);

    char *token;
    char **attributes;

    buffer = (char *) malloc(bufferSize * sizeof(char));
    propBuffer = (char *) malloc(bufferSize * sizeof(char));
    valueBuffer = (char *) malloc(bufferSize * sizeof(char));
    if (buffer == NULL) {
        perror("Unable to allocate buffer");
        return 0;
    }
    config = (Config *) malloc(sizeof(Config));
    config->desc = NULL;
    configFile = fopen(configFilePath, "r");

    if (!configFile) {
        perror("Unable to open model cofig file");
        return 0;
    }

    attributes = (char **) malloc(BUFFER_SIZE * sizeof(char *));

    while ((charactersRead = getline(&buffer, &bufferSize, configFile)) != EOF) {
        if (charactersRead > 0) {
            if (sscanf(buffer, "%[abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ]=%[^\t\n]", propBuffer, valueBuffer) ==
                2) {
                if (strcasecmp(propBuffer, "modelFile") == 0) {
                    config->modelFile = malloc(strlen(valueBuffer) + strlen(dirPath) + 2);
                    trimwhitespace(config->modelFile, strlen(valueBuffer) + 1, valueBuffer);
                    strcpy(valueBuffer, config->modelFile);

                    if(valueBuffer[0] == '/'){
                        strcpy(config->modelFile, valueBuffer);
                    } else {
                        strcpy(config->modelFile, dirPath);
                        strcat(config->modelFile, "/");
                        strcat(config->modelFile, valueBuffer);
                    }

                } else if (strcasecmp(propBuffer, "attributes") == 0) {
                    tokenCount = 0;

                    token = strtok(valueBuffer, ",");
                    while (token != NULL && (tokenSize = strlen(token)) > 0) {

                        attributes[tokenCount] = (char *) malloc(tokenSize + 1);
                        trimwhitespace(attributes[tokenCount], tokenSize + 1, token);
                        tokenCount++;
                        token = strtok(NULL, " ,.-");
                    }
                    config->attributeVectorSize = tokenCount;
                    config->attributes = malloc(tokenCount * sizeof(char *));
                    memcpy(config->attributes, attributes, tokenCount * sizeof(char *));
                } else if (strcasecmp(propBuffer, "desc") == 0) {
                    config->desc = malloc(sizeof(char) * strlen(valueBuffer) + 1);
                    strcpy(config->desc, valueBuffer);
                }
            }
        }
    }

    free(attributes);
    free(buffer);
    free(propBuffer);
    free(valueBuffer);


    return config;
}

void configDestroy(Config *config) {
    size_t i;
    for ( i = 0; i < config->attributeVectorSize; ++i) {
        free(config->attributes[i]);
    }
    free(config->desc);
    free(config->attributes);
    free(config->modelFile);
    free(config);
}
