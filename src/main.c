#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <zlib.h>
#include <math.h>

#include "kseq.h"
#include "suffixArray.h"
#include "configFileParser.h"
#include "svm.h"


#define CODON_LENGTH 3

KSEQ_INIT(gzFile, gzread)


void usage(int argc, char *argv[]) {
    fprintf(stderr, "Usage: %s -c <model config file> -i <fasta file> [-d <both|-|+> ] [-s <size>]\n"
            "\t-c\tModel config file. File specifying the libSVM model and the used attributes\n"
            "\t-i\tInput Fasta file for prediction\n"
            "\t-d\tDirection. The direction of the fasta file sequences for prediction. '+' will\n"
            "\t\tread the sequences as is. '-' will use the complementary sequence. And both will\n"
            "\t\tuse both options. This option is used for all sequences in the file. Default '+'\n"
            "\t-s\tSize limit. This attribute ignore sequences shorter than this limit. Default 200\n", argv[0]);
    exit(EXIT_FAILURE);
}

int fractionalDigits(double d) {
    char buf[
            1 + // sign, '-' or '+'
            (sizeof(d) * CHAR_BIT + 3) / 4 + // mantissa hex digits max
            1 + // decimal point, '.'
            1 + // mantissa-exponent separator, 'p'
            1 + // mantissa sign, '-' or '+'
            (sizeof(d) * CHAR_BIT + 2) / 3 + // exponent decimal digits max
            1 // string terminator, '\0'
    ];
    int n;
    char *pp, *p;
    int e, lsbFound, lsbPos;

    // convert d into "+/- 0x h.hhhh p +/- ddd" representation and check for errors
    if ((n = snprintf(buf, sizeof(buf), "%+a", d)) < 0 ||
        (unsigned) n >= sizeof(buf))
        return -1;

//printf("{%s}", buf);

    // make sure the conversion didn't produce something like "nan" or "inf"
    // instead of "+/- 0x h.hhhh p +/- ddd"
    if (strstr(buf, "0x") != buf + 1 ||
        (pp = strchr(buf, 'p')) == NULL)
        return 0;

    // extract the base-2 exponent manually, checking for overflows
    e = 0;
    p = pp + 1 + (pp[1] == '-' || pp[1] == '+'); // skip the exponent sign at first
    for (; *p != '\0'; p++) {
        if (e > INT_MAX / 10)
            return -2;
        e *= 10;
        if (e > INT_MAX - (*p - '0'))
            return -2;
        e += *p - '0';
    }
    if (pp[1] == '-') // apply the sign to the exponent
        e = -e;

//printf("[%s|%d]", buf, e);

    // find the position of the least significant non-zero bit
    lsbFound = lsbPos = 0;
    for (p = pp - 1; *p != 'x'; p--) {
        if (*p == '.')
            continue;
        if (!lsbFound) {
            int hdigit = (*p >= 'a') ? (*p - 'a' + 10) : (*p - '0'); // assuming ASCII chars
            if (hdigit) {
                static const int lsbPosInNibble[16] = {0, 4, 3, 4, 2, 4, 3, 4, 1, 4, 3, 4, 2, 4, 3, 4};
                lsbFound = 1;
                lsbPos = -lsbPosInNibble[hdigit];
            }
        } else {
            lsbPos -= 4;
        }
    }
    lsbPos += 4;

    if (!lsbFound)
        return 0; // d is 0 (integer)

    // adjust the least significant non-zero bit position
    // by the base-2 exponent (just add them), checking
    // for overflows

    if (lsbPos >= 0 && e >= 0)
        return 0; // lsbPos + e >= 0, d is integer

    if (lsbPos < 0 && e < 0)
        if (lsbPos < INT_MIN - e)
            return -2; // d isn't integer and needs too many fractional digits

    if ((lsbPos += e) >= 0)
        return 0; // d is integer

    if (lsbPos == INT_MIN && -INT_MAX != INT_MIN)
        return -2; // d isn't integer and needs too many fractional digits

    return -lsbPos;
}

int compare(const void *a, const void *b) {
    size_t size_t_a = *((size_t *) a);
    size_t size_t_b = *((size_t *) b);

    if (size_t_a == size_t_b) return 0;
    else if (size_t_a < size_t_b) return -1;
    else return 1;
}


size_t firstOrfSize(SuffixArray sa, const char *seq) {
    size_t first;
    size_t count;
    size_t length;
    size_t maxLenght;
    const char *suffix;
    size_t i;
    size_t j;
    count = suffixArraySearch(sa, "atg", &first);

    if (count) {
        maxLenght = 0;
        for (i = first; i < first + count; i++) {
            length = strlen(sa->suffix[i]);
            if (maxLenght < length) {
                maxLenght = length;
                suffix = sa->suffix[i];
            }
        }
        for (j = CODON_LENGTH; j + CODON_LENGTH < strlen(suffix); j += CODON_LENGTH) {
            if (strncmp(suffix + j, "taa", CODON_LENGTH) == 0 ||
                strncmp(suffix + j, "tga", CODON_LENGTH) == 0 ||
                strncmp(suffix + j, "tag", CODON_LENGTH) == 0) {
                return j + CODON_LENGTH;
            }
        }
        return strlen(suffix);


    }
    return 0;
}

double predictData(const char *seqName, const char *seq, const Config *config, FILE *output,
                   struct svm_model *model) {

    size_t count;
    size_t attr_length;
    size_t total;
    size_t length;
    ssize_t orfLength;
    SuffixArray sa;
    double predict_label;
    double predict_estimates;

    length = strlen(seq);

    sa = suffixArrayCreate(seq);

    orfLength = -1;
    struct svm_node *x = malloc(sizeof(struct svm_node) * (config->attributeVectorSize + 1));
    size_t i, j;
    for (i = 0; i < config->attributeVectorSize; ++i) {
        char *attr = config->attributes[i];

        x[i].index = i + 1;
        if (strcasecmp(attr, "ol") == 0) {
            if (orfLength == -1) {
                orfLength = firstOrfSize(sa, seq);
            }
            x[i].value = orfLength;

        } else if (strcasecmp(attr, "op") == 0) {
            if (orfLength == -1) {
                orfLength = firstOrfSize(sa, seq);
            }
            x[i].value = (double) orfLength / (double) length;
        } else {
            count = suffixArraySearch(sa, attr, NULL);
            attr_length = strlen(attr);

            total = 0;
            for (j = 0; j < attr_length; ++j) {
                total += (length - j) / attr_length;
            }

            x[i].value = (double) count / (double) total;

        }
    }
    x[i].index = -1;
//    fprintf(output, ">%s", seqName);
//    for (int i = 0; x[i].index != -1; ++i) {
//        fprintf(output, " %d:%.16f", i + 1, x[i].value);
//
//    }
//    fprintf(output, "\n");
//    fprintf(output, "%s\n", seq);

    predict_label = svm_predict_probability(model, x, &predict_estimates);
    fprintf(output, "\"%s\",%zu,%s,%.16f,%.16f\n", seqName, length,
            predict_label > 0 ? "lncRNA" : "PCT", predict_estimates, 1 - predict_estimates);


    free(x);
    suffixArrayDestroy(sa);

    return predict_label;


}

void reverse(char *str) {
    if (str) {
        char *end = str + strlen(str) - 1;

        // swap the values in the two given variables
        // XXX: fails when a and b refer to same memory location
#   define XOR_SWAP(a, b) do\
    {\
      a ^= b;\
      b ^= a;\
      a ^= b;\
    } while (0)

        // walk inwards from both ends of the string,
        // swapping until we get to the middle
        while (str < end) {
            XOR_SWAP(*str, *end);
            str++;
            end--;
        }
#   undef XOR_SWAP
    }
}

void complementary(char *seq) {
    size_t i;
    for ( i = 0; i < strlen(seq); ++i) {
        switch (seq[i]) {
            case 'a':
                seq[i] = 't';
                break;
            case 't':
                seq[i] = 'a';
                break;
            case 'g':
                seq[i] = 'c';
                break;
            case 'c':
                seq[i] = 'g';
                break;
            default:
                break;

        }

    }
}

int main(int argc, char *argv[]) {

    int opt;

    char *confFile = NULL;
    char *fastaFile = NULL;
    char *direction = NULL;

    char *querySequence;

    SuffixArray sa;

    gzFile fp;
    kseq_t *seq;
    int l;
    ssize_t sizeLimit;

    FILE *directData = NULL;
    FILE *reverseData = NULL;

    char *directDataPath;
    char *complementaryDataPath;

    sizeLimit = -1;
    while ((opt = getopt(argc, argv, "c:i:d:")) != -1) {
        switch (opt) {
            case 'c':
                confFile = optarg;
                break;
            case 'i':
                fastaFile = optarg;
                break;
            case 'd':
                direction = optarg;
                break;
            case 's':
                if (sscanf(optarg, "%zd", &sizeLimit) < 1 || sizeLimit < 0) {
                    sizeLimit = 200;
                }
                break;
            default:
                usage(argc, argv);
        }
    }

    if (!confFile) {
        usage(argc, argv);
    }

    if (!fastaFile) {
        usage(argc, argv);
    }

    if (sizeLimit < 0) {
        sizeLimit = 200;
    }

    Config *config = parseConfigFile(confFile);
    if (!config) {
        usage(argc, argv);
    }

    directDataPath = malloc(sizeof(char) * (strlen(fastaFile) + 10));
    strcpy(directDataPath, fastaFile);
    strcat(directDataPath, ".pred.csv");
    complementaryDataPath = malloc(sizeof(char) * (strlen(fastaFile) + 18));
    strcpy(complementaryDataPath, fastaFile);
    strcat(complementaryDataPath, ".reverse.pred.csv");

    if (direction) {
        if (strcasecmp(direction, "both") == 0) {
            directData = fopen(directDataPath, "w");
            reverseData = fopen(complementaryDataPath, "w");
        } else if (strcasecmp(direction, "-") == 0) {
            reverseData = fopen(complementaryDataPath, "w");
        } else if (strcasecmp(direction, "+") == 0) {
            directData = fopen(directDataPath, "w");
        } else {
            usage(argc, argv);
        }
    } else {
        directData = fopen(directDataPath, "w");
    }

    if (access(fastaFile, F_OK | R_OK) == -1) {
        printf("Cannot read fasta file '%s'.\n", fastaFile);
        usage(argc, argv);
    }

    fp = gzopen(fastaFile, "r");
    seq = kseq_init(fp);

    struct svm_model *model = svm_load_model(config->modelFile);
    if (directData) {
        fprintf(directData, "ID,Size,Classification,Probability lncRNA,Probability PCT\n");
    }
    if (reverseData) {
        fprintf(reverseData, "ID,Size,Classification,Probability lncRNA,Probability PCT\n");
    }
    size_t count = 0;
    size_t reversePositive = 0;
    size_t positive = 0;
    while ((l = kseq_read(seq)) >= 0) {
        if (strlen(seq->seq.s) >= sizeLimit) {
            count++;
            if (directData) {
                if (predictData(seq->name.s, seq->seq.s, config, directData, model) > 0) {
                    positive++;
                }
            }

            if (reverseData) {
                reverse(seq->seq.s);
                complementary(seq->seq.s);
                if (predictData(seq->name.s, seq->seq.s, config, reverseData, model)) {
                    reversePositive++;
                }
            }
        }
    }

    configDestroy(config);
    kseq_destroy(seq);
    gzclose(fp);

    if (directData) {
        printf("Predicted lncRNAs: %.2f%% (%zu/%zu)\n", ((double) positive / (double) count) * 100.0, positive, count);
        fclose(directData);
    }

    if (reverseData) {
        printf("Predicted lncRNAs in complementary sequence: %.2f%% (%zu/%zu)\n",
               ((double) reversePositive / (double) count) * 100.0, reversePositive, count);
        fclose(reverseData);
    }

    svm_free_and_destroy_model(&model);

    free(directDataPath);
    free(complementaryDataPath);


    return 0;
}


