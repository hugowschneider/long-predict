#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <zlib.h>
#include <math.h>
#include <curses.h>

#include "kseq.h"
#include "suffixArray.h"
#include "configFileParser.h"
#include "svm.h"


#define CODON_LENGTH 3

#define MAX_PATTERNS 340
const char possible_attributes[MAX_PATTERNS][5] = {
        "ol",
        "op",
        "ll",
        "lp",
        "aa",
        "ac",
        "ag",
        "at",
        "ca",
        "cc",
        "cg",
        "ct",
        "ga",
        "gc",
        "gg",
        "gt",
        "ta",
        "tc",
        "tg",
        "tt",
        "aaa",
        "aac",
        "aag",
        "aat",
        "aca",
        "acc",
        "acg",
        "act",
        "aga",
        "agc",
        "agg",
        "agt",
        "ata",
        "atc",
        "atg",
        "att",
        "caa",
        "cac",
        "cag",
        "cat",
        "cca",
        "ccc",
        "ccg",
        "cct",
        "cga",
        "cgc",
        "cgg",
        "cgt",
        "cta",
        "ctc",
        "ctg",
        "ctt",
        "gaa",
        "gac",
        "gag",
        "gat",
        "gca",
        "gcc",
        "gcg",
        "gct",
        "gga",
        "ggc",
        "ggg",
        "ggt",
        "gta",
        "gtc",
        "gtg",
        "gtt",
        "taa",
        "tac",
        "tag",
        "tat",
        "tca",
        "tcc",
        "tcg",
        "tct",
        "tga",
        "tgc",
        "tgg",
        "tgt",
        "tta",
        "ttc",
        "ttg",
        "ttt",
        "aaaa",
        "aaac",
        "aaag",
        "aaat",
        "aaca",
        "aacc",
        "aacg",
        "aact",
        "aaga",
        "aagc",
        "aagg",
        "aagt",
        "aata",
        "aatc",
        "aatg",
        "aatt",
        "acaa",
        "acac",
        "acag",
        "acat",
        "acca",
        "accc",
        "accg",
        "acct",
        "acga",
        "acgc",
        "acgg",
        "acgt",
        "acta",
        "actc",
        "actg",
        "actt",
        "agaa",
        "agac",
        "agag",
        "agat",
        "agca",
        "agcc",
        "agcg",
        "agct",
        "agga",
        "aggc",
        "aggg",
        "aggt",
        "agta",
        "agtc",
        "agtg",
        "agtt",
        "ataa",
        "atac",
        "atag",
        "atat",
        "atca",
        "atcc",
        "atcg",
        "atct",
        "atga",
        "atgc",
        "atgg",
        "atgt",
        "atta",
        "attc",
        "attg",
        "attt",
        "caaa",
        "caac",
        "caag",
        "caat",
        "caca",
        "cacc",
        "cacg",
        "cact",
        "caga",
        "cagc",
        "cagg",
        "cagt",
        "cata",
        "catc",
        "catg",
        "catt",
        "ccaa",
        "ccac",
        "ccag",
        "ccat",
        "ccca",
        "cccc",
        "cccg",
        "ccct",
        "ccga",
        "ccgc",
        "ccgg",
        "ccgt",
        "ccta",
        "cctc",
        "cctg",
        "cctt",
        "cgaa",
        "cgac",
        "cgag",
        "cgat",
        "cgca",
        "cgcc",
        "cgcg",
        "cgct",
        "cgga",
        "cggc",
        "cggg",
        "cggt",
        "cgta",
        "cgtc",
        "cgtg",
        "cgtt",
        "ctaa",
        "ctac",
        "ctag",
        "ctat",
        "ctca",
        "ctcc",
        "ctcg",
        "ctct",
        "ctga",
        "ctgc",
        "ctgg",
        "ctgt",
        "ctta",
        "cttc",
        "cttg",
        "cttt",
        "gaaa",
        "gaac",
        "gaag",
        "gaat",
        "gaca",
        "gacc",
        "gacg",
        "gact",
        "gaga",
        "gagc",
        "gagg",
        "gagt",
        "gata",
        "gatc",
        "gatg",
        "gatt",
        "gcaa",
        "gcac",
        "gcag",
        "gcat",
        "gcca",
        "gccc",
        "gccg",
        "gcct",
        "gcga",
        "gcgc",
        "gcgg",
        "gcgt",
        "gcta",
        "gctc",
        "gctg",
        "gctt",
        "ggaa",
        "ggac",
        "ggag",
        "ggat",
        "ggca",
        "ggcc",
        "ggcg",
        "ggct",
        "ggga",
        "gggc",
        "gggg",
        "gggt",
        "ggta",
        "ggtc",
        "ggtg",
        "ggtt",
        "gtaa",
        "gtac",
        "gtag",
        "gtat",
        "gtca",
        "gtcc",
        "gtcg",
        "gtct",
        "gtga",
        "gtgc",
        "gtgg",
        "gtgt",
        "gtta",
        "gttc",
        "gttg",
        "gttt",
        "taaa",
        "taac",
        "taag",
        "taat",
        "taca",
        "tacc",
        "tacg",
        "tact",
        "taga",
        "tagc",
        "tagg",
        "tagt",
        "tata",
        "tatc",
        "tatg",
        "tatt",
        "tcaa",
        "tcac",
        "tcag",
        "tcat",
        "tcca",
        "tccc",
        "tccg",
        "tcct",
        "tcga",
        "tcgc",
        "tcgg",
        "tcgt",
        "tcta",
        "tctc",
        "tctg",
        "tctt",
        "tgaa",
        "tgac",
        "tgag",
        "tgat",
        "tgca",
        "tgcc",
        "tgcg",
        "tgct",
        "tgga",
        "tggc",
        "tggg",
        "tggt",
        "tgta",
        "tgtc",
        "tgtg",
        "tgtt",
        "ttaa",
        "ttac",
        "ttag",
        "ttat",
        "ttca",
        "ttcc",
        "ttcg",
        "ttct",
        "ttga",
        "ttgc",
        "ttgg",
        "ttgt",
        "ttta",
        "tttc",
        "tttg",
        "tttt"
};


KSEQ_INIT(gzFile, gzread)


void usage(int argc, char *argv[]) {
    fprintf(stderr, "Usage: %s [-o] [-c <model config file>] -i <fasta file> [-d <both|-|+> ] [-s <size>]\n"
            "\t-c\tModel config file. File specifying the libSVM model and the used attributes\n"
            "\t-d\tDirection. The direction of the fasta file sequences for prediction. '+' will\n"
            "\t\tread the sequences as is. '-' will use the complementary sequence. And both will\n"
            "\t\tuse both options. This option is used for all sequences in the file. Default '+'\n"
            "\t-i\tInput Fasta file for prediction. This file should be a plain text or a gzip fasta\n"
            "\t\tfile\n"
            "\t-o\tData only. Output the frequencies of nucleotides patterns and orf size and relation to \n"
            "\t\ttranscript size.\n"
            "\t\tThe model file will determine the nucleotides patterns will be used, if not present, all\n"
            "\t\tpatterns will be calculated\n"
            "\t-s\tSize limit. This attribute ignore sequences shorter than this limit. Default 200\n"
            "\n"
            "Model Config File is a plain text file containing the following attributes:\n"
            "\tmodelFile\tThe path to the model file. It can be relative to the config file or\n"
            "\tdect\tThe model description\n"
            "\tattributes\tThe list of attributes used in the model training. This attributes are valid\n"
            "\t\t\t\tnucleotide frequencies, for example 'aa' and 'atc', and the values 'ol' for\n"
            "\t\t\t\tfirst ORF lenght and 'op' for first ORF percentage of the corresponding transcript\n"
            "\t\t\t\tlength\n"
            "Ex.:\n"
            "modelFile=human.model\n"
            "attributes=aa,aaa,ac,aca,acg,op\n"
            "", argv[0]);
    exit(EXIT_FAILURE);
}


size_t firstOrfSize(SuffixArray sa) {
    size_t first;
    size_t count;
    size_t length;
    size_t maxLenght;
    const char *suffix;
    size_t i;
    size_t j;
    count = suffixArraySearch(sa, "atg", &first);
    suffix = NULL;
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

size_t longestOrfSize(SuffixArray sa) {
    size_t first;
    size_t count;
    size_t maxLenght;
    const char *suffix;
    size_t i;
    size_t j;
    bool found;

    count = suffixArraySearch(sa, "atg", &first);


    maxLenght = 0;
    if (count) {
        for (i = first; i < first + count; i++) {
            suffix = sa->suffix[i];

            found = false;
            for (j = CODON_LENGTH; j + CODON_LENGTH < strlen(suffix); j += CODON_LENGTH) {
                if (strncmp(suffix + j, "taa", CODON_LENGTH) == 0 ||
                    strncmp(suffix + j, "tga", CODON_LENGTH) == 0 ||
                    strncmp(suffix + j, "tag", CODON_LENGTH) == 0) {
                    if (j + CODON_LENGTH > maxLenght) {
                        maxLenght = j + CODON_LENGTH;
                        found = true;
                        break;
                    }
                }
            }
            if (!found && strlen(suffix) > maxLenght) {
                maxLenght = strlen(suffix);
            }
        }
    }
    return maxLenght;
}

void compute(const char *seqName, const char *seq, const Config *config, FILE *output) {

    size_t count;
    size_t attr_length;
    size_t total;
    size_t length;
    ssize_t orfLength;
    ssize_t longestOrfLength;
    SuffixArray sa;

    length = strlen(seq);

    size_t attributeVectorSize;
    char **attributes;

    sa = suffixArrayCreate(seq);

    if (config) {
        attributeVectorSize = config->attributeVectorSize;
        attributes = config->attributes;
    } else {
        attributeVectorSize = MAX_PATTERNS;
        attributes = malloc(attributeVectorSize * sizeof(char *));
        size_t i;
        for (i = 0; i < attributeVectorSize; ++i) {
            attributes[i] = (char *) &possible_attributes[i];
        }

    }

    orfLength = -1;
    longestOrfLength = -1;

    int i, j;

    fprintf(output, "%s", seqName);
    double value = 0;
    for (i = 0; i < attributeVectorSize; ++i) {
        char *attr = attributes[i];

        if (strcasecmp(attr, "ol") == 0) {
            if (orfLength == -1) {
                orfLength = firstOrfSize(sa);
            }
            value = orfLength;

        } else if (strcasecmp(attr, "op") == 0) {
            if (orfLength == -1) {
                orfLength = firstOrfSize(sa);
            }
            value = (double) orfLength / (double) length;

        } else if (strcasecmp(attr, "ll") == 0) {
            if (longestOrfLength == -1) {
                longestOrfLength = longestOrfSize(sa);
            }
            value = longestOrfLength;

        } else if (strcasecmp(attr, "lp") == 0) {
            if (longestOrfLength == -1) {
                longestOrfLength = longestOrfSize(sa);
            }
            value = (double) longestOrfLength / (double) length;

        } else {
            count = suffixArraySearch(sa, attr, NULL);
            attr_length = strlen(attr);

            total = 0;
            for (j = 0; j < attr_length; ++j) {
                total += (length - j) / attr_length;
            }

            value = (double) count / (double) total;

        }
        fprintf(output, ",%.16f", value);
    }

    fprintf(output, "\n");
    fflush(output);

    suffixArrayDestroy(sa);

}

double predictData(const char *seqName, const char *seq, const Config *config, FILE *output,
                   struct svm_model *model) {

    size_t count;
    size_t attr_length;
    size_t total;
    size_t length;
    ssize_t orfLength;
    ssize_t longestOrfLength;
    SuffixArray sa;
    double predict_label;
    double predict_estimates;

    length = strlen(seq);

    sa = suffixArrayCreate(seq);

    orfLength = -1;
    longestOrfLength = -1;
    struct svm_node *x = malloc(sizeof(struct svm_node) * (config->attributeVectorSize + 1));
    int i, j;
    for (i = 0; i < config->attributeVectorSize; ++i) {
        char *attr = config->attributes[i];

        x[i].index = i + 1;
        if (strcasecmp(attr, "ol") == 0) {
            if (orfLength == -1) {
                orfLength = firstOrfSize(sa);
            }
            x[i].value = orfLength;

        } else if (strcasecmp(attr, "op") == 0) {
            if (orfLength == -1) {
                orfLength = firstOrfSize(sa);
            }
            x[i].value = (double) orfLength / (double) length;
        } else if (strcasecmp(attr, "ll") == 0) {
            if (longestOrfLength == -1) {
                longestOrfLength = longestOrfSize(sa);
            }
            x[i].value = longestOrfLength;

        } else if (strcasecmp(attr, "lp") == 0) {
            if (longestOrfLength == -1) {
                longestOrfLength = longestOrfSize(sa);
            }
            x[i].value = (double) longestOrfLength / (double) length;
        } else{
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
    for (i = 0; i < strlen(seq); ++i) {
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

char *join_strings(char *strings[], int count) {
    char *str = NULL;             /* Pointer to the joined strings  */
    size_t total_length = 0;      /* Total length of joined strings */
    size_t length = 0;            /* Length of a string             */
    int i = 0;                    /* Loop counter                   */

    /* Find total length of joined strings */
    for (i = 0; i < count; i++) {
        total_length += strlen(strings[i]);
        ++total_length; /* For newline to be added */
    }
    // ++total_length;     /* For joined string terminator */

    str = (char *) malloc(total_length);  /* Allocate memory for joined strings */
    str[0] = '\0';                      /* Empty string we can append to      */

    /* Append all the strings */
    for (i = 0; i < count; i++) {
        strcat(str, strings[i]);
        length = strlen(str);

        /* Check if we need to insert newline */
        if (i + 1 < count) {
            str[length] = ',';             /* Append a newline       */
        }
        str[length + 1] = '\0';           /* followed by terminator */

    }
    return str;
}

int main(int argc, char *argv[]) {

    int opt;
    int computeOnly = FALSE;

    char *confFile = NULL;
    char *fastaFile = NULL;
    char *direction = NULL;

    gzFile fp;
    kseq_t *seq;
    int l;
    ssize_t sizeLimit;

    FILE *directData = NULL;
    FILE *reverseData = NULL;

    char *directDataPath;
    char *complementaryDataPath;

    struct svm_model *model = NULL;

    sizeLimit = -1;
    while ((opt = getopt(argc, argv, "oc:i:d:")) != -1) {
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
            case 'o':
                computeOnly = TRUE;
                break;
            case 's':
                if (sscanf(optarg, "%zu", &sizeLimit) < 1 || sizeLimit < 0) {
                    sizeLimit = 200;
                }
                break;
            default:
                usage(argc, argv);
        }
    }

    if (!confFile && !computeOnly) {
        usage(argc, argv);
    }

    if (!fastaFile) {
        usage(argc, argv);
    }

    if (sizeLimit < 0) {
        sizeLimit = 200;
    }

    Config *config = NULL;
    if (confFile) {
        config = parseConfigFile(confFile);
    }
    if (!config && !computeOnly) {
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
        fprintf(stderr, "Cannot read fasta file '%s'.\n", fastaFile);
        free(directDataPath);
        free(complementaryDataPath);
        usage(argc, argv);
    }


    if (!computeOnly && access(config->modelFile, F_OK | R_OK) == -1) {
        fprintf(stderr, "Cannot read model file '%s'. Please review your config file\n", config->modelFile);
        free(directDataPath);
        free(complementaryDataPath);
        usage(argc, argv);
    }


    fp = gzopen(fastaFile, "r");
    seq = kseq_init(fp);

    if (config) {
        printf("Model file: %s\n", config->modelFile);
        if(config->desc) {
            printf("Model description: %s\n", config->desc);
        }

    }

    if (computeOnly) {
        size_t attributeVectorSize;
        char **attributes;
        if (config) {
            attributeVectorSize = config->attributeVectorSize;
            attributes = config->attributes;
        } else {
            attributeVectorSize = MAX_PATTERNS;

            attributes = malloc(attributeVectorSize * sizeof(char *));
            size_t i;
            for (i = 0; i < attributeVectorSize; ++i) {
                attributes[i] = (char *) &possible_attributes[i];
            }

        }

        char *string = join_strings(attributes, attributeVectorSize);

        if (!config) {
            free(attributes);
        }

        if (directData) {
            fprintf(directData, "ID,");
            fprintf(directData, "%s", string);
            fprintf(directData, "\n");
        }
        if (reverseData) {
            fprintf(reverseData, "ID,");
            fprintf(reverseData, "%s", string);
            fprintf(reverseData, "\n");
        }

    } else {
        model = svm_load_model(config->modelFile);
        if (directData) {
            fprintf(directData, "ID,Size,Classification,Probability lncRNA,Probability PCT\n");
        }
        if (reverseData) {
            fprintf(reverseData, "ID,Size,Classification,Probability lncRNA,Probability PCT\n");
        }
    }
    size_t count = 0;
    size_t reversePositive = 0;
    size_t positive = 0;
    if (computeOnly) {
        while ((l = kseq_read(seq)) >= 0) {
            if (strlen(seq->seq.s) >= sizeLimit) {
                if (directData) {
                    compute(seq->name.s, seq->seq.s, config, directData);

                }

                if (reverseData) {
                    reverse(seq->seq.s);
                    complementary(seq->seq.s);
                    compute(seq->name.s, seq->seq.s, config, reverseData);

                }
            }
        }
    } else {


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
    }
    if (config) {
        configDestroy(config);
    }
    kseq_destroy(seq);
    gzclose(fp);

    if (directData) {
        printf("Predicted lncRNAs: %.2f%% (%zu/%zu)\n", ((double) positive / (double) count) * 100.0, positive, count);
        fclose(directData);
    }

    if (reverseData && !computeOnly) {
        printf("Predicted lncRNAs in complementary sequence: %.2f%% (%zu/%zu)\n",
               ((double) reversePositive / (double) count) * 100.0, reversePositive, count);
        fclose(reverseData);
    }

    svm_free_and_destroy_model(&model);

    free(directDataPath);
    free(complementaryDataPath);


    return 0;
}


