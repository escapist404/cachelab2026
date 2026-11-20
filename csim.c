// student id: 2025201731
// please change the above line to your student id

#include <assert.h>
#include <getopt.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void printSummary(int hits, int misses, int evictions) {
    printf("hits:%d misses:%d evictions:%d\n", hits, misses, evictions);
    FILE* output_fp = fopen(".csim_results", "w");
    assert(output_fp);
    fprintf(output_fp, "%d %d %d\n", hits, misses, evictions);
    fclose(output_fp);
}

void printHelp(const char* name) {
    printf(
        "Usage: %s [-hv] -s <num> -E <num> -b <num> -t <file>\n"
        "Options:\n"
        "  -h         Print this help message.\n"
        "  -v         Optional verbose flag.\n"
        "  -s <num>   Number of set index bits.\n"
        "  -E <num>   Number of lines per set.\n"
        "  -b <num>   Number of block offset bits.\n"
        "  -t <file>  Trace file.\n\n"
        "Examples:\n"
        "  linux>  %s -s 4 -E 1 -b 4 -t traces/yi.trace\n"
        "  linux>  %s -v -s 8 -E 2 -b 4 -t traces/yi.trace\n",
        name, name, name);
}

typedef struct CacheLine {
    bool valid;
    unsigned long long tag;
    unsigned long long lru;
} CacheLine;

int main(int argc, char* argv[]) {
    int s = -1, E = -1, b = -1;
    const char* trace_file = NULL;
    bool verbose = false;

    int opt;

    while ((opt = getopt(argc, argv, "hvs:E:b:t:")) != -1) {
        switch (opt) {
            case 'h':
                printHelp(argv[0]);
                return 0;
            case 'v':
                verbose = true;
                break;
            case 's':
                s = atoi(optarg);
                break;
            case 'E':
                E = atoi(optarg);
                break;
            case 'b':
                b = atoi(optarg);
                break;
            case 't':
                trace_file = optarg;
                break;
            default:
                printHelp(argv[0]);
                return 1;
        }
    }

    if (s < 0 || E <= 0 || b < 0 || trace_file == NULL) {
        fprintf(stderr, "Missing or invalid arguments.\n");
        printHelp(argv[0]);
        return 1;
    }
    (void)verbose;

    FILE* fp = fopen(trace_file, "r");
    if (fp == NULL) {
        perror("fopen");
        return 1;
    }

    int S = 1 << s;

    CacheLine* cache = calloc(S * E, sizeof(CacheLine));

    char line[256];
    char op;
    unsigned long long addr;

    unsigned long long hits = 0;
    unsigned long long misses = 0;
    unsigned long long evictions = 0;

    unsigned long long tick = 0;

    while (fgets(line, sizeof(line), fp) != NULL) {
        if (sscanf(line, " %c %llx", &op, &addr) == 2) {
            addr >>= b;
            unsigned long long set = addr & (S - 1);
            addr >>= s;

            bool ok = false;
            int ok_idx = -1;
            int min_lru_idx = -1;
            unsigned min_lru = -1;
            int empty_idx = -1;
            unsigned long long start = set * E;

            for (int i = 0; i < E; ++i) {
                if (cache[start + i].valid) {
                    if (cache[start + i].tag == addr) {
                        ok = true;
                        ok_idx = i;
                        break;
                    } else if (min_lru_idx == -1 || cache[start + i].lru < min_lru) {
                        min_lru_idx = i;
                        min_lru = cache[start + i].lru;
                    }
                } else if (empty_idx == -1) {
                    empty_idx = i;
                }
            }

            if (ok) {
                hits++;
                cache[start + ok_idx].lru = tick;
            } else {
                misses++;
                if (empty_idx == -1) {
                    evictions++;
                    cache[start + min_lru_idx].tag = addr;
                    cache[start + min_lru_idx].lru = tick;
                } else {
                    cache[start + empty_idx].valid = true;
                    cache[start + empty_idx].tag = addr;
                    cache[start + empty_idx].lru = tick;
                }
            }
            tick++;
        }
    }

    fclose(fp);
    printSummary(hits, misses, evictions);
    free(cache);
    return 0;
}
