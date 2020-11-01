#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parse-args.h"
#include "execution-settings.h"
#include "community-development.h"
#include "utilities.h"
#include "version-parallel-sort-select-chunks.h"


void set_default(execution_settings *s) {
	s->input_file = NULL;
	s->graph_type = WEIGHTED;
	s->minimum_phase_improvement = DEFAULT_MINIMUM_PHASE_IMPROVEMENT;
	s->minimum_iteration_improvement = DEFAULT_MINIMUM_ITERATION_IMPROVEMENT;
	s->output_communities_file = NULL;
	s->output_graphs_file = NULL;
	s->number_of_threads = DEFAULT_NUMBER_OF_THREADS;
	s->benchmark_runs = DEFAULT_BENCHMARK_RUNS;
	s->verbose = 0;
	s->input_file_format = DEFAULT_FILE_FORMAT;
	s->algorithm_version = DEFAULT_ALGORITHM_VERSION;
	s->execution_settings_sort_select_chunks_chunk_size = DEFINE_EXECUTION_SETTINGS_SORT_SELECT_CHUNKS_CHUNK_SIZE;
	s->execution_settings_vertex_following = DEFAULT_EXECUTION_SETTINGS_VERTEX_FOLLOWING;

	s->execution_settings_parallel_partitions_higher_power_of_2 = DEFAULT_EXECUTION_SETTINGS_PARALLEL_PARTITIONS_HIGHER_POWER_OF_2;
}

void print_help(char *prog_name) {

	printf("\nUsage:\n"
			"%s input-file\n", prog_name);

	printf(PRINTING_UTILITY_STARS);
}

int parse_args(int argc, char *argv[], execution_settings *s){
	int valid = 1;
	int i, execution_option, file_format, algorithm_version;

	set_default(s);

	if(argc < 2) {
		printf("Wrong number of arguments!\n");

		valid = 0;
	}

	if(valid && strcmp(argv[1],"-h") == 0) {
		valid = 0;
	}

	if(valid) {

		s->input_file = argv[1];
		s->algorithm_version = ALGORITHM_VERSION_PARALLEL_1_SORT_SELECT;
		if(i + 1 < argc) {
			i++;
			file_format = atoi(argv[i]);
		}
		s->input_file_format = FILE_FORMAT_EDGE_LIST_WEIGHTED;
	}

	if(!valid)
		print_help(argv[0]);

	return valid;
}
