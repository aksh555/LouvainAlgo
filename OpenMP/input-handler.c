#include "parse-args.h"
#include "execution-settings.h"
#include "dynamic-graph.h"
#include "dynamic-weighted-graph.h"
#include "utilities.h"
#include <stdio.h>
#include <string.h>
#include "shared-graph.h"

int dynamic_weighted_graph_parse_file(dynamic_weighted_graph *dg, char *filename) {
	FILE * graph_file;
	int src, dest, weight;

	if(!dynamic_weighted_graph_init(dg,DEFAULT_INIT_NODES_SIZE))
		return 0;

	// Read file - Insert edges
	if(graph_file = fopen(filename,"r")) {
		while(fscanf(graph_file, "%d %d %d", &src, &dest, &weight) == 3)
			if(!dynamic_weighted_graph_insert(dg,src,dest, weight)) {
				printf("Could not insert edge %d - %d!\n", src, dest);
				return 0;
			}
	} else {
		printf("Could not open file: %s!\n", filename);

		return 0;
	}

	// Reduce to optimal size
	if(!dynamic_weighted_graph_reduce(dg)) {
		printf("Could not reduce graph size!\n");
		return 0;
	}
	return 1;
}

int parse_input(dynamic_graph *dg, dynamic_weighted_graph *dwg, execution_settings *settings) {
	int valid = 1;
	settings->graph_type = WEIGHTED;

	if(!dynamic_weighted_graph_parse_file(dwg, settings->input_file))
			valid = 0;

	if(!valid){
		printf("Failed to read input file: %s\n", settings->input_file);
		return 0;
	}
	return 1;
}

