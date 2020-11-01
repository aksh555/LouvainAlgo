#include "parse-args.h"
#include "execution-settings.h"
#include "dynamic-graph.h"
#include "dynamic-weighted-graph.h"
#include "input-handler.h"
#include "execution-handler.h"
#include "utilities.h"
#include "community-development.h"
#include "execution-briefing.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char * argv[]) {
	execution_settings settings;
	dynamic_graph input_dg;
	dynamic_weighted_graph input_dwg;
	execution_briefing briefing;

	dynamic_weighted_graph *output_dwg;
	int *output_communities;

	if(!parse_args(argc, argv, &settings))
		return -1;

	if(!parse_input(&input_dg, &input_dwg, &settings))
		return -1;

	if(!execute_community_detection(&input_dg, &input_dwg, &settings, &output_dwg, &output_communities, &briefing))
		printf("Community detection exited with errors!\n");
}

