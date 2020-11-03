#include "dynamic-weighted-graph.h"
#include "input-handler.h"
#include "execution-handler.h"
#include "utilities.h"
#include "community-development.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
	execution_settings settings;
	dynamic_graph input_dg;
	dynamic_weighted_graph input_dwg;
	execution_briefing briefing;

	dynamic_weighted_graph *output_dwg;
	int *output_communities;

	if (!parse_args(argc, argv, &settings))
		return -1;

	if (!parse_input(&input_dg, &input_dwg, &settings))
		return -1;

	if (!execute_community_detection(&input_dg, &input_dwg, &settings, &output_dwg, &output_communities, &briefing))
		printf("Community detection exited with errors!\n");

	printf("\nExecution Summary\n\n");
	printf("\tAverage modularity obtained:                         %f\n", briefing.output_modularity);
	printf("\tAverage execution time:                              %fs\n", briefing.execution_time);
	printf("\tAverage sum of execution time over all threads:      %fs\n", briefing.clock_execution_time);
	printf("\tTotal execution time:                                %fs\n", briefing.global_execution_time);
}
