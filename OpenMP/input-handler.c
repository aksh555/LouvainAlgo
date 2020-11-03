#include "dynamic-weighted-graph.h"
#include "input-handler.h"
#include "utilities.h"
#include <stdio.h>
#include <string.h>

int dynamic_weighted_graph_parse_file(dynamic_weighted_graph *dg, char *filename)
{
	FILE *graph_file;
	int src, dest, weight;

	if (!dynamic_weighted_graph_init(dg, DEFAULT_INIT_NODES_SIZE))
		return 0;

	if (graph_file = fopen(filename, "r"))
	{
		while (fscanf(graph_file, "%d %d %d", &src, &dest, &weight) == 3)
			if (!dynamic_weighted_graph_insert(dg, src, dest, weight))
			{
				printf("Could not insert edge %d - %d!\n", src, dest);
				return 0;
			}
	}
	else
	{
		printf("Could not open file: %s!\n", filename);

		return 0;
	}

	// Reduce to optimal size
	if (!dynamic_weighted_graph_reduce(dg))
	{
		printf("Could not reduce graph size!\n");
		return 0;
	}
	return 1;
}

int parse_input(dynamic_graph *dg, dynamic_weighted_graph *dwg, execution_settings *settings)
{
	int valid = 1;
	settings->graph_type = WEIGHTED;

	if (!dynamic_weighted_graph_parse_file(dwg, settings->input_file))
		valid = 0;

	if (!valid)
	{
		printf("Failed to read input file: %s\n", settings->input_file);
		return 0;
	}
	return 1;
}

void set_default(execution_settings *s)
{
	s->input_file = NULL;
	s->graph_type = WEIGHTED;
	s->minimum_phase_improvement = 0;
	s->minimum_iteration_improvement = 0;
	s->output_communities_file = NULL;
	s->output_graphs_file = NULL;
	s->number_of_threads = 4;
	s->algorithm_version = 1;
	s->execution_settings_parallel_partitions_higher_power_of_2 = 0;
}

void print_help(char *prog_name)
{

	printf("\nUsage:\n"
		   "%s input-file\n",
		   prog_name);

	printf(PRINTING_UTILITY_STARS);
}

int parse_args(int argc, char *argv[], execution_settings *s)
{
	int valid = 1;
	int i, execution_option, file_format, algorithm_version;

	set_default(s);

	if (argc < 2)
	{
		printf("Wrong number of arguments!\n");

		valid = 0;
	}

	if (valid && strcmp(argv[1], "-h") == 0)
	{
		valid = 0;
	}

	if (valid)
	{

		s->input_file = argv[1];
		s->algorithm_version = 1;
		if (i + 1 < argc)
		{
			i++;
			file_format = atoi(argv[i]);
		}
	}

	if (!valid)
		print_help(argv[0]);

	return valid;
}
