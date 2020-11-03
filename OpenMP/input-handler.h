#ifndef INPUT_HANDLER_H
#define INPUT_HANDLER_H

#include <stdio.h>

typedef struct dynamic_graph dynamic_graph;
typedef struct dynamic_weighted_graph dynamic_weighted_graph;
typedef struct execution_settings execution_settings;
typedef struct phase_execution_briefing phase_execution_briefing;

#define WEIGHTED 1

typedef struct execution_settings
{
    char *input_file;
    int graph_type;
    double minimum_phase_improvement;
    double minimum_iteration_improvement;
    char *output_communities_file;
    char *output_graphs_file;
    int number_of_threads;
    int execution_settings_parallel_partitions_higher_power_of_2;
    int algorithm_version;

    int (*phase_executor_weighted)(dynamic_weighted_graph *, execution_settings *, dynamic_weighted_graph **, int **, phase_execution_briefing *);
    int (*phase_executor_not_weighted)(dynamic_graph *, execution_settings *, dynamic_weighted_graph **, int **, phase_execution_briefing *);

} execution_settings;

void set_default(execution_settings *s);
void print_help(char *prog_name);

int parse_args(int argc, char *argv[], execution_settings *s);

int dynamic_weighted_graph_parse_file(dynamic_weighted_graph *dg, char *filename);
int parse_input(dynamic_graph *dg, dynamic_weighted_graph *dwg, execution_settings *settings);
#endif
