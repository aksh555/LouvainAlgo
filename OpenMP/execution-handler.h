#ifndef EXECUTION_HANDLER_H
#define EXECUTION_HANDLER_H

typedef struct execution_briefing {
	int performed_runs;
	int execution_successful;
	double output_modularity;
	double execution_time;
	double clock_execution_time;
	// Precompute time is NOT included in execution times
	double precompute_time;

	double minimum_execution_time;
	double minimum_clock_execution_time;

	// Includes file IO, not relevant for scalability measures!
	double global_execution_time;
} execution_briefing;

typedef struct algorithm_execution_briefing {
	int execution_successful;
	int number_of_phases;
	double output_modularity;
	double execution_time;
	double clock_execution_time;
	// Precompute time is NOT included in execution times
	double precompute_time;

	// Includes file IO, not relevant for scalability measures!
	double global_execution_time;
} algorithm_execution_briefing;

typedef struct phase_execution_briefing {
	int execution_successful;
	int number_of_iterations;
	double average_iteration_duration;
	double output_modularity;
} phase_execution_briefing;

int execute_community_detection(dynamic_graph *input_dg, dynamic_weighted_graph *input_dwg, execution_settings *settings, dynamic_weighted_graph **community_graph, int **community_vector, execution_briefing *briefing);
#endif
