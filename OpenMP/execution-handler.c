#include "dynamic-weighted-graph.h"
#include "input-handler.h"
#include "utilities.h"
#include "community-computation-weighted.h"
#include "execution-handler.h"
#include <stdio.h>
#include <omp.h>
#include <time.h>
#include <stdlib.h>

void merge_briefings_average(execution_briefing *briefing, algorithm_execution_briefing *internal_briefing, int number_of_previous_runs)
{
	briefing->clock_execution_time = merge_average(briefing->clock_execution_time, number_of_previous_runs, internal_briefing->clock_execution_time, 1);
	briefing->execution_time = merge_average(briefing->execution_time, number_of_previous_runs, internal_briefing->execution_time, 1);
	briefing->precompute_time = merge_average(briefing->precompute_time, number_of_previous_runs, internal_briefing->precompute_time, 1);
	briefing->output_modularity = merge_average(briefing->output_modularity, number_of_previous_runs, internal_briefing->output_modularity, 1);
}

int find_communities(dynamic_graph *dg, dynamic_weighted_graph *dwg, execution_settings *settings, dynamic_weighted_graph **community_graph, int **community_vector, algorithm_execution_briefing *briefing)
{
	int phase_counter;

	dynamic_weighted_graph *phase_output_community_graph;
	dynamic_weighted_graph *input_graph = dwg;

	phase_execution_briefing phase_briefing;

	double initial_phase_modularity, final_phase_modularity;

	FILE *output_communities_file;
	FILE *output_graphs_file;

	// For timing
	clock_t begin, end;
	double global_begin_wtime, global_end_wtime;
	double begin_wtime, end_wtime;
	double begin_pre_compute_wtime, end_pre_compute_wtime;
	double clock_time, wtime_omp;
	double global_wtime_omp;

	double minimum_phase_improvement = settings->minimum_phase_improvement;
	double minimum_iteration_improvement = settings->minimum_iteration_improvement;

	char *output_communities_filename = settings->output_communities_file;
	char *output_graphs_filename = settings->output_graphs_file;

	omp_set_num_threads(settings->number_of_threads);

	global_begin_wtime = omp_get_wtime();

	output_graphs_file = output_communities_file = NULL;

	if (!dwg || !valid_minimum_improvement(minimum_phase_improvement) || !valid_minimum_improvement(minimum_iteration_improvement))
	{
		printf("Invalid algorithm parameters!");
		briefing->execution_successful = 0;

		return 0;
	}

	*community_vector = NULL;
	phase_output_community_graph = NULL;

	phase_counter = 0;

	clock_time = 0;
	wtime_omp = 0;
	briefing->precompute_time = 0;

	final_phase_modularity = compute_modularity_init_weighted_reference_implementation_method(dwg);

	do
	{
		free(*community_vector);

		initial_phase_modularity = final_phase_modularity;

		begin = clock();
		begin_wtime = omp_get_wtime();

		if (!(settings->phase_executor_weighted(dwg, settings, &phase_output_community_graph, community_vector, &phase_briefing)))
		{
			printf("Bad phase #%d computation!\n", phase_counter);
			briefing->execution_successful = 0;

			return 0;
		}

		final_phase_modularity = phase_briefing.output_modularity;

		// Just for performance measurement
		end_wtime = omp_get_wtime();
		end = clock();
		clock_time += (double)(end - begin) / CLOCKS_PER_SEC;
		wtime_omp += end_wtime - begin_wtime;

		if (output_communities_file && !output_save_communities(output_communities_file, *community_vector, dwg->size))
		{
			printf("Couldn't save communities output of phase #%d!\n", phase_counter);
			briefing->execution_successful = 0;

			return 0;
		}

		if (output_graphs_file && !output_save_community_graph(output_graphs_file, phase_output_community_graph, phase_counter))
		{
			printf("Couldn't save graph output of phase #%d!\n", phase_counter);
			briefing->execution_successful = 0;

			return 0;
		}
		printf("\nEnd of Phase #%d\n\n", phase_counter);
		printf("\tInitial modularity:                 %f\n", initial_phase_modularity);
		printf("\tFinal modularity:                   %f\n", final_phase_modularity);
		printf("\tGain:                               %f\n", final_phase_modularity - initial_phase_modularity);
		printf("\tExecution time:                     %fs\n", end_wtime - begin_wtime);
		printf("\tExecution time over all threads:    %fs\n", (double)(end - begin) / CLOCKS_PER_SEC);
		printf("\tNumber of iterations:               %d\n", phase_briefing.number_of_iterations);
		// Clean memory
		// Avoids freeing initial input graph
		if (dwg != input_graph)
			dynamic_weighted_graph_free(dwg);

		// Prepare for next phase
		phase_counter++;
		dwg = phase_output_community_graph;
	} while (final_phase_modularity - initial_phase_modularity > minimum_phase_improvement);

	*community_graph = phase_output_community_graph;

	if (output_communities_file)
		fclose(output_communities_file);

	if (output_graphs_file)
		fclose(output_graphs_file);

	global_end_wtime = omp_get_wtime();

	global_wtime_omp = global_end_wtime - global_begin_wtime;

	briefing->execution_successful = 1;
	briefing->execution_time = wtime_omp;
	briefing->number_of_phases = phase_counter;
	briefing->clock_execution_time = clock_time;
	briefing->output_modularity = final_phase_modularity;
	briefing->global_execution_time = global_wtime_omp;

	return 1;
}

int select_phase_executors(execution_settings *settings)
{
	settings->phase_executor_weighted = parallel_phase_weighted;
	return 1;
}

int execute_community_detection(dynamic_graph *input_dg, dynamic_weighted_graph *input_dwg, execution_settings *settings, dynamic_weighted_graph **community_graph, int **community_vector, execution_briefing *briefing)
{
	algorithm_execution_briefing internal_briefing;

	clock_t global_begin, global_end;
	double global_clock_time;

	double global_begin_wtime, global_end_wtime;
	double global_wtime_omp;

	briefing->performed_runs = 0;

	// Select proper phase executors, depending on selected algorithm version
	if (!select_phase_executors(settings))
	{
		printf("Could not set phase executors!\n");
		return 0;
	}

	global_begin_wtime = omp_get_wtime();
	if (!find_communities(input_dg, input_dwg, settings, community_graph, community_vector, &internal_briefing))
	{
		printf("Could not complete execution!\n");
		return 0;
	}

	briefing->execution_time = internal_briefing.execution_time;
	briefing->clock_execution_time = internal_briefing.clock_execution_time;
	briefing->minimum_execution_time = internal_briefing.execution_time;
	briefing->minimum_clock_execution_time = internal_briefing.clock_execution_time;
	briefing->output_modularity = internal_briefing.output_modularity;

	briefing->performed_runs++;

	global_end_wtime = omp_get_wtime();

	global_wtime_omp = global_end_wtime - global_begin_wtime;

	briefing->global_execution_time = global_wtime_omp;

	briefing->execution_successful = 1;

	return 1;
}
