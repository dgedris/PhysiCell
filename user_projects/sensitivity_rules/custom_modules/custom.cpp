/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

#include "../addons/libRoadrunner/src/librr_intracellular.h"
#include <cmath>
#include <cfloat>
#include <array>
#include <cassert>

// Place random_int in an anonymous namespace to avoid multiple definition errors
namespace {
	int random_int(int min, int max) {
		static std::random_device rd;
		static std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis(min, max);
		return dis(gen);
	}
}

void create_cell_types( void )
{
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
	   Cell rule definitions 
	*/

	setup_behavior_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/

	cell_defaults.functions.update_phenotype = phenotype_function;
	cell_defaults.functions.custom_cell_rule = end_migration_pial_function;
	cell_defaults.functions.contact_function = contact_function;

	find_cell_definition("apical")->is_movable = false;

	find_cell_definition("layer_6")->functions.update_phenotype = migrating_phenotype_function;
	find_cell_definition("layer_5")->functions.update_phenotype = migrating_phenotype_function;
	find_cell_definition("layer_4")->functions.update_phenotype = migrating_phenotype_function;
	find_cell_definition("layer_3")->functions.update_phenotype = migrating_phenotype_function;
	find_cell_definition("layer_2")->functions.update_phenotype = migrating_phenotype_function;
	find_cell_definition("dying")->functions.update_phenotype = migrating_phenotype_function;
	find_cell_definition("G1_dead")->functions.update_phenotype = migrating_phenotype_function;

	find_cell_definition("layer_6")->functions.custom_cell_rule = end_migration_pial_function;
	find_cell_definition("layer_5")->functions.custom_cell_rule = end_migration_pial_function;
	find_cell_definition("layer_4")->functions.custom_cell_rule = end_migration_pial_function;
	find_cell_definition("layer_3")->functions.custom_cell_rule = end_migration_pial_function;
	find_cell_definition("layer_2")->functions.custom_cell_rule = end_migration_pial_function;
	find_cell_definition("dying")->functions.custom_cell_rule = end_migration_pial_function;
	find_cell_definition("G1_dead")->functions.custom_cell_rule = end_migration_pial_function;
	find_cell_definition("subplate")->functions.custom_cell_rule = basement_membrane_interaction_side;
	find_cell_definition("rgc")->functions.custom_cell_rule = basement_membrane_interaction_bottom;


	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
	// Set the spring anchoring of RGCs and subplate cells
	//find_cell_definition("rgc")->functions.custom_cell_rule = basement_membrane_interaction_bottom;

	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{

	// create some of each type of cell 
	
	cell_definitions_by_name["rgc"]->functions.pre_update_intracellular = rgc_pre_update_intracellular;
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();
	set_parameters_from_distributions();


	for (Cell* pCell : (*all_cells))
	{
		// Setup phase indexes
		int start_phase_index = 0;
		int end_phase_index = static_cast<int>(pCell->phenotype.cycle.model().phases.size()) - 1;

		if (pCell->type_name == "subplate")
		{
			double x_min = default_microenvironment_options.X_range[0];
			double x_max = default_microenvironment_options.X_range[1];
			double buffer = parameters.doubles("buffer_distance");
			double x_min_buf = x_min + buffer;
			double x_max_buf = x_max - buffer;

			if (pCell->position[0] >= x_min_buf && pCell->position[0] <= x_max_buf)
			{
				/*if ((pCell->position[0] >= 10.0 && pCell->position[0] <= 40.0) ||(pCell->position[0] >= 370.0 && pCell->position[0] <= 400.0))
				{
					pCell->functions.custom_cell_rule = basement_membrane_interaction_bottom;
					continue;
				}*/
				pCell->functions.custom_cell_rule = NULL;
			}
		}
		if (pCell->type_name != "rgc")
		{ continue; }

		pCell->phenotype.intracellular->set_parameter_value("notch", pCell->custom_data["notch"]);
		pCell->phenotype.intracellular->set_parameter_value("delta", pCell->custom_data["delta"]);
		pCell->phenotype.intracellular->set_parameter_value("nicd", pCell->custom_data["nicd"]);
		pCell->phenotype.intracellular->set_parameter_value("hesM", pCell->custom_data["hesM"]);
		pCell->phenotype.intracellular->set_parameter_value("hesC", pCell->custom_data["hesC"]);
		pCell->phenotype.intracellular->set_parameter_value("hesN", pCell->custom_data["hesN"]); 
		pCell->phenotype.intracellular->set_parameter_value("ngn", pCell->custom_data["ngn"]);
		// Randomly select the initial phase, only if valid
		if (end_phase_index >= start_phase_index) {
			pCell->phenotype.cycle.data.current_phase_index = random_int(start_phase_index, end_phase_index);
		} else {
			pCell->phenotype.cycle.data.current_phase_index = 0;
		}
		// Randomly set the elapsed time in phase
		pCell->phenotype.cycle.data.elapsed_time_in_phase = UniformRandom() / get_single_base_behavior(pCell, "exit from cycle phase "+ pCell->phenotype.cycle.data.current_phase_index);
		
	}
	
	return; 
}

void rgc_pre_update_intracellular( Cell* pCell, Phenotype& phenotype, double dt )
{
	// first use up the next of the delay terms
	RoadRunnerIntracellular* rri = getRoadRunnerModel(pCell);
	double neighbor_notch = 0.0;
	double neighbor_delta = 0.0;
	const std::vector<Cell *> *neighbors = &pCell->state.neighbors;
	for (auto neighbor : (*neighbors))
	{
		if (neighbor->type_name == "rgc")
		{neighbor_notch += neighbor->custom_data["notch"];
		neighbor_delta += neighbor->custom_data["delta"];}
		
		if (neighbor->type_name == "ipc")
		{neighbor_notch += neighbor->custom_data["notch"];
		neighbor_delta += neighbor->custom_data["delta"];}

		continue;
	}
	rri->set_parameter_value("neighbor_notch", neighbor_notch);
	rri->set_parameter_value("neighbor_delta", neighbor_delta);
	return;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void migrating_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	if( phenotype.motility.migration_speed > 0.2 )
	{
		phenotype.mechanics.cell_cell_repulsion_strength = 0;
	}
	else
	{
		phenotype.mechanics.cell_cell_repulsion_strength = find_cell_definition( pCell->type_name )->phenotype.mechanics.cell_cell_repulsion_strength; 
	}
	return;
}

void end_migration_pial_function( Cell* pCell, Phenotype& phenotype , double dt )
{
	// check if any neighbors are pial cells
	// if so, stop migration
	// if not, continue migration
	for (int i = 0; i < pCell->state.neighbors.size(); i++)
	{
		if (pCell->state.neighbors[i]->type_name == "pial")
		{
			phenotype.motility.migration_speed = 0;
			pCell->functions.custom_cell_rule = NULL;
			return;
		}
	}
	return;
}

void end_migration_subplate_function( Cell* pCell, Phenotype& phenotype , double dt )
{
	// check if any neighbors are pial cells
	// if so, stop migration
	// if not, continue migration
	for (int i = 0; i < pCell->state.neighbors.size(); i++)
	{
		if (pCell->state.neighbors[i]->type_name == "subplate")
		{
			phenotype.motility.migration_speed = 0;
			pCell->functions.custom_cell_rule = NULL;
			return;
		}
	}
	return;
}


// ################### Function that implements Cell-to-BM force ####################
void basement_membrane_interaction_side(Cell* pCell, Phenotype& phenotype, double dt)
{
	double cell_x = pCell->position[0];
    double x_min = default_microenvironment_options.X_range[0];
    double x_max = default_microenvironment_options.X_range[1];
	double px = 0.0;
	double del_x = 0.0;
	double x_min_buf = 0.0;
	double x_max_buf = 0.0;
	double buffer = parameters.doubles("buffer_distance");

	x_min_buf = x_min + buffer;
	x_max_buf = x_max - buffer;
	if (cell_x < x_min_buf){px = x_min;}
	else if (cell_x > x_max_buf){px = x_max;}
	else {
		std::cout << "Subplate Boundary Warning: Subplate cell (id=" << pCell->ID << ") is not in the buffer zone any longer. Currently at " << cell_x << "." <<std::endl;
		return;
	}	
	double R = pCell->phenotype.geometry.radius;

	del_x = px - cell_x;
    del_x = del_x - (del_x < 0 ? R : -R);

	double strength = parameters.doubles("membrane_adhesion_strength");

	pCell->velocity[0] += strength * del_x;
}

void basement_membrane_interaction_bottom(Cell* pCell, Phenotype& phenotype, double dt)
{
    double cell_y = pCell->position[1];
	double py = 8.41;
	double del_y = 0.0;
	double R = pCell->phenotype.geometry.radius;

	del_y = py - cell_y;
    del_y = del_y - (del_y < 0 ? R : -R);

	double strength = parameters.doubles("membrane_adhesion_strength_bottom");

	pCell->velocity[1] += strength * del_y;
}

/*void basement_membrane_interaction_bottom(Cell* pCell, Phenotype& phenotype, double dt)
{
	double cell_x = pCell->position[0];
    double cell_y = pCell->position[1];
	double py = 0.0;
	double px = cell_x;
	double dist = 0.0;
	
    dist = cell_y - py;

    //bool inside = is_inside(cell_x, cell_y, boundary_membrane_pts);
    //double d = inside ? -dist : dist;

    double R = pCell->phenotype.geometry.radius;
    double de = dist - (dist < 0 ? -R : R);

	// Implement replusion here
	double cell_deadzone = parameters.doubles("cell_deadzone");

	if (de>0){

		double displacement_needed = de + cell_deadzone; // Make cells only move to deadzone

		double nx = cell_x - px; // away from boundary
		double ny = cell_y - py;
		double norm = sqrt(nx * nx + ny * ny);
		if (norm > 1e-16)
		{
			nx /= norm;
			ny /= norm;
		}

		// Calculate the mag of the corrective velocity
		double correction_rate = parameters.doubles("membrane_correction_rate");
		double mag = correction_rate * displacement_needed;

		pCell->velocity[0] += mag * nx;
		pCell->velocity[1] += mag * ny;
	}

	else{ 
		double L = parameters.doubles("membrane_interaction_length");
		if (fabs(de) >= L) return;

		if (fabs(de) < cell_deadzone) return;

		double strength = parameters.doubles("membrane_adhesion_strength");
		double mag = strength * fabs(de);

		double nx = px - cell_x;
		double ny = py - cell_y;
		double norm = sqrt(nx * nx + ny * ny);
		if (norm > 0)
		{
			nx /= norm;
			ny /= norm;
		}

		double sign = (dist < 0.0 ? 1.0 : -1.0);

		pCell->velocity[0] += mag * nx;
		pCell->velocity[1] += mag * ny;
	}
}*/

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 
