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
	cell_defaults.functions.custom_cell_rule = custom_function;
	cell_defaults.functions.contact_function = contact_function;

	find_cell_definition("apical")->is_movable = false;

	find_cell_definition("layer_6")->functions.update_phenotype = migrating_phenotype_function;
	find_cell_definition("layer_5")->functions.update_phenotype = migrating_phenotype_function;
	find_cell_definition("layer_4")->functions.update_phenotype = migrating_phenotype_function;
	find_cell_definition("layer_3")->functions.update_phenotype = migrating_phenotype_function;
	find_cell_definition("layer_2")->functions.update_phenotype = migrating_phenotype_function;

	find_cell_definition("layer_6")->functions.custom_cell_rule = custom_function;
	find_cell_definition("layer_5")->functions.custom_cell_rule = custom_function;
	find_cell_definition("layer_4")->functions.custom_cell_rule = custom_function;
	find_cell_definition("layer_3")->functions.custom_cell_rule = custom_function;
	find_cell_definition("layer_2")->functions.custom_cell_rule = custom_function;

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
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
		if (pCell->type_name != "rgc")
		{ continue; }

		pCell->phenotype.intracellular->set_parameter_value("notch", pCell->custom_data["notch"]);
		pCell->phenotype.intracellular->set_parameter_value("delta", pCell->custom_data["delta"]);
		pCell->phenotype.intracellular->set_parameter_value("nicd", pCell->custom_data["nicd"]);
		pCell->phenotype.intracellular->set_parameter_value("hesM", pCell->custom_data["hesM"]);
		pCell->phenotype.intracellular->set_parameter_value("hesC", pCell->custom_data["hesC"]);
		pCell->phenotype.intracellular->set_parameter_value("hesN", pCell->custom_data["hesN"]); 
		pCell->phenotype.intracellular->set_parameter_value("ngn", pCell->custom_data["ngn"]);
		pCell->phenotype.cycle.data.elapsed_time_in_phase = UniformRandom() / get_single_base_behavior(pCell, "cycle entry");
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
		if (neighbor->type_name != "rgc")
		{ continue; }

		neighbor_notch += neighbor->custom_data["notch"];
		neighbor_delta += neighbor->custom_data["delta"];
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

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
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

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 
