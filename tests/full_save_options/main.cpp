/*
 Test driver for the full save options set as attributes of
 <save><full_data><enable ...>. It sets up a microenvironment and a few cells
 (with regular and spring attachments), then writes two full saves so that
 both the "create the XML DOM" and "update the XML DOM" paths are exercised.

 Usage: ./full_save_options <config.xml>
*/

#include <cstdio>
#include <iostream>

#include "../../core/PhysiCell.h"
#include "../../modules/PhysiCell_standard_modules.h"

using namespace BioFVM;
using namespace PhysiCell;

int main( int argc, char* argv[] )
{
	if( argc < 2 )
	{
		std::cout << "Usage: " << argv[0] << " <config.xml>" << std::endl;
		return 1;
	}
	if( !load_PhysiCell_config_file( argv[1] ) )
	{ return 1; }

	omp_set_num_threads( 1 );

	initialize_microenvironment();
	create_cell_container_for_microenvironment( microenvironment, 30 );

	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	initialize_cell_definitions_from_pugixml();
	build_cell_definitions_maps();

	Cell* pCells[4];
	for( int i=0 ; i < 4 ; i++ )
	{
		pCells[i] = create_cell( *cell_definitions_by_index[0] );
		pCells[i]->assign_position( -15.0 + 10.0*i , 0.0 , 0.0 );
	}
	attach_cells( pCells[0] , pCells[1] );
	attach_cells_as_spring( pCells[2] , pCells[3] );

	set_save_biofvm_mesh_as_matlab( true );
	set_save_biofvm_data_as_matlab( true );
	set_save_biofvm_cell_data( true );
	set_save_biofvm_cell_data_as_custom_matlab( true );

	char filename[1024];
	for( int n=0 ; n < 2 ; n++ )
	{
		sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str() , n );
		save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time );
		PhysiCell_globals.current_time += PhysiCell_settings.full_save_interval;
	}

	return 0;
}
