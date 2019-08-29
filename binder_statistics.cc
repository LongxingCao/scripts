// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// headers


#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/motif/reference_frames.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/id/AtomID_Map.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

// residue selectors
#include <core/select/util.hh>
#include <core/select/residue_selector/ResidueNameSelector.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/ResiduePDBInfoHasLabelSelector.hh>

#include <devel/init.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <numeric/xyzVector.hh>
#include <protocols/idealize/IdealizeMover.hh>
// #include <protocols/sicdock/Assay.hh>
#include <utility/pointer/memory.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>

#include <boost/foreach.hpp>

#include <json.hpp>

// protocols
#include <core/pack/task/TaskFactory.hh>
#include <protocols/task_operations/ProteinInterfaceDesignOperation.hh>
#include <protocols/simple_moves/AlignChainMover.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.fwd.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <protocols/simple_filters/InterfaceSasaFilter.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/TaskAwareMinMover.hh>
#include <protocols/simple_ddg/DdgFilter.hh>

#include <protocols/simple_filters/ShapeComplementarityFilter.cc>

#include <fstream>
#include <algorithm>
#include <map>

#include <core/kinematics/FoldTree.hh>
#include <core/pose/subpose_manipulation_util.hh>

static basic::Tracer TR( "pilot.longxing.residue_statistics" );


//OPT_1GRP_KEY( String, longxing, ref_pdb )
//OPT_1GRP_KEY( Real, longxing, ddg_threshold )
//OPT_1GRP_KEY( Boolean, longxing, multi_segs)
OPT_1GRP_KEY( Boolean, longxing, pretty )
OPT_1GRP_KEY( Boolean, longxing, compute_ddg )
OPT_1GRP_KEY( Boolean, longxing, compute_sasa )
OPT_1GRP_KEY( Boolean, longxing, compute_sc )
OPT_1GRP_KEY( String, longxing, binder_chain )
OPT_1GRP_KEY( String, longxing, prune_residues )
OPT_1GRP_KEY( Real, longxing, residue_sasa_probe_radius )
OPT_1GRP_KEY( Boolean, longxing, residue_sasa_include_H )
OPT_1GRP_KEY( Real, longxing, neighborhood_distance )
OPT_1GRP_KEY( Boolean, longxing, compute_sidechain_neighbors )
OPT_1GRP_KEY( Real, longxing, scnb_angle_exponent )
OPT_1GRP_KEY( Real, longxing, scnb_angle_shift_factor )
OPT_1GRP_KEY( Real, longxing, scnb_dist_exponent )
OPT_1GRP_KEY( Real, longxing, scnb_dist_midpoint )
OPT_1GRP_KEY( Real, longxing, scnb_neighbor_denominator )
OPT_1GRP_KEY( String, longxing, scaffold_info_fname )
OPT_1GRP_KEY(  StringVector, longxing, pdbinfo_labels )
void register_options() {
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    //NEW_OPT( longxing::ref_pdb                 , "the reference pdb file"                                      , ""    );
    //NEW_OPT( longxing::ddg_threshold           , "the ddg cutoff value of the motif"                           , -15.0 );
    //NEW_OPT( longxing::multi_segs              , "dump multi segments or not if they are connected by a loop"  , false );
    NEW_OPT( longxing::pretty                       , "bueatify the json format"                                            , false );
    NEW_OPT( longxing::compute_ddg                  , "compute the real ddg, which is very time consuming"                  , false );
    NEW_OPT( longxing::compute_sasa                 , "compute the interface buried surface area"                           , false );
    NEW_OPT( longxing::compute_sc                   , "compute the interface shape complementarity"                         , false );
    NEW_OPT( longxing::binder_chain                 , "chain identifier of the binder"                                      , "A"   );
    NEW_OPT( longxing::prune_residues               , "purne these residues to ALA before analysis"                         , ""    );
    NEW_OPT( longxing::residue_sasa_probe_radius    , "probe radius used to calculate the residue level sasa"               , 2.0   );
    NEW_OPT( longxing::residue_sasa_include_H       , "include hydrogen atom for residue sasa calculateion  "               , false );
    NEW_OPT( longxing::neighborhood_distance        , "neighborhood distance for interface residues selection"              , 8.0   );
    NEW_OPT( longxing::compute_sidechain_neighbors  , "do the sidechain neighbor calculation or not"                        , false );
    NEW_OPT( longxing::scnb_angle_exponent          , "see layselector"                                                     , 2.0   );
    NEW_OPT( longxing::scnb_angle_shift_factor      , "see layselector"                                                     , 0.5   );
    NEW_OPT( longxing::scnb_dist_exponent           , "see layselector"                                                     , 1.0   );
    NEW_OPT( longxing::scnb_dist_midpoint           , "see layselector"                                                     , 9.0   );
    NEW_OPT( longxing::scnb_neighbor_denominator    , "see layselector"                                                     , 1.0   );
    NEW_OPT( longxing::scaffold_info_fname          , "chain identifier of the binder"                                      , ""   );
		NEW_OPT( longxing::pdbinfo_labels               , "pdb info labels that you want to get the residue index              ", utility::vector1<std::string>());
}

struct PilotOptions
{
		bool pretty;
		bool compute_ddg;
		bool compute_sasa;
		bool compute_sc;
    utility::vector1< std::string > input_pdbs;
    char binder_chain;
    std::string prune_residues;
    core::Real residue_sasa_probe_radius;
    bool residue_sasa_include_H;
    core::Real neighborhood_distance;
    
    // for layer definition
    bool       compute_sidechain_neighbors;
    core::Real scnb_angle_exponent;
    core::Real scnb_angle_shift_factor;
    core::Real scnb_dist_exponent;
    core::Real scnb_dist_midpoint;
    core::Real scnb_neighbor_denominator;

		std::string scaffold_info_fname;
		std::map<std::string, std::string> scaffold_info;

		utility::vector1<std::string> pdbinfo_labels;
    
    void init_from_cli()
    {
        using namespace basic::options;
        using namespace basic::options::OptionKeys;
        
				pretty        = option[ longxing::pretty       ]();
				compute_ddg   = option[ longxing::compute_ddg  ]();
				compute_sasa  = option[ longxing::compute_sasa ]();
				compute_sc    = option[ longxing::compute_sc   ]();

        // parse input pdbs
        if ( option[ in::file::s ].user() ) {
            input_pdbs = option[ in::file::s ]();
        } else if ( option[ in::file::l ].user() ) {
            
            utility::vector1<utility::file::FileName> list = basic::options::option[ in::file::l ]();
            std::string fname;
            for ( unsigned int h=1; h<=list.size(); h++ ) {
                utility::io::izstream pdbs(list[h]);
                if ( pdbs ) {
                    while ( pdbs >> fname ) {
                        if ( fname.empty() || fname.find('#') == 0 ) continue;
                        input_pdbs.push_back(fname);
                    }
                }
            }
        } else {
            TR.Error << "You must specific input pdbs using the name of the pdb or just a list file" << std::endl;
            utility_exit_with_message("Exit program because no input structures were specificed!!!");
        }
        
        // define the binder chain
        binder_chain = option[ longxing::binder_chain ]()[0];
        prune_residues = option[ longxing::prune_residues ]();
        residue_sasa_probe_radius = option[ longxing::residue_sasa_probe_radius ]();
        residue_sasa_include_H     = option[ longxing::residue_sasa_include_H ]();
        neighborhood_distance      = option[ longxing::neighborhood_distance  ]();
        
        //sidechain neighbor calculation
        compute_sidechain_neighbors = option[ longxing::compute_sidechain_neighbors ]();
        scnb_angle_exponent         = option[ longxing::scnb_angle_exponent ]();
        scnb_angle_shift_factor     = option[ longxing::scnb_angle_shift_factor ]();
        scnb_dist_exponent          = option[ longxing::scnb_dist_exponent ]();
        scnb_dist_midpoint          = option[ longxing::scnb_dist_midpoint ]();
        scnb_neighbor_denominator   = option[ longxing::scnb_neighbor_denominator ]();

				// initialize the scaffold info map
				if ( option[ longxing::scaffold_info_fname ].user() ) {
						scaffold_info_fname = option[ longxing::scaffold_info_fname ]();
						if ( false == utility::file::file_exists( scaffold_info_fname ) ) {
								utility_exit_with_message("Exit program because the scaffold info file doesn't exist!!!!!!!!!!!!");
						}
						std::string line;
						std::ifstream f( scaffold_info_fname );
						while( std::getline( f, line ) ) {
								if ( line.empty() || line.find('#') == 0 ) continue;
								utility::vector1<std::string> splited_line = utility::quoted_split( line );
								runtime_assert_msg( splited_line.size() == 3, "Currently the format of the scaffold info file should be \"# scaffold_tag scaffold_path scaffold_sequence\"");
								scaffold_info[ splited_line[1] ] = line;
						}
				}

				pdbinfo_labels = option[ longxing::pdbinfo_labels ]();
    }
    
};

void pose_to_polyXXX( core::pose::Pose & pose, std::string name3, utility::vector1<core::Size> const & res_sel, bool ignore_gpc = false ){
	core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
	core::conformation::ResidueOP ala = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map(name3) );
	for( int iri = 1; iri <= res_sel.size(); ++iri ){
		int ir = res_sel[iri];
		if( ! pose.residue(ir).is_protein()   ) continue;
		if( ignore_gpc && pose.residue(ir).name3()=="GLY" ) continue;
		if( ignore_gpc && pose.residue(ir).name3()=="PRO" ) continue;
		if( ignore_gpc && pose.residue(ir).name3()=="CYD" ) continue;
		pose.replace_residue( ir, *ala, true );
	}
}

utility::vector1<core::Real> residue_sasa( core::pose::Pose const & pose, PilotOptions const & opt )
{
    using core::Real;
    
    Real total_sasa;
    core::id::AtomID_Map<Real> atom_sasa;
    utility::vector1<Real> rsd_sasa;
    Real const probe_radius = opt.residue_sasa_probe_radius;
    bool const sasa_include_H = opt.residue_sasa_include_H;
    core::id::AtomID_Map<bool> atom_subset;
    atom_subset.resize( pose.size() );
    // create atom mask, check include hydrogens or not??
    for ( unsigned int ires = 1; ires <= pose.size(); ++ires )
    {
        atom_subset.resize( ires, pose.residue_type(ires).natoms(), false);
        for ( unsigned int iatom = 1; iatom<= pose.residue_type(ires).nheavyatoms(); ++iatom)
        {
            atom_subset[ires][iatom] = true;
        }
        if ( sasa_include_H )
        {
            for ( unsigned int iatom = pose.residue_type(ires).nheavyatoms()+1; iatom <= pose.residue_type(ires).natoms(); ++iatom )
            {
                atom_subset[ires][iatom] = true;
            }
        }
    }
    // let's use the reduce radii here, as it seems to be the default radii.'
    // no big H
    total_sasa = core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, false, atom_subset, false);
    
    TR << "Total sasa of the complex computed with calc_per_atom_sasa method: " << total_sasa << std::endl;
    
    return rsd_sasa;
}

utility::vector1<core::Size> get_interface_residues( core::pose::Pose const & pose, PilotOptions const & opt )
{
    using namespace core::select::residue_selector;
    runtime_assert_msg( pose.num_chains() == 2, "I could only deal with 2 chains" );
    NotResidueSelectorOP target( utility::pointer::make_shared<NotResidueSelector>( utility::pointer::make_shared<ChainSelector>( opt.binder_chain ) ) );
    NeighborhoodResidueSelectorOP interface( utility::pointer::make_shared<NeighborhoodResidueSelector>( target, opt.neighborhood_distance, false ) );
    
    utility::vector1<core::Size> interface_residues =  core::select::get_residues_from_subset( interface->apply(pose) );
    core::Size index_shift = pose.pdb_info()->chain(1) == opt.binder_chain ? 0: pose.chain_end(1);
    for ( core::Size ii = 1; ii <= interface_residues.size(); ++ii )
    {
        interface_residues[ii] -= index_shift;
    }
    
    return interface_residues;
}

// copied from burial_utilities.hh
inline core::Real calculate_point_in_cone(
                                          numeric::xyzVector<core::Real> const &point_coordinates,
                                          numeric::xyzVector<core::Real> const &conevect,
                                          numeric::xyzVector<core::Real> const &conevect_coordinate_2,
                                          core::Real const angle_exponent,
                                          core::Real const angle_shift_factor,
                                          core::Real const dist_exponent,
                                          core::Real const dist_midpoint
                                          ) {
		using core::Real;
    numeric::xyzVector< Real > vect(point_coordinates - conevect_coordinate_2);
    core::Real const dist_term(1.0 / (1.0 + exp( dist_exponent*(vect.length() - dist_midpoint)  ))); // dist_term = 1/(1+exp( n*(d - m))); sigmoidal falloff with midpoint at m, with sharpness controlled by n
    core::Real angle_term( ( conevect.dot(vect.normalize()) + angle_shift_factor ) / (1 + angle_shift_factor ) );
    if ( angle_term < 0 ) {
        angle_term = 0.0;
    }
    return (dist_term * pow(angle_term, angle_exponent) );
}

// copied directly from SelectResidueByLayer.cc. The function is a private function in that file and I cann't use it directly.
utility::vector1< core::Real > const
calc_sc_neighbors( core::pose::Pose const & pose, PilotOptions const & opt ) {
    
    using core::Real;
    using core::Size;
    
    utility::vector1< Real > rsd_sc_neighbors;
    
		core::Size binder_resi_start(0), binder_resi_end(0);
    if ( pose.pdb_info()->chain(1) == opt.binder_chain ) {
        binder_resi_start = 1;
        binder_resi_end   = pose.chain_end(1);
    } else {
				binder_resi_start = pose.chain_end(1)+1;
        binder_resi_end   = pose.size();
    }
    for ( Size i = binder_resi_start; i <= binder_resi_end; ++i ) {
        Real my_neighbors(0.0);
        
        numeric::xyzVector< Real > my_sc_coordinates;
        numeric::xyzVector< Real > my_bb_coordinates;
        
        if ( pose.residue( i ).name3() == "GLY" ) {
            my_sc_coordinates = pose.residue(i).atom(pose.residue(i).atom_index("2HA")).xyz() ;
            my_bb_coordinates = pose.residue(i).atom(pose.residue(i).atom_index("CA")).xyz() ;
        } else {
            if ( pose.residue(i).is_polymer() ) {
                my_sc_coordinates = pose.residue(i).atom(pose.residue(i).first_sidechain_atom()).xyz() ;
                core::Size parent_atom_index = pose.residue(i).icoor( pose.residue(i).first_sidechain_atom() ).stub_atom1().atomno();
                my_bb_coordinates = pose.residue(i).atom( parent_atom_index ).xyz() ;
            } else {
                rsd_sc_neighbors.push_back(0); //For now, ligands do not have their neighbours counted.  This could change in the future.
                continue;
            }
        }
        
        numeric::xyzVector< Real > my_sc_vector = (my_sc_coordinates - my_bb_coordinates).normalize() ;
        
        for ( Size j = 1; j <= pose.size(); ++j ) {
            
            if ( i != j ) {
                
                numeric::xyzVector< Real > other_bb_coordinates;
                if ( pose.residue(j).name3() == "GLY" ) {
                    other_bb_coordinates = pose.residue(j).atom(pose.residue(j).atom_index("CA")).xyz();
                } else {
                    if ( pose.residue(j).is_polymer() ) { //If this is a polymer atom, use the parent of the first sidechain atom.
                        core::Size parent_atom_index = pose.residue(j).icoor( pose.residue(j).first_sidechain_atom() ).stub_atom1().atomno();
                        other_bb_coordinates = pose.residue(j).atom( parent_atom_index ).xyz();
                    } else { //If this is not a polymer residue, use the nbr_atom:
                        core::Size nbr_atom_index = pose.residue(j).nbr_atom();
                        other_bb_coordinates = pose.residue(j).atom( nbr_atom_index ).xyz();
                    }
                }
                
                my_neighbors += calculate_point_in_cone( other_bb_coordinates, my_sc_vector, my_sc_coordinates, opt.scnb_angle_exponent, opt.scnb_angle_shift_factor, opt.scnb_dist_exponent, opt.scnb_dist_midpoint );
            }
        }
        rsd_sc_neighbors.push_back(my_neighbors / opt.scnb_neighbor_denominator);
        //TR << i << "  " << my_neighbors << std::endl;
    }
    return rsd_sc_neighbors;
}

void get_scaffold_information( core::pose::Pose const & pose, PilotOptions const & opt, std::string & scaffold_path, std::string & scaffold_sequence )
{
		std::string description = pose.pdb_info()->name();

		// someday I should just rename all the scaffolds, it's strange to have such long names.'
		utility::vector1<std::string> possible_scaffolds;
		//for ( std::map<std::string, std::string>::iterator it = opt.scaffold_info.begin(); it != opt.scaffold_info.end(); ++it )
		for ( auto it = opt.scaffold_info.begin(); it != opt.scaffold_info.end(); ++it )
		{
				if ( description.find( it->first ) != std::string::npos ) {
						possible_scaffolds.push_back( it->second );
				}
		}

		if ( possible_scaffolds.size() == 0 ) {
        utility_exit_with_message("Can not find the original scaffold!!!!");
		} else if ( possible_scaffolds.size() == 1 ) {
				utility::vector1<std::string> splited_line = utility::quoted_split( possible_scaffolds[1] );
				scaffold_path     = splited_line[2];
				scaffold_sequence = splited_line[3];
				return;
		} else {
				core::Size mutations = 9999;
				std::string design_sequence = pose.sequence();
				for ( std::string & line : possible_scaffolds ) {
						utility::vector1<std::string> splited_line = utility::quoted_split( line );
						// same length, grafting ?? same length ??
						if ( splited_line[3].size() != design_sequence.size() ) continue;
						core::Size curr_mutations(0);
						for ( core::Size ires = 0; ires < design_sequence.size(); ++ires ) {
								if ( design_sequence[ires] != splited_line[3][ires] ) ++curr_mutations;
						}
						if ( curr_mutations < mutations ) {
								mutations = curr_mutations;
								scaffold_path = splited_line[2];
								scaffold_sequence = splited_line[3];
						}
				}

				if ( mutations > design_sequence.size() ) {
						utility_exit_with_message( "Even though I find some possible scaffold, but non of them work!!");
				}
				return;
		}
}


void compute_binder_info( std::string fname, core::scoring::ScoreFunctionOP sfxn, PilotOptions const & opt)
{
    using core::Real;
    using core::Size;
    
    std::string ftag;
    utility::vector1<std::string> temp_split = utility::string_split(fname, '/');
    std::string basename = temp_split.back();
    if   ( fname.substr(fname.size()-4,4)==".pdb" )         ftag = basename.substr(0,basename.size()-4);
    else if ( fname.substr(fname.size()-7,7)==".pdb.gz" )   ftag = basename.substr(0,basename.size()-7);
    else                                                    utility_exit_with_message("Unknown PDB format!!!");
    
    nlohmann::json json_results;
    
    core::pose::PoseOP pPose     = core::import_pose::pose_from_file(fname);
    
    if ( 0 != opt.prune_residues.size() ) {
        core::select::residue_selector::ResidueNameSelectorCOP resname_selector( utility::pointer::make_shared< core::select::residue_selector::ResidueNameSelector >(opt.prune_residues, false /* dummy */ ) );
				core::select::residue_selector::ResidueSubset resi_to_prune = resname_selector->apply( *pPose );
				TR << "prune this residues to ALA before analysis: ";
				for( int resi : core::select::get_residues_from_subset(resi_to_prune) ) {
						TR << resi << "  ";
				}
				TR << std::endl;
        protocols::simple_moves::MutateResidueOP mut( utility::pointer::make_shared<protocols::simple_moves::MutateResidue>() );
        mut->set_selector(resname_selector);
        mut->set_res_name( "ALA" );
        mut->apply( *pPose );
    }
    
    core::pose::PoseOP pbinder, ptarget;
    runtime_assert_msg( pPose->num_chains() == 2, "I could only deal with 2 chains" );
    utility::vector1<core::pose::PoseOP> splited_complex = pPose->split_by_chain();
    unsigned int binder_chain_index(0);
    unsigned int binder_length(0);
    core::pose::PoseOP pose_binder_polyGLY = pPose->clone();
    utilitY::vector1<int> residues_to_mutate;
    if ( pPose->pdb_info()->chain(1) == opt.binder_chain ) {
        binder_chain_index = 1;
        binder_length = pPose->chain_end(1);
        for( int ii = 1; ii <= pPose->chain_end(1); ++ii )residues_to_mutate.push_back(ii);
    } else {
        binder_chain_index = 2;
        binder_length = pPose->size() - pPose->chain_end(1);
        for( int ii = pPose->chain_end(1) + 1; ii <= pPose->size(); ++ii )residues_to_mutate.push_back(ii);
    }
    pose_to_polyXXX(pose_binder_polyGLY, "GLY", residues_to_mutate);
    pbinder = splited_complex[binder_chain_index];
		// make sure there if pdb_info object
		pbinder->pdb_info()->name(ftag);
    TR << "The " << binder_chain_index << " chain is selected as the binder." << std::endl;
    
    Real score = sfxn->score(*pPose);
    TR << "Total score of the complex: " << score << std::endl;
    core::scoring::Energies    const & energies       ( pPose->energies() );
    core::scoring::EnergyGraph const & energy_graph   ( energies.energy_graph() );
    core::scoring::EMapVector  const & energy_weights ( energies.weights() );
    utility::vector1<Real> residue_ddg;
    Real ddg_norepack(0.0);
    
    for ( Size ires = 1; ires <=pPose->size(); ++ires ){
        if ( pPose->chain(ires) != binder_chain_index ) continue;
        Real sum = 0;
        for ( utility::graph::Graph::EdgeListConstIter
             iru  = energy_graph.get_node(ires)->const_edge_list_begin(),
             irue = energy_graph.get_node(ires)->const_edge_list_end();
             iru != irue; ++iru
             ) {
            core::scoring::EnergyEdge const & edge( static_cast< core::scoring::EnergyEdge const & > (**iru) );
            Size const jres( edge.get_other_ind(ires) );
            if ( pPose->chain(jres) == binder_chain_index ) continue;
            
            sum += edge.dot( energy_weights );
        }
        residue_ddg.push_back( sum );
        ddg_norepack += sum;
    }
    TR << "The ddg_norepack score of the binder is: " << ddg_norepack << std::endl;
    // TODO: lots of duplicate code, clean the code
    sfxn->score(*pose_binder_polyGLY);
    core::scoring::Energies    const & pose_binder_polyGLY_energies       ( pose_binder_polyGLY->energies() );
    core::scoring::EnergyGraph const & pose_binder_polyGLY_energy_graph   ( pose_binder_polyGLY_energies.energy_graph() );
    core::scoring::EMapVector  const & pose_binder_polyGLY_energy_weights ( pose_binder_polyGLY_energies.weights() );
    utility::vector1<Real> residue_backbone_ddg;
    Real ddg_backbone_norepack(0.0);
    
    for ( Size ires = 1; ires <=pPose->size(); ++ires ){
        if ( pPose->chain(ires) != binder_chain_index ) continue;
        Real sum = 0;
        for ( utility::graph::Graph::EdgeListConstIter
             iru  = pose_binder_polyGLY_energy_graph.get_node(ires)->const_edge_list_begin(),
             irue = pose_binder_polyGLY_energy_graph.get_node(ires)->const_edge_list_end();
             iru != irue; ++iru
             ) {
            core::scoring::EnergyEdge const & edge( static_cast< core::scoring::EnergyEdge const & > (**iru) );
            Size const jres( edge.get_other_ind(ires) );
            if ( pPose->chain(jres) == binder_chain_index ) continue;
            
            sum += edge.dot( pose_binder_polyGLY_energy_weights );
        }
        residue_backbone_ddg.push_back( sum );
        ddg_backbone_norepack += sum;
    }

		// real ddg
		if ( opt.compute_ddg ) {
				protocols::task_operations::ProteinInterfaceDesignOperationOP pack_long( new protocols::task_operations::ProteinInterfaceDesignOperation() );
				pack_long->repack_chain1( true );
				pack_long->repack_chain2( true );
				pack_long->design_chain1( false );
				pack_long->design_chain2( false );
				pack_long->interface_distance_cutoff( 15.0 );
				core::pack::task::TaskFactoryOP task_factory( new core::pack::task::TaskFactory() );
				task_factory->push_back( pack_long );
				protocols::minimization_packing::MinMoverOP minmover( new protocols::minimization_packing::MinMover() );
				protocols::minimization_packing::TaskAwareMinMoverOP taskaware_minmover( new protocols::minimization_packing::TaskAwareMinMover(minmover, task_factory) );
				taskaware_minmover->bb( false );
				taskaware_minmover->chi( true );
				taskaware_minmover->jump( false );
				protocols::simple_ddg::DdgFilterOP ddg_filter( new protocols::simple_ddg::DdgFilter( 999 /*Threshold*/, sfxn /*scorefunction*/, 1 /*jump*/, 5/*repeats*/) );
				ddg_filter->relax_mover( taskaware_minmover );
				core::Real ddg = ddg_filter->compute( *pPose );
				json_results["ddg"] = ddg;
				TR << "The ddg score is: " << ddg << std::endl;
		}
    
    utility::vector1<core::Size> interface_residues = get_interface_residues( *pPose, opt );
    
    // interface sasa and interface sc
		if ( opt.compute_sasa ) {
				protocols::simple_filters::InterfaceSasaFilterOP interface_sasa( new protocols::simple_filters::InterfaceSasaFilter() );
				Real interface_sasa_value = interface_sasa->compute(*pPose);
				json_results["interface_sasa"] = interface_sasa_value;
				TR << "interface_sasa: " << interface_sasa_value << std::endl;

		}
		if ( opt.compute_sc ) {
				protocols::simple_filters::ShapeComplementarityFilterOP interface_sc( new protocols::simple_filters::ShapeComplementarityFilter() );
				Real interface_sc_value   = interface_sc->compute(*pPose).sc;
				Real interface_median_dist_value = interface_sc->compute(*pPose).distance;
				json_results["interface_median_dist"] = interface_median_dist_value;
				json_results["interface_sc"] = interface_sc_value;
				TR << "interface_sc: " << interface_sc_value << std::endl;
				TR << "median dist: " << interface_median_dist_value << std::endl;
		}
    
    // residue sasa info
    utility::vector1<Real> rsd_sasa_all = residue_sasa( *pPose, opt );
    utility::vector1<Real> rsd_sasa_monomer = residue_sasa( *pbinder, opt );
		utility::vector1<Real> rsd_sasa_complex;
    for ( unsigned int ires = 1; ires <= pPose->size(); ++ires )
    {
        if ( pPose->chain(ires) == binder_chain_index ) rsd_sasa_complex.push_back( rsd_sasa_all[ires] );
    }
    
    Real score_monomer = sfxn->score( *pbinder );
    TR << "The monomer score of the binder is: " << score_monomer << std::endl;
    TR << "Score_per_res of the monomer is: " << score_monomer / binder_length << std::endl;
    
    // sequence info of the binder
    std::string binder_seq = pbinder->sequence();
    TR << "The sequence of the binder is: " << binder_seq << std::endl;
    
    // dssp
    core::scoring::dssp::Dssp dssp(*pbinder);
    std::string dssp_str = dssp.get_dssp_secstruct();
    // fix the N-ter and C-ter problems
    assert( dssp_str.length() >=2 );
    if( dssp_str[0] == 'L' && dssp_str[0] != dssp_str[1] ) dssp_str[0] = dssp_str[1];
    if( dssp_str[ dssp_str.length() -1 ] == 'L' && dssp_str[ dssp_str.length() -1 ] != dssp_str[ dssp_str.length() -2 ] ) dssp_str[ dssp_str.length() -1 ] = dssp_str[ dssp_str.length() -2 ];
    TR << "Secondary structure of the binder: " << dssp_str << std::endl;
    
    json_results["description"] = ftag;
    json_results["ddg_norepack"] = ddg_norepack;
    json_results["ddg_backbone_norepack"] = ddg_backbone_norepack;
    json_results["score_monomer"] = score_monomer;
    json_results["length"] = binder_length;
    json_results["score_per_res"] = score_monomer / binder_length;
    json_results["sequence"] = binder_seq;
    json_results["scondary_structure"] = dssp_str;
    // has to do this, wired bug!
		json_results["ddg_per_res"] = residue_ddg;
    json_results["ddg_backbone_per_res"] = residue_backbone_ddg;
    json_results["residue_sasa_complex"] = rsd_sasa_complex;
		json_results["residue_sasa_monomer"] = rsd_sasa_monomer;
    json_results["interface_residues"] = interface_residues;
    json_results["pdb_index_shift"] = pPose->pdb_info()->chain(1) == opt.binder_chain ? 0: pPose->chain_end(1);

		if ( opt.compute_sidechain_neighbors ) {
				utility::vector1<core::Real> num_sidechain_neighbors = calc_sc_neighbors( *pPose, opt );
				json_results["num_sidechain_neighbors"] = num_sidechain_neighbors;
		}

		if ( opt.scaffold_info_fname.size() != 0 ) {
				std::string scaffold_path;
				std::string scaffold_sequence;
				get_scaffold_information( *pbinder, opt, scaffold_path, scaffold_sequence);
				json_results["scaffold_path"] = scaffold_path;
				json_results["scaffold_sequence"] = scaffold_sequence;
		}

		if ( opt.pdbinfo_labels.size() != 0 ) {
				core::select::residue_selector::ResiduePDBInfoHasLabelSelector pdbinfo_selector;
				for ( std::string const & label : opt.pdbinfo_labels ) {
						pdbinfo_selector.set_label( label );
						utility::vector1<core::Size> residues_on_binder_with_label;
						utility::vector1<core::Size> selected_residues =  core::select::get_residues_from_subset( pdbinfo_selector.apply(*pPose) );
						core::Size index_shift = pPose->pdb_info()->chain(1) == opt.binder_chain ? 0: pPose->chain_end(1);
						for ( core::Size ii = 1; ii <= selected_residues.size(); ++ii )
						{
								core::Size shifted_index = selected_residues[ii] - index_shift;
								if ( shifted_index > 0 ) {
										residues_on_binder_with_label.push_back( shifted_index );
								}
						}
						json_results[label] = residues_on_binder_with_label;
				}
		}
    
    std::ofstream myfile( ftag + ".json");
		if ( opt.pretty ) myfile << std::setw(4);
    myfile << json_results << std::endl;
    myfile.close();
    
}

int main(int argc, char *argv[]) {
    try{
        
        using core::Real;
        using core::Size;
        
        register_options();
        devel::init(argc,argv);
        PilotOptions opt;
        opt.init_from_cli();
        
        
        // register the sasa calculator
        core::pose::metrics::CalculatorFactory & calculator_factory = core::pose::metrics::CalculatorFactory::Instance();
        if ( ! calculator_factory.check_calculator_exists("sasa") ) {
            core::pose::metrics::PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy() );
            calculator_factory.register_calculator("sasa", sasa_calculator);
        }
        
        // scoring function
        core::scoring::ScoreFunctionOP sf = core::scoring::ScoreFunctionFactory::create_score_function("beta_nov16");
        core::scoring::methods::EnergyMethodOptions myopt = sf->energy_method_options();
        myopt.hbond_options().decompose_bb_hb_into_pair_energies(true);
        sf->set_energy_method_options(myopt);
        
        for( std::string pdb : opt.input_pdbs ) {
            compute_binder_info( pdb, sf, opt );
        }
        
        TR << "Job done!" << std::endl;
        
    } catch ( utility::excn::Exception const & e ) {
        TR.Error << "caught exception " << e.msg() << std::endl;
        std::exit( 1 );
    }
    return 0;
}
