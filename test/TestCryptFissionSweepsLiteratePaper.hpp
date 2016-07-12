/*
 * TestCryptFissionLiteratePaper.hpp
 *
 * Created on: Jul 08 2016
 * Last modified: Jul 12 2016
 * 		Author: Axel A. Almet
 */

#ifndef TESTCRYPTFISSIONSWEEPSLITERATEPAPER_HPP_
#define TESTCRYPTFISSIONSWEEPSLITERATEPAPER_HPP_
/*
 * = Analysis of fission in epithelial layer =
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this test we show how Chaste can be used to analyse the incidence of buckling in a layer of epithelial cells.
 * Details of the computational model can be found in
 * Almet et al (2016) "Paneth cell-rich regions separated by a cluster of Lgr5+ cells initiate
 * fission in the intestinal stem cell niche".
 *
 * == Including header files ==
 *
 * EMPTYLINE
 *
 * We begin by including the necessary header files. The first ones are common to all cell_based Chaste simulations
 */

#include <cxxtest/TestSuite.h> //Needed for all test files
#include "CellBasedSimulationArchiver.hpp" //Needed if we would like to save/load simulations
#include "AbstractCellBasedTestSuite.hpp" //Needed for cell-based tests: times simulations, generates random numbers and has cell properties
#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "SmartPointers.hpp" //Enables macros to save typing

/* The next set of classes are needed specifically for the simulation, which can be found in the core code. */

#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "DifferentiatedCellProliferativeType.hpp" //Non-proliferative cell type
#include "TransitCellProliferativeType.hpp" //Proliferative cell type
#include "WildtypeCellMutationState.hpp" //Standard wildtype state for all cells
#include "FakePetscSetup.hpp" //Forbids tests running in parallel

/* This set of classes are not part of the core code and were made for this model. */

#include "StochasticTargetProportionBasedCellCycleModel.hpp" //Asymmetric-division-based cell cycle model
#include "HardCellMutationState.hpp" //Mutation class that defines hard cells

#include "EpithelialLayerLinearSpringForce.hpp" //Spring force law to account for different cell type pairs
#include "EpithelialLayerBasementMembraneForce.hpp" //Basement membrane force, as defined in Dunn et al. (2012)
#include "EpithelialLayerAnoikisCellKiller.hpp" //Cell killer to remove proliferative cells that have fallen into the lumen

#include "EpithelialLayerDataTrackingModifier.hpp" //Modifier for all the necessary data

#include "FixedRegionPlaneBoundaryCondition.hpp" //Boundary condition that fixes cells past a given plane


/*
 * Define the Chaste simulation as a test class. This is how all simulations
 * in Chaste are defined.
 */
class TestCryptFissionSweepsLiteratePaper : public AbstractCellBasedTestSuite
{
public:
	void TestEpithelialLayerSweeps() throw(Exception)
	{
		/* We first set all the simulation parameters. */

		//Simulation time parameters
		double dt = 0.005; //Set dt
		double end_time = 100.0; //Set end time
		double sampling_timestep = end_time/dt; //Set sampling timestep multiple (we only need the final state for analysis purposes

		//Set all the spring stiffness variables
		double epithelial_epithelial_stiffness = 15.0; //Epithelial-epithelial spring connections
		double epithelial_nonepithelial_stiffness = 15.0; //Epithelial-non-epithelial spring connections
		double nonepithelial_nonepithelial_stiffness = 15.0; //Non-epithelial-non-epithelial spring connections

		//Set the relative adhesiveness of hard cells. The values used in the paper were 1.0, 1.3 and 2.0
		double hard_cell_drag_multiplier = 1.0;

		//Set the BM force parameters
		double bm_force = 10.0; //Set the basement membrane stiffness
		double target_curvature = 0.2; //Set the target curvature, i.e. how circular the layer wants to be

		/* Set the domain of the model. */
		unsigned cells_across = 24; //Desired width + a few more layers, to avoid the box collapsing
		unsigned cells_up = 27; //Since height of each cell is 0.5*sqrt(3), we need to specify the no. of cells such that #cells*0.5*sqrt(3) = desired height
		unsigned ghosts = 4; //Define a sufficient layer of ghost nodes to avoid infinite tessellations and hence excessively large forces

		//Translate mesh 1.5 units left and 1.5 units down, so that we have a sufficient layer of fixed cells around the boundary.
		c_vector<double, 2> translate_left = zero_vector<double>(2);
		translate_left(0) = -1.5;

		c_vector<double, 2> translate_down = zero_vector<double>(2);
		translate_down(1) = -1.5;

		/* Define the initially circular lumen by centre and radius */
		c_vector<double,2> circle_centre;
		circle_centre(0) = 10.5;
		circle_centre(1) = 10.0;

		double circle_radius = 2.5; //Size of hole
		assert(circle_radius > 0); //Just in case someone does something crazy.

		double layer_radius = circle_radius + 2.0; //Radius of the layer of cells. This isn't the actual radius, just has to be large enough for later
		assert((layer_radius <= cells_across)&&(layer_radius <= cells_up)); //Again, just in case.

		/* Generate the initial mesh of cells. */
		HoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
		MutableMesh<2,2>* p_mesh = generator.GetMesh();

		//Translate mesh appropriately
		p_mesh->Translate(translate_left);
		p_mesh->Translate(translate_down);

		/* Define the lumen as an inner region of ghost nodes. */

		std::vector<unsigned> initial_real_indices = generator.GetCellLocationIndices(); //Obtain the locations of real nodes

		std::vector<unsigned> real_indices; //Vector used to define the locations of non-ghost nodes

		//Sweep over the initial real indices
		for (unsigned i = 0; i < initial_real_indices.size(); i++)
		{
			unsigned cell_index = initial_real_indices[i];
			double x = p_mesh->GetNode(cell_index)->rGetLocation()[0];
			double y = p_mesh->GetNode(cell_index)->rGetLocation()[1];

			// If the location of the node falls inside the defined lumen region, then it becomes a ghost node.
			if (pow(x-circle_centre[0],2) + pow(y-circle_centre[1],2) > pow(circle_radius,2))
			{
				real_indices.push_back(cell_index);
			}
		}

		//This is for the purpose of tracking cells. Weird computer science stuff.
		boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_soft_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_hard_state = CellPropertyRegistry::Instance()->Get<HardCellMutationState>();
		boost::shared_ptr<AbstractCellProperty> p_wildtype_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();

		/* Reseed the random number generator via a command line argument */
		unsigned index = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-myintval");
		RandomNumberGenerator::Instance()->Reseed(100*index);

		/* Perform a parameter sweep over the stiffness ratio ond target proportions */
		for (double hstiff = 0.0; hstiff < 9.0; hstiff ++) //Sweep over 3 <= \mu_H/\mu_S <= 5 in increments of 0.25
		{
			for (double tpropn = 0.0; tpropn < 15.0; tpropn ++) //Sweep over 0.25 <= p <= 0.95 in increments of 0.05
			{

				double stiffness_ratio = 3.0 + 0.25*hstiff;
				double target_proportion = 0.25 + 0.05*tpropn;

				//Create vector of cells
				std::vector<CellPtr> cells;

				/*
				 * Create asymmetric-division-based cell cycle for each cell. However, we initially set all cells
		 	 	 * to be non-epithelial cells before defining our layer of epithelial cells.
				 */

				for (unsigned i = 0; i<real_indices.size(); i++)
				{
					//Set cell cycle
					StochasticTargetProportionBasedCellCycleModel* p_cycle_model = new StochasticTargetProportionBasedCellCycleModel();
					p_cycle_model->SetTargetProportion(target_proportion); //Set the division parameter
					p_cycle_model->SetDimension(2);

					//To avoid a 'pulsing' behaviour with birth events, we set each cell's initial age to be
					// ~U(-12, 0) in the past, as each cell cycle duration is U(11, 13).
					double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf();
					p_cycle_model->SetBirthTime(-birth_time);

					CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
					p_cell->SetCellProliferativeType(p_diff_type); //Set the cell to be differentiated and hence non-epithelial
					p_cell->InitialiseCellCycleModel();

					cells.push_back(p_cell);
				}

				//Create cell population
				MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, real_indices);

				//Create layer of proliferative cells
				for (unsigned i = 0; i < real_indices.size(); i++)
				{
					unsigned cell_index = real_indices[i];
					CellPtr cell_iter = cell_population.GetCellUsingLocationIndex(cell_index);
					double x = cell_population.GetLocationOfCellCentre(cell_iter)[0];
					double y = cell_population.GetLocationOfCellCentre(cell_iter)[1];


					/* We 'un-differentiate' any cells adjacent to the lumen into epithelial cells.
					 * We only consider cells within the pre-defined layer radius. Not only does this narrow down
					 * our search, but it also means we won't accidentally consider any cells on the outside and
					 * turn them into epithelial cells.
					 */

					if (pow(x-circle_centre[0],2) + pow(y-circle_centre[1],2) <= pow(layer_radius,2))
					{
						Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(cell_iter);
						unsigned node_index = p_node->GetIndex();

						//Iterate over all possible neighbours of the node
						for (Node<2>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
								iter != p_node->ContainingElementsEnd();
								++iter)
						{
							bool element_contains_ghost_nodes = false;

							// Get a pointer to the element
							Element<2,2>* p_element = cell_population.rGetMesh().GetElement(*iter);

							// Check whether it's triangulation contains a ghost node
							for (unsigned local_index=0; local_index<3; local_index++)
							{
								unsigned nodeGlobalIndex = p_element->GetNodeGlobalIndex(local_index);

								if (cell_population.IsGhostNode(nodeGlobalIndex) == true)
								{
									element_contains_ghost_nodes = true;
									break; 				// This should break out of the inner for loop
								}
							}

							//If a cell has a ghost node as a neighbour, we make it an epithelial cells.
							if(element_contains_ghost_nodes)
							{
								cell_iter->SetCellProliferativeType(p_soft_type);
							}
						}
					}
				}

				/* Iterate again and check that proliferative cells are also attached to non-epithelial
				 * cells. If they are not, remove them from the simulation.
				 */

				for (unsigned i = 0; i < real_indices.size(); i++)
				{
					unsigned cell_index = real_indices[i];
					CellPtr cell_iter = cell_population.GetCellUsingLocationIndex(cell_index);
					double x = cell_population.GetLocationOfCellCentre(cell_iter)[0];
					double y = cell_population.GetLocationOfCellCentre(cell_iter)[1];

					//Only consider this inside the pre-defined layer to narrow down our search

					if (pow(x-circle_centre[0],2) + pow(y-circle_centre[1],2) <= pow(layer_radius,2))
					{
						Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(cell_iter);
						unsigned node_index = p_node->GetIndex();

						//Only iterate over the initial layer of transit cells
						if (cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>() == false)
						{
							bool element_contains_nonepithelial_nodes = false;

							//Iterate over elements (triangles) containing the node
							for (Node<2>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
									iter != p_node->ContainingElementsEnd();
									++iter)
							{
								// Get a pointer to the element (triangle)
								Element<2,2>* p_element = cell_population.rGetMesh().GetElement(*iter);

								// Check if its triangulation contains a gel node
								for (unsigned local_index=0; local_index<3; local_index++)
								{
									unsigned nodeGlobalIndex = p_element->GetNodeGlobalIndex(local_index);
									bool is_ghost_node = cell_population.IsGhostNode(nodeGlobalIndex);

									if (is_ghost_node == false) //Make sure we're not dealing with ghost nodes (otherwise this stuff will fail)
									{
										CellPtr p_local_cell = cell_population.GetCellUsingLocationIndex(nodeGlobalIndex);
										if (p_local_cell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()==true)
										{
											element_contains_nonepithelial_nodes = true;
											break; 				// This should break out of the inner for loop
										}
									}
								}
							}

							if(!element_contains_nonepithelial_nodes)
							{
								cell_iter->Kill();
							}
						}
					}
				}

				/*
				 * Randomly assign cells in the layer to be hard cells.
				 */

				std::vector<unsigned> cells_in_layer; //Initialise vector

				//Obtain the proliferative cells
				for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
						cell_iter != cell_population.End();
						++cell_iter)
				{
					unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);

					//If the cell is an epithelial cell
					if (!cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>() )
					{
						cells_in_layer.push_back(node_index); //Add the angle and node index
					}
				}

				/* For each cell in the layer, we draw a random number and assign cells to be soft cells with
				 * a probability equal to the target proportion, as defined above.
				 */

				for (unsigned i = 0; i < cells_in_layer.size(); i++)
				{
					unsigned node_index = cells_in_layer[i];

					CellPtr cell = cell_population.GetCellUsingLocationIndex(node_index);

					//Randomly generate number
					double random_number = RandomNumberGenerator::Instance()->ranf();

					if(random_number >= target_proportion) //Assign cells to be Paneth with 1 - target_proportion
					{
						cell->SetMutationState(p_hard_state);
					}
				}

				/* Define the simulation class. */
				OffLatticeSimulation<2> simulator(cell_population);

				/* Set the drag constant of hard cells relative to soft cells */
				double normal_damping_constant = cell_population.GetDampingConstantNormal();
				cell_population.SetDampingConstantMutant(hard_cell_drag_multiplier*normal_damping_constant);

				//Simulation pre-amble

				//Set output directory
				std::stringstream out;
				out << "Drag_" << hard_cell_drag_multiplier;
				out << "/BM_" << bm_force << "_TC_" << target_curvature;
				std::stringstream params;
				params << "SRatio_" << stiffness_ratio << "/TPropn_" << target_proportion;
				std::stringstream run;
				run << index;
				std::string output_directory = "CryptFissionSweepsLiteratePaper/" + out.str() + "/"
						+ run.str() + "/" + params.str();
				simulator.SetOutputDirectory(output_directory);

				simulator.SetDt(dt);
				simulator.SetSamplingTimestepMultiple(sampling_timestep); //Sample the simulation at every hour
				simulator.SetEndTime(end_time); //Hopefully this is long enough for a steady state

				/* We add a modifier class to track relevant cell population numbers and shape measurements.*/
				MAKE_PTR(EpithelialLayerDataTrackingModifier<2>, p_data_tracking_modifier);
				simulator.AddSimulationModifier(p_data_tracking_modifier);

				/* Add linear spring force which has different spring stiffness constants, depending
				 * on the pair of cells it is connecting.
				 */
				MAKE_PTR(EpithelialLayerLinearSpringForce<2>, p_spring_force);
				p_spring_force->SetCutOffLength(1.5);
				//Set the spring stiffnesses
				p_spring_force->SetEpithelialEpithelialSpringStiffness(epithelial_epithelial_stiffness);
				p_spring_force->SetEpithelialNonepithelialSpringStiffness(epithelial_nonepithelial_stiffness);
				p_spring_force->SetNonepithelialNonepithelialSpringStiffness(nonepithelial_nonepithelial_stiffness);
				p_spring_force->SetHardCellStiffnessRatio(stiffness_ratio);
				simulator.AddForce(p_spring_force);

				/* Add the basement membrane force. */
				MAKE_PTR(EpithelialLayerBasementMembraneForce, p_bm_force);
				p_bm_force->SetBasementMembraneParameter(bm_force); //Equivalent to beta in SJD's papers
				p_bm_force->SetTargetCurvature(target_curvature); //This is equivalent to 1/R in SJD's papers
				simulator.AddForce(p_bm_force);

				/* Add an anoikis-based cell killer. */
				MAKE_PTR_ARGS(EpithelialLayerAnoikisCellKiller, p_anoikis_killer, (&cell_population));
				simulator.AddCellKiller(p_anoikis_killer);

				/* We fix all cells outside of the 20 x 20 box. */

				c_vector<double,2> point = zero_vector<double>(2);
				c_vector<double,2> normal = zero_vector<double>(2);

				//Fix cells in the region x < 0
				normal(0) = -1.0;
				MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
				simulator.AddCellPopulationBoundaryCondition(p_bc1);

				//Fix cells in the region x > 20
				point(0) = 20.0;
				normal(0) = 1.0;
				MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
				simulator.AddCellPopulationBoundaryCondition(p_bc2);

				//Fix cells in the region y < 0
				point(0) = 0.0;
				point(1) = 0.0;
				normal(0) = 0.0;
				normal(1) = -1.0;
				MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
				simulator.AddCellPopulationBoundaryCondition(p_bc3);

				//Fix cells in the region y > 20
				point(1) = 20.0;
				normal(1) = 1.0;
				MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal));
				simulator.AddCellPopulationBoundaryCondition(p_bc4);

				/* Run the simulation. */
				simulator.Solve();

				/* We now collect statistics on the positions and areas of the cells in the epithelial layer that are
				 * used for subsequent analysis.
				 *
				 * We first need the epithelial cells ordered around the layer.
				 */

				std::vector<unsigned> indices_of_cells_in_layer; //Initialise vector

				c_vector<double,2> centre_of_mass; //Initialise centre of mass

				unsigned num_epithelial_cells = 0; //Initialise count

				for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
						cell_iter != simulator.rGetCellPopulation().End(); ++cell_iter)
				{
					//Only consider epithelial cells
					if(!cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>() )
					{
						//Get node inex
						unsigned node_index = simulator.rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter);

						//Get location of cell
						c_vector<double,2> cell_location = simulator.rGetCellPopulation().GetLocationOfCellCentre(*cell_iter);

						indices_of_cells_in_layer.push_back(node_index); //Add node index

						centre_of_mass += cell_location; //Update the centre of mass

						num_epithelial_cells += 1; //Update count of epithelial cells in layer
					}
				}

				//Normalise centre_of_mass
				centre_of_mass /= num_epithelial_cells;

				//Initialise vector of angles and the cell indices
				std::vector<std::pair<double, unsigned> > angles_and_cell_indices;

				/* Get the angle of each epithelial cell with respect to the centre of mass.
				 * That way we can sort the list in counter-clockwise order before we create
				 *  the file.
				 */
				for (unsigned k = 0; k < indices_of_cells_in_layer.size(); k++)
				{
					unsigned cell_index = indices_of_cells_in_layer[k]; //Get the node index

					CellPtr cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(cell_index); //Get the cell

					//Get the cell location
					double x = simulator.rGetCellPopulation().GetLocationOfCellCentre(cell)[0];
					double y = simulator.rGetCellPopulation().GetLocationOfCellCentre(cell)[1];

					//Make the pair of the angle and cell
					std::pair<double, unsigned> angle_cell;

					//Get point relative to the circle centre
					double rel_x = x - centre_of_mass[0];
					double rel_y = y - centre_of_mass[1];

					double circle_angle = atan(rel_y/rel_x); //Get initial angle argument

					if (rel_x<0.0) //If the point is in the second quadrant or third quadrant
					{
						circle_angle += M_PI;
					}
					else if ((rel_x>=0.0)&&(rel_y<0.0)) //Fourth quadrant
					{
						circle_angle += 2*M_PI;
					}

					angle_cell = std::make_pair(circle_angle, cell_index);

					angles_and_cell_indices.push_back(angle_cell); //Add the angle and node index
				}

				//Sort the vector of cell indices by their angle
				std::sort(angles_and_cell_indices.begin(), angles_and_cell_indices.end());

				/* Create file to store the location of the cells in the layer and their cell type. */
				OutputFileHandler results_handler(output_directory + "/", false); //Create output file handler

				std::string cells_in_layer_filename = "cellpositionsandareasinlayer.dat";
				out_stream cells_in_layer_file = results_handler.OpenOutputFile(cells_in_layer_filename);

				/* Add each cell location and area to the file, along with a 1 if the cell is
				 * a hard cell and 0 otherwise.
				 */
				for (unsigned k=0; k<angles_and_cell_indices.size(); k++)
				{
					std::pair<double, unsigned> angle_and_cell_index = angles_and_cell_indices[k]; //Get pair
					unsigned cell_index = angle_and_cell_index.second; //Get cell index

					CellPtr cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(cell_index); //Get the cell

					//Get the area
					double cell_area = simulator.rGetCellPopulation().GetVolumeOfCell(cell);

					//Get the cell location
					double x = simulator.rGetCellPopulation().GetLocationOfCellCentre(cell)[0];
					double y = simulator.rGetCellPopulation().GetLocationOfCellCentre(cell)[1];

					unsigned is_hard_cell = 1;

					if(!cell->GetMutationState()->IsType<HardCellMutationState>() )
					{
						is_hard_cell = 0;
					}

					*cells_in_layer_file << x << "\t" << y << "\t" << cell_area << is_hard_cell << "\n";
				}

				//Tidying up
				SimulationTime::Instance()->Destroy();
				SimulationTime::Instance()->SetStartTime(0.0);
			}
		}
	}

};

#endif /* TESTCRYPTFISSIONSWEEPSLITERATEPAPER_HPP_ */
