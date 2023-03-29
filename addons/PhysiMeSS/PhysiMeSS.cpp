#include "PhysiMeSS.h"

#include <algorithm>
#include <iterator> 

// !!! PHYSIMESS CODE BLOCK START !!! //
std::map<Cell*, std::vector<Cell*> > fibres_crosslinkers;
std::map<Cell*, std::vector<double> > fibres_crosslink_point;
    
void initialize_crosslinkers(Cell* pCell) { fibres_crosslinkers[pCell] = std::vector<Cell*>(); }
void initialize_crosslink_points(Cell* pCell) { fibres_crosslink_point[pCell] = std::vector<double>(3); }
std::vector<Cell*>& get_crosslinkers(Cell* pCell) { return fibres_crosslinkers[pCell]; }
std::vector<double>& get_crosslink_point(Cell* pCell) { return fibres_crosslink_point[pCell]; }

// void Fibre::degrade_fibre(Cell *fibre_to_degrade) {

//     // don't degrade anything that is not a fibre
//     if ((*fibre_to_degrade).type_name != "fibre") { return; }
//     // de-allocate (delete) the cell;
//     fibre_to_degrade->flag_for_removal();
//     return;
// }

std::vector<double> nearest_point_on_fibre(std::vector<double> point, Cell *fibre_agent, std::vector<double> &displacement) 
{

    // don't bother if the "fibre_agent" is not a fibre
    if ((*fibre_agent).type_name != "fibre") { return displacement; }

    double fibre_length = 2 * (*fibre_agent).custom_data["mLength"];
    // vector pointing from one endpoint of "fibre_agent" to "point"
    std::vector<double> fibre_to_agent(3, 0.0);
    // |fibre_to_agent| squared
    double fibre_to_agent_length_squared = 0;
    // scalar product fibre_to_agent * fibre_vector
    double fibre_to_agent_dot_fibre_vector = 0;

    double distance = 0;
    for (unsigned int i = 0; i < 3; i++) {
        fibre_to_agent[i] = point[i] - ((*fibre_agent).position[i]
                                        - (*fibre_agent).custom_data["mLength"] * (*fibre_agent).state.orientation[i]);
        fibre_to_agent_length_squared += fibre_to_agent[i] * fibre_to_agent[i];
        fibre_to_agent_dot_fibre_vector += fibre_to_agent[i] * fibre_length * (*fibre_agent).state.orientation[i];
    }

    // "point" is closest to the selected endpoint of "fibre_agent"
    if (fibre_to_agent_dot_fibre_vector < 0.) {
        for (int i = 0; i < 3; i++) {
            displacement[i] = fibre_to_agent[i];
        }
        //std::cout << "The point is closest to the start of the fibre" << std::endl;
        //std::cout << " Displacement: " << displacement << std::endl;
    }
        // “point” is closest to the other endpoint of “fibre_agent”
    else if (fibre_to_agent_dot_fibre_vector > fibre_length * fibre_length) {
        for (unsigned int i = 0; i < 3; i++) {
            displacement[i] = point[i] - ((*fibre_agent).position[i]
                                            + (*fibre_agent).custom_data["mLength"] * (*fibre_agent).state.orientation[i]);
        }
        //std::cout << "The point is closest to the end of the fibre" << std::endl;
        //std::cout << " Displacement: " << displacement << std::endl;
    }
        // “point” is closest to a point along “fibre_agent”
    else {
        double fibre_to_agent_length_cos_alpha_squared =
                fibre_to_agent_dot_fibre_vector * fibre_to_agent_dot_fibre_vector /
                (fibre_length * fibre_length);
        double l = sqrt(fibre_to_agent_length_cos_alpha_squared);
        for (unsigned int i = 0; i < 3; i++) {
            displacement[i] = fibre_to_agent[i] - l * (*fibre_agent).state.orientation[i];
        }
        //std::cout << "The point is closest to a point along the fibre" << std::endl;
        //std::cout << " Displacement: " << displacement << std::endl;
    }

    // the function returns the displacement vector
    return displacement;
}

void check_fibre_crosslinks(Cell* fibre, Cell *fibre_neighbor) {

    if (fibre == fibre_neighbor) { return; }

    if (fibre->type_name == "fibre" && (*fibre_neighbor).type_name == "fibre") {

        // fibre endpoints
        std::vector<double> point1(3, 0.0);
        std::vector<double> point2(3, 0.0);
        std::vector<double> point3(3, 0.0);
        std::vector<double> point4(3, 0.0);
        for (int i = 0; i < 3; i++) {
            // endpoints of "this" fibre
            point1[i] = fibre->position[i]
                        - fibre->custom_data["mLength"] * fibre->state.orientation[i];
            point2[i] = fibre->position[i]
                        + fibre->custom_data["mLength"] * fibre->state.orientation[i];
            // endpoints of "neighbor" fibre
            point3[i] = (*fibre_neighbor).position[i]
                        - fibre->custom_data["mLength"] * (*fibre_neighbor).state.orientation[i];
            point4[i] = (*fibre_neighbor).position[i]
                        + fibre->custom_data["mLength"] * (*fibre_neighbor).state.orientation[i];
        }

        //vectors between fibre endpoints
        std::vector<double> p1_to_p2(3, 0.0);
        std::vector<double> p3_to_p4(3, 0.0);
        std::vector<double> p1_to_p3(3, 0.0);
        std::vector<double> centre_to_centre(3, 0.0);
        for (int i = 0; i < 3; i++) {
            // "fibre" fibre vector
            p1_to_p2[i] = point2[i] - point1[i];
            // "neighbor" fibre vector
            p3_to_p4[i] = point4[i] - point3[i];
            // vector from "fibre" to "neighbor"
            p1_to_p3[i] = point3[i] - point1[i];
            // vector between fibre centres
            centre_to_centre[i] = (*fibre_neighbor).position[i] - fibre->position[i];
        }

        double co_radius = fibre->custom_data["mRadius"] + (*fibre_neighbor).custom_data["mRadius"];
        double co_length = fibre->custom_data["mLength"] + (*fibre_neighbor).custom_data["mLength"];
        std::vector<double> zero(3, 0.0);
        double distance = dist(zero, centre_to_centre);
        normalize(&centre_to_centre);

        /* test if fibres intersect
            (1) if fibres are coplanar and parallel:
            the cross product of the two fibre vectors is zero
            [(P2 - P1) x (P4 - P3)].[(P2 - P1) x (P4 - P3)] = 0 */
        std::vector<double> FCP(3, 0.0);
        CrossProduct(p1_to_p2, p3_to_p4, FCP);
        /* parallel fibres could intersect if colinear
            we need to take into account the length of fibres
            possibly could put a tolerance on this to take into
            account the radius of fibres */
        if (DotProduct(FCP,FCP) == 0 &&
            (centre_to_centre == fibre->state.orientation ||
            centre_to_centre == -1.0 * fibre->state.orientation) &&
            distance <= co_length) {
            //std::cout << "fibre " << fibre->ID << " crosslinks end to end with parallel fibre " <<  (*fibre_neighbor).ID << std::endl;
            if (std::find(fibres_crosslinkers[fibre].begin(), fibres_crosslinkers[fibre].end(), (fibre_neighbor)) ==
                fibres_crosslinkers[fibre].end()) {
                fibres_crosslinkers[fibre].push_back(fibre_neighbor);
            }
            if (std::find(fibres_crosslinkers[fibre_neighbor].begin(), fibres_crosslinkers[fibre_neighbor].end(), (fibre)) ==
                fibres_crosslinkers[fibre_neighbor].end()) {
                fibres_crosslinkers[fibre_neighbor].push_back(fibre);
            }
            fibres_crosslink_point[fibre] = fibre->position + fibre->custom_data["mLength"] * centre_to_centre;
            fibres_crosslink_point[fibre_neighbor] = fibres_crosslink_point[fibre];
        }
        std::vector<double> displacement(3, 0.0);
        nearest_point_on_fibre(point1, fibre_neighbor, displacement);
        double test_point1 = dist(zero, displacement);
        nearest_point_on_fibre(point2, fibre_neighbor, displacement);
        double test_point2 = dist(zero, displacement);
        if (test_point1 <= co_radius &&
            centre_to_centre != fibre->state.orientation &&
            centre_to_centre != -1.0 * fibre->state.orientation) {
            //std::cout << "fibre " << fibre->ID << " crosslinks in parallel plane with fibre " <<  (*fibre_neighbor).ID << std::endl;
            if (std::find(fibres_crosslinkers[fibre].begin(), fibres_crosslinkers[fibre].end(), (fibre_neighbor)) ==
                fibres_crosslinkers[fibre].end()) {
                fibres_crosslinkers[fibre].push_back(fibre_neighbor);
            }
            if (std::find(fibres_crosslinkers[fibre_neighbor].begin(), fibres_crosslinkers[fibre_neighbor].end(), (fibre)) ==
                fibres_crosslinkers[fibre_neighbor].end()) {
                fibres_crosslinkers[fibre_neighbor].push_back(fibre);
            }
            fibres_crosslink_point[fibre] = point1;
            fibres_crosslink_point[fibre_neighbor] = fibres_crosslink_point[fibre];
        }
        if (test_point2 <= co_radius &&
            centre_to_centre != fibre->state.orientation &&
            centre_to_centre != -1.0 * fibre->state.orientation) {
            //std::cout << "fibre " << fibre->ID << " crosslinks in parallel plane with fibre " <<  (*fibre_neighbor).ID << std::endl;
            if (std::find(fibres_crosslinkers[fibre].begin(), fibres_crosslinkers[fibre].end(), (fibre_neighbor)) ==
                fibres_crosslinkers[fibre].end()) {
                fibres_crosslinkers[fibre].push_back(fibre_neighbor);
            }
            if (std::find(fibres_crosslinkers[fibre_neighbor].begin(), fibres_crosslinkers[fibre_neighbor].end(), (fibre)) ==
                fibres_crosslinkers[fibre_neighbor].end()) {
                fibres_crosslinkers[fibre_neighbor].push_back(fibre);
            }
            fibres_crosslink_point[fibre] = point2;
            fibres_crosslink_point[fibre_neighbor] = fibres_crosslink_point[fibre];
        }
        /*  (2) if fibres are skew (in parallel planes):
            the scalar triple product (P3 - P1) . [(P2 - P1) x (P4 - P3)] != 0
            we include a tolerance on this to allow for fibre radius */
        double test2_tolerance = co_radius;
        double test2 = DotProduct(p1_to_p3, FCP);
        /*  For fibres to intersect neither (1) not (2) must be true
            and the "t" values for both line equations must lie in [0,1]
            again we include a tolerance to allow for fibre radius */
        if (std::abs(test2) < test2_tolerance) {
            double a = DotProduct(p1_to_p2, p1_to_p3) / DotProduct(p1_to_p2, p1_to_p2);
            double b = DotProduct(p1_to_p2, p3_to_p4) / DotProduct(p1_to_p2, p1_to_p2);
            std::vector<double> c(3, 0.0);
            std::vector<double> n(3, 0.0);
            for (int i = 0; i < 3; i++) {
                c[i] = b * p1_to_p2[i] - p3_to_p4[i];
                n[i] = p1_to_p3[i] - a * p1_to_p2[i];
            }
            double t_neighbor = DotProduct(c, n) / DotProduct(c, c);
            double t_this = a + b * t_neighbor;
            std::vector<double> crosslink_point(3, 0.0);
            for (int i = 0; i < 2; i++) {
                crosslink_point[i] = point1[i] + t_this * p1_to_p2[i];
            }
            double tolerance = (*fibre_neighbor).custom_data["mRadius"] / (2 * fibre->custom_data["mLength"]);
            double lower_bound = 0.0 - tolerance;
            double upper_bound = 1.0 + tolerance;
            if (lower_bound <= t_neighbor && t_neighbor <= upper_bound &&
                lower_bound <= t_this && t_this <= upper_bound) {
                // if (std::find(fibre->state.crosslinkers.begin(), fibre->state.crosslinkers.end(), (fibre_neighbor)) ==
                //     fibre->state.crosslinkers.end()) {
                //     //std::cout << "fibre " << fibre->ID << " crosslinks with fibre " << (*fibre_neighbor).ID << std::endl;
                //     fibre->state.crosslinkers.push_back(fibre_neighbor);
                // }
                // if (std::find(fibre_neighbor->state.crosslinkers.begin(), fibre_neighbor->state.crosslinkers.end(), (fibre)) ==
                //     fibre_neighbor->state.crosslinkers.end()) {
                //     fibre_neighbor->state.crosslinkers.push_back(fibre);
                // }
                if (std::find(fibres_crosslinkers[fibre].begin(), fibres_crosslinkers[fibre].end(), (fibre_neighbor)) ==
                    fibres_crosslinkers[fibre].end()) {
                    //std::cout << "fibre " << fibre->ID << " crosslinks with fibre " << (*fibre_neighbor).ID << std::endl;
                    fibres_crosslinkers[fibre].push_back(fibre_neighbor);
                }
                if (std::find(fibres_crosslinkers[fibre_neighbor].begin(), fibres_crosslinkers[fibre_neighbor].end(), (fibre)) ==
                    fibres_crosslinkers[fibre_neighbor].end()) {
                    fibres_crosslinkers[fibre_neighbor].push_back(fibre);
                }
                // fibre->state.crosslink_point = crosslink_point;
                // fibre_neighbor->state.crosslink_point = crosslink_point;
                fibres_crosslink_point[fibre] = crosslink_point;
                fibres_crosslink_point[fibre_neighbor] = crosslink_point;
            }
        }
    } else { return; }
}


// !!! PHYSIMESS CODE BLOCK START !!! //
void force_update_motility_vector(Cell* pCell, double dt_) {

    if (!pCell->phenotype.motility.is_motile) {
        pCell->phenotype.motility.motility_vector.assign(3, 0.0);
        return;
    }

    // force cell to update its motility because it is stuck
    // choose a uniformly random unit vector
    double temp_angle = 6.28318530717959 * UniformRandom();
    double temp_phi = 3.1415926535897932384626433832795 * UniformRandom();

    double sin_phi = sin(temp_phi);
    double cos_phi = cos(temp_phi);

    if (pCell->phenotype.motility.restrict_to_2D) {
        sin_phi = 1.0;
        cos_phi = 0.0;
    }

    std::vector<double> randvec;
    randvec.resize(3, sin_phi);

    randvec[0] *= cos(temp_angle); // cos(theta)*sin(phi)
    randvec[1] *= sin(temp_angle); // sin(theta)*sin(phi)
    randvec[2] = cos_phi; //  cos(phi)

    // if the update_bias_vector function is set, use it
    /*if (functions.update_migration_bias) {
    functions.update_migration_bias(this, phenotype, dt_);
    }*/

    //phenotype.motility.motility_vector *= -1.0;//phenotype.motility.migration_bias_direction; // motiltiy = bias_vector
    //phenotype.motility.motility_vector *= phenotype.motility.migration_bias; // motility = bias*bias_vector

    double one_minus_bias = 1.0;// - phenotype.motility.migration_bias;

    axpy(&(pCell->phenotype.motility.motility_vector), one_minus_bias,randvec); // motility = (1-bias)*randvec + bias*bias_vector

    normalize(&(pCell->phenotype.motility.motility_vector));

    pCell->phenotype.motility.motility_vector *= pCell->phenotype.motility.migration_speed;

    return;
}

void add_potentials_cell_to_fibre(Cell* pCell, Cell* other_agent)
{
    double distance = 0.0;
    nearest_point_on_fibre(pCell->position, other_agent, pCell->displacement);
    for (int index = 0; index < 3; index++) {
        distance += pCell->displacement[index] * pCell->displacement[index];
    }
    distance = std::max(sqrt(distance), 0.00001);
    /*if( this->phenotype.motility.is_motile) {
        std::cout << " determining distance from " << this->type_name << " " << this->ID << " to "
                    << (*other_agent).type_name << " " << (*other_agent).ID
                    << "   the distance is " << distance << std::endl;
    }*/

    // as per PhysiCell
    static double simple_pressure_scale = 0.027288820670331;

    // check distance relative repulsion and adhesion distances
    // cell should repel from a fibre if it comes within cell radius plus fibre radius (note fibre radius ~2 micron)
    double R = pCell->phenotype.geometry.radius + (*other_agent).custom_data["mRadius"];
    // cell should feel adhesion over
    double max_interactive_distance =
            pCell->phenotype.mechanics.relative_maximum_adhesion_distance * pCell->phenotype.geometry.radius +
            (*other_agent).phenotype.mechanics.relative_maximum_adhesion_distance *
            (*other_agent).custom_data["mRadius"];

    // First Repulsion as per PhysiCell
    double temp_r = 0;
    if (distance > R) {
        temp_r = 0;
    } else {
        // temp_r = 1 - distance/R;
        temp_r = -distance; // -d
        temp_r /= R; // -d/R
        temp_r += 1.0; // 1-d/R
        temp_r *= temp_r; // (1-d/R)^2

        // add the relative pressure contribution NOT SURE IF NEEDED
        pCell->state.simple_pressure += (temp_r / simple_pressure_scale);

        double effective_repulsion = sqrt(pCell->phenotype.mechanics.cell_cell_repulsion_strength *
                                            (*other_agent).phenotype.mechanics.cell_cell_repulsion_strength);
        temp_r *= effective_repulsion;
    }

    if (fabs(temp_r) < 1e-16) { return; }
    temp_r /= distance;

    axpy(&(pCell->velocity), temp_r, pCell->displacement);

    //Then additional repulsion/adhesion as per Cicely's code
    double fibre_adhesion = 0;
    double fibre_repulsion = 0;
    if (distance < max_interactive_distance) {
        double cell_velocity_dot_fibre_direction = 0.;
        for (unsigned int j = 0; j < 3; j++) {
            cell_velocity_dot_fibre_direction += (*other_agent).state.orientation[j] * pCell->previous_velocity[j];
        }
        double cell_velocity = 0;
        for (unsigned int j = 0; j < pCell->velocity.size(); j++) {
            cell_velocity += pCell->previous_velocity[j] * pCell->previous_velocity[j];
        }
        cell_velocity = std::max(sqrt(cell_velocity), 1e-8);

        double p_exponent = 1.0;
        double q_exponent = 1.0;
        double xi = fabs(cell_velocity_dot_fibre_direction) / (cell_velocity);
        double xip = pow(xi, p_exponent);
        double xiq = pow((1 - xi * xi), q_exponent);

        fibre_adhesion = parameters.doubles("vel_adhesion") * xip *
                            (1 - cell_velocity / parameters.doubles("cell_velocity_max"));

        fibre_repulsion = parameters.doubles("vel_contact") * xiq;

        axpy(&(pCell->velocity), fibre_adhesion, (*other_agent).state.orientation);
        naxpy(&(pCell->velocity), fibre_repulsion, pCell->previous_velocity);

        // Fibre degradation by cell - switched on by flag fibre_degradation
        int stuck_threshold = parameters.doubles("fibre_stuck");
        if (parameters.bools("fibre_degradation") && pCell->custom_data["stuck_counter"] >= stuck_threshold) {
            //std::cout << "Cell " << this->ID << " is stuck at time " << PhysiCell_globals.current_time
            //<< " near fibre " << (*other_agent).ID  << std::endl;
            pCell->displacement *= -1.0/distance;
            double dot_product = DotProduct(pCell->displacement, pCell->phenotype.motility.motility_vector);
            if (dot_product >= 0) {
                double rand_degradation = UniformRandom();
                double prob_degradation = parameters.doubles("fibre_deg_rate");
                if (rand_degradation <= prob_degradation) {
                    //std::cout << " --------> fibre " << (*other_agent).ID << " is flagged for degradation " << std::endl;
                    // (*other_agent).parameters.degradation_flag = true;
                    other_agent->flag_for_removal();
                    std::cout << "Killing a cell using flag for removal !!" << std::endl;
                    pCell->custom_data["stuck_counter"] = 0;
                }
            }
        }
    }

}

void add_potentials_fibre_to_cell(Cell* pCell, Cell* other_agent) 
{
    // fibres only get pushed or rotated by motile cells
    if (!other_agent->phenotype.motility.is_motile || pCell->custom_data["X_crosslink_count"] >= 2) {
        return;
    }

    double distance = 0.0;
    nearest_point_on_fibre((*other_agent).position, pCell, pCell->displacement);
    for (int index = 0; index < 3; index++) {
        distance += pCell->displacement[index] * pCell->displacement[index];
    }
    distance = std::max(sqrt(distance), 0.00001);
    // fibre should only interact with cell if it comes within cell radius plus fibre radius (note fibre radius ~2 micron)
    double R = pCell->phenotype.geometry.radius + pCell->custom_data["mRadius"];
    if (distance <= R) {
        std::vector<double> point_of_impact(3, 0.0);
        for (int index = 0; index < 3; index++) {
            point_of_impact[index] = (*other_agent).position[index] - pCell->displacement[index];
        }
        // cell-fibre pushing only if fibre no crosslinks
        if (pCell->custom_data["X_crosslink_count"] == 0) {
            //fibre pushing turned on
            if (parameters.bools("fibre_pushing")) {
                // as per PhysiCell
                static double simple_pressure_scale = 0.027288820670331;
                // temp_r = 1 - distance/R;
                double temp_r = 0;
                temp_r = -distance;
                temp_r /= R;
                temp_r += 1.0;
                temp_r *= temp_r;
                // add the relative pressure contribution NOT SURE IF NEEDED
                pCell->state.simple_pressure += (temp_r / simple_pressure_scale);

                double effective_repulsion = sqrt(pCell->phenotype.mechanics.cell_cell_repulsion_strength *
                                                    (*other_agent).phenotype.mechanics.cell_cell_repulsion_strength);
                temp_r *= effective_repulsion;

                if (fabs(temp_r) < 1e-16) { return; }
                temp_r /= distance;
                naxpy(&pCell->velocity, temp_r, pCell->displacement);
            }

            // fibre rotation turned on (2D)
            if (parameters.bools("fibre_rotation")) {
                std::vector<double> old_orientation(3, 0.0);
                for (int i = 0; i < 2; i++) {
                    old_orientation[i] = pCell->state.orientation[i];
                }

                double moment_arm_magnitude = sqrt(
                        point_of_impact[0] * point_of_impact[0] + point_of_impact[1] * point_of_impact[1]);
                double impulse = parameters.doubles("mFibreStickiness")*(*other_agent).phenotype.motility.migration_speed * moment_arm_magnitude;
                double fibre_length = 2 * pCell->custom_data["mLength"];
                double angular_velocity = impulse / (0.5 * fibre_length * fibre_length);
                double angle = angular_velocity;
                pCell->state.orientation[0] = old_orientation[0] * cos(angle) - old_orientation[1] * sin(angle);
                pCell->state.orientation[1] = old_orientation[0] * sin(angle) + old_orientation[1] * cos(angle);
                normalize(&pCell->state.orientation);
            }
        }

        // fibre rotation around other fibre (2D only and fibres intersect at a single point)
        if (parameters.bools("fibre_rotation") && pCell->custom_data["X_crosslink_count"] == 1) {
            double distance_fibre_centre_to_crosslink = 0.0;
            std::vector<double> fibre_centre_to_crosslink(3, 0.0);
            for (int i = 0; i < 2; i++) {
                fibre_centre_to_crosslink[i] = fibres_crosslink_point[pCell][i]-pCell->position[i];
                distance_fibre_centre_to_crosslink += fibre_centre_to_crosslink[i]*fibre_centre_to_crosslink[i];
            }
            distance_fibre_centre_to_crosslink = sqrt(distance_fibre_centre_to_crosslink);

            std::vector<double> old_orientation(3, 0.0);
            for (int i = 0; i < 2; i++) {
                old_orientation[i] = pCell->state.orientation[i];
            }
            double moment_arm_magnitude = sqrt(
                    point_of_impact[0] * point_of_impact[0] + point_of_impact[1] * point_of_impact[1]);
            double impulse = parameters.doubles("mFibreStickiness")*(*other_agent).phenotype.motility.migration_speed * moment_arm_magnitude;
            double fibre_length = 2 * pCell->custom_data["mLength"];
            double angular_velocity = impulse / (0.5 * fibre_length * fibre_length);
            double angle = angular_velocity;
            pCell->state.orientation[0] = old_orientation[0] * cos(angle) - old_orientation[1] * sin(angle);
            pCell->state.orientation[1] = old_orientation[0] * sin(angle) + old_orientation[1] * cos(angle);
            normalize(&pCell->state.orientation);
            pCell->position[0] = fibres_crosslink_point[pCell][0]-distance_fibre_centre_to_crosslink*pCell->state.orientation[0];
            pCell->position[1] = fibres_crosslink_point[pCell][1]-distance_fibre_centre_to_crosslink*pCell->state.orientation[1];
        }
    }

    return;
}
void add_potentials_fibre_to_fibre(Cell* pCell, Cell* other_agent) 
{
    /* probably want something here to model tension along fibres
        * this will be strong tension along the fibre for those fibres with a crosslink
        * and weak tension with background ECM */
    return;
}
    
// !!! PHYSIMESS CODE BLOCK END !!! //


void add_crosslinks( Cell* pCell) {

    if (pCell->type_name != "fibre") { return; }

    // Determine all crosslinks between pCell and its neighbors provided both are fibres
    else if (pCell->type_name == "fibre") {
        std::vector<Cell *>::iterator f_neighbor;
        std::vector<Cell *>::iterator f_end = pCell->state.neighbors.end();
        for (f_neighbor = pCell->state.neighbors.begin(); f_neighbor != f_end; ++f_neighbor) {
            if ((*f_neighbor)->type_name == "fibre") {
                check_fibre_crosslinks(pCell, *f_neighbor);
            }
        }
    }
}
// !!! PHYSIMESS CODE BLOCK END !!! //

// !!! PHYSIMESS CODE BLOCK START !!! //
std::list<int> register_fibre_voxels(Cell *pCell) {

    std::list<int> agent_voxels;
    int voxel;

    //only do this for fibres
    if (pCell->type_name != "fibre") {
        voxel = pCell->get_container()->underlying_mesh.nearest_voxel_index(pCell->position);
        agent_voxels.push_back(voxel);
    }

    if (pCell->type_name == "fibre") {
        int voxel_size = 30; // note this must be the same as the mechanics_voxel_size
        int test = 2.0 * pCell->custom_data["mLength"] / voxel_size; //allows us to sample along the fibre

        std::vector<double> fibre_start(3, 0.0);
        std::vector<double> fibre_end(3, 0.0);
        for (unsigned int i = 0; i < 3; i++) {
            fibre_start[i] = pCell->position[i] - pCell->custom_data["mLength"] * pCell->state.orientation[i];
            fibre_end[i] = pCell->position[i] + pCell->custom_data["mLength"] * pCell->state.orientation[i];
        }
        // first add the voxel of the fibre end point
        voxel = pCell->get_container()->underlying_mesh.nearest_voxel_index(fibre_end);
        agent_voxels.push_back(voxel);
        if (std::find(pCell->get_container()->agent_grid[voxel].begin(),
                        pCell->get_container()->agent_grid[voxel].end(),
                        pCell) == pCell->get_container()->agent_grid[voxel].end()) {
            pCell->get_container()->agent_grid[voxel].push_back(pCell);
        }
        // then walk along the fibre from fibre start point sampling and adding voxels as we go
        std::vector<double> point_on_fibre(3, 0.0);
        for (unsigned int j = 0; j < test + 1; j++) {
            for (unsigned int i = 0; i < 3; i++) {
                point_on_fibre[i] = fibre_start[i] + j * voxel_size * pCell->state.orientation[i];
            }
            voxel = pCell->get_container()->underlying_mesh.nearest_voxel_index(point_on_fibre);
            agent_voxels.push_back(voxel);
            if (std::find(pCell->get_container()->agent_grid[voxel].begin(),
                            pCell->get_container()->agent_grid[voxel].end(),
                            pCell) == pCell->get_container()->agent_grid[voxel].end()) {
                pCell->get_container()->agent_grid[voxel].push_back(pCell);
            }
        }

        agent_voxels.sort();
        agent_voxels.unique();
    }

    return agent_voxels;
}

void deregister_fibre_voxels(Cell *pCell) {

    //only do this for fibres
    if (pCell->type_name != "fibre") { return;}

    if (pCell->type_name == "fibre") {
        std::list<int> voxels = register_fibre_voxels(pCell);
        int centre_voxel =
                (*pCell).get_container()->underlying_mesh.nearest_voxel_index(pCell->position);
        for (int x: voxels) {
            if (x != centre_voxel) {
                (*pCell).get_container()->remove_agent_from_voxel(pCell, x);
            }
        }
    }
}

std::list<int> find_agent_voxels(Cell *pCell) {

    /* this code is for creating a list of all voxels which either contain the agent
        * or are neighboring voxels of the voxel containing the agent */
    std::list<int> all_agent_voxels_to_test;
    std::list<int> voxels = register_fibre_voxels(pCell);
    for (int x: voxels) {
        all_agent_voxels_to_test.push_back(x);
        std::vector<int>::iterator xx;
        std::vector<int>::iterator x_end =
                pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[x].end();
        for (xx = pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[x].begin();
                xx != x_end; ++xx) {
            all_agent_voxels_to_test.push_back(*xx);
        }
    }
    // get rid of any duplicated voxels
    all_agent_voxels_to_test.sort();
    all_agent_voxels_to_test.unique();

    return all_agent_voxels_to_test;
}

void find_agent_neighbors(Cell *pCell) {

    /* this code is for finding all neighbors of an agent: first we call find_agent_voxels
        *  to create a list of all the voxels to test, then we search for agents in those voxels */
    std::list<int> voxels_to_test = find_agent_voxels(pCell);

    //std::cout << "Agent " << pCell->ID << " is tested in voxels: " ;
    for (int x: voxels_to_test) {
        //std::cout << x << " " ;
        std::vector<Cell *>::iterator neighbor;
        std::vector<Cell *>::iterator end = pCell->get_container()->agent_grid[x].end();
        for (neighbor = pCell->get_container()->agent_grid[x].begin(); neighbor != end; ++neighbor) {
            // do not include the neighbor if it is the agent itself or if it is in the list already
            if (std::find(pCell->state.neighbors.begin(), pCell->state.neighbors.end(), (*neighbor)) ==
                pCell->state.neighbors.end() && pCell != (*neighbor)) {
                pCell->state.neighbors.push_back(*neighbor);
            }
        }
    }
    //std::cout << std::endl;
}
// !!! PHYSIMESS CODE BLOCK END !!! //

void fibre_agent_SVG(std::ofstream& os, PhysiCell::Cell* pC, double z_slice, std::vector<std::string> (*cell_coloring_function)(Cell*), double X_lower, double Y_lower) {
	// !!! PHYSIMESS CODE BLOCK START !!! //
	// place a rod if it's a fibre (note fibre already renamed here)
	const auto agentname = std::string(pC->type_name);
	const auto fibre = std::string("fibre");

	if (agentname.find(fibre) != std::string::npos ){

		int crosslinks = pC->custom_data["X_crosslink_count"];
		if (crosslinks >= 3){
			// if fibre has cross-links different colour than if not
			Write_SVG_line(os, (pC->position)[0] - (pC->custom_data["mLength"]) * (pC->state.orientation)[0] - X_lower,
							(pC->position)[1] - (pC->custom_data["mLength"]) * (pC->state.orientation)[1] - Y_lower,
							(pC->position)[0] + (pC->custom_data["mLength"]) * (pC->state.orientation)[0] - X_lower,
							(pC->position)[1] + (pC->custom_data["mLength"]) * (pC->state.orientation)[1] - Y_lower,
							4.0, "darkblue");
		}
		else if (crosslinks == 2){
			// if fibre has cross-links different colour than if not
			Write_SVG_line(os, (pC->position)[0] - (pC->custom_data["mLength"]) * (pC->state.orientation)[0] - X_lower,
							(pC->position)[1] - (pC->custom_data["mLength"]) * (pC->state.orientation)[1] - Y_lower,
							(pC->position)[0] + (pC->custom_data["mLength"]) * (pC->state.orientation)[0] - X_lower,
							(pC->position)[1] + (pC->custom_data["mLength"]) * (pC->state.orientation)[1] - Y_lower,
							4.0, "blue");
		}
		else if (crosslinks == 1){
			// if fibre has cross-links different colour than if not
			Write_SVG_line(os, (pC->position)[0] - (pC->custom_data["mLength"]) * (pC->state.orientation)[0] - X_lower,
							(pC->position)[1] - (pC->custom_data["mLength"]) * (pC->state.orientation)[1] - Y_lower,
							(pC->position)[0] + (pC->custom_data["mLength"]) * (pC->state.orientation)[0] - X_lower,
							(pC->position)[1] + (pC->custom_data["mLength"]) * (pC->state.orientation)[1] - Y_lower,
							4.0, "steelblue");
		}
		else {
			//Write_SVG_fibre( os, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower, plot_radius , 0.5, Colors[1], Colors[0] );
			Write_SVG_line(os, (pC->position)[0] - (pC->custom_data["mLength"]) * (pC->state.orientation)[0] - X_lower,
							(pC->position)[1] - (pC->custom_data["mLength"]) * (pC->state.orientation)[1] - Y_lower,
							(pC->position)[0] + (pC->custom_data["mLength"]) * (pC->state.orientation)[0] - X_lower,
							(pC->position)[1] + (pC->custom_data["mLength"]) * (pC->state.orientation)[1] - Y_lower,
							4.0, "lightskyblue");
		}

	}
	else{
		double r = pC->phenotype.geometry.radius ; 
		double rn = pC->phenotype.geometry.nuclear_radius ; 
		double z = fabs( (pC->position)[2] - z_slice) ; 

		std::vector<std::string> Colors = cell_coloring_function( pC ); 
        
		os << "   <g id=\"cell" << pC->ID << "\" " 
		<< "type=\"" << pC->type_name << "\" "; // new April 2022  
		if( pC->phenotype.death.dead == true )
		{ os << "dead=\"true\" " ; } 
		else
		{ os << "dead=\"false\" " ; } 
		os << ">" << std::endl; 
		
		// figure out how much of the cell intersects with z = 0
		double plot_radius = sqrt( r*r - z*z );

		// then normal cell, plot sphere if it intersects z = 0;
		Write_SVG_circle( os, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower,
							plot_radius , 0.5, Colors[1], Colors[0] );
		// plot the nucleus if it, too intersects z = 0;
		if( fabs(z) < rn && PhysiCell_SVG_options.plot_nuclei == true )
		{
			plot_radius = sqrt( rn*rn - z*z );
			Write_SVG_circle( os, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower,
								plot_radius, 0.5, Colors[3],Colors[2]);
		}
            			os << "   </g>" << std::endl;

	}
}

void fibre_agent_legend(std::ofstream& os, Cell_Definition* cell_definition, double& cursor_x, double& cursor_y, std::vector<std::string> (*cell_coloring_function)(Cell*), double temp_cell_radius, double padding, double font_size) {
	
	// switch to the cell type 
	Cell C; 
	C.convert_to_cell_definition( *cell_definition );
	
	// get the colors using the current coloring function 
	std::vector<std::string> colors = cell_coloring_function(&C); 	

	// place the label 
	// !!! PHYSIMESS CODE BLOCK START !!! //
	// place a rod if it's a fibre (note fibre not yet renamed)
	const auto agentname = std::string(C.type_name);
	const auto ecm = std::string("ecm");
	const auto matrix = std::string("matrix");
	const auto fiber = std::string("fiber");
	const auto fibre = std::string("fibre");
	const auto rod = std::string("rod");

	if (agentname.find(ecm) != std::string::npos ||
		agentname.find(matrix) != std::string::npos ||
		agentname.find(fiber) != std::string::npos ||
		agentname.find(fibre) != std::string::npos ||
		agentname.find(rod) != std::string::npos) {
		//Write_SVG_fibre(os, cursor_x, cursor_y , 0.5*temp_cell_radius , 1.0 , colors[1] , colors[0] );
		Write_SVG_line(os, cursor_x, cursor_y-20.0 , cursor_x , cursor_y+20.0 , 4.0 , "lightskyblue" );
	}
	else {
		// place a big circle with cytoplasm colors 
		Write_SVG_circle(os,cursor_x, cursor_y , temp_cell_radius , 1.0 , colors[1] , colors[0] ); 
		// place a small circle with nuclear colors 
		Write_SVG_circle(os,cursor_x, cursor_y , 0.5*temp_cell_radius , 1.0 , colors[3] , colors[2] ); 
	}
	// !!! PHYSIMESS CODE BLOCK END !!! //
	
	// place the label
	cursor_x += temp_cell_radius + 2*padding; 
	cursor_y += 0.3*font_size; 
	
	Write_SVG_text( os , cell_definition->name.c_str() , cursor_x , cursor_y, font_size , 
		PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	
	// move the cursor down to the next row 
	
	cursor_y -= 0.3*font_size; 
	cursor_y += ( 2.0 * padding + 2.0*temp_cell_radius ); 
	cursor_x = padding + temp_cell_radius;
}