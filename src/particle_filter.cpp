
/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    num_particles=101;
    default_random_engine gen;

    particles.resize(num_particles);
    weights.resize(num_particles);


    // This line creates a normal (Gaussian) distribution for x.
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);


    for (int i = 0; i < num_particles; i++){
        double sample_x, sample_y, sample_theta;
        sample_x = dist_x(gen);
        sample_y = dist_y(gen);
        sample_theta=dist_theta(gen);
        particles[i].id=i;
        particles[i].x=sample_x;
        particles[i].y=sample_y;
        particles[i].theta=sample_theta;
        particles[i].weight=1;
        weights[i]=1;
    }
    is_initialized = true;




}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/

    default_random_engine gen;

    normal_distribution<double> dist_x(0, std_pos[0]);
    normal_distribution<double> dist_y(0, std_pos[0]);
    normal_distribution<double> dist_theta(0, std_pos[0]);


    int N = particles.size();
    for (int i=0; i<N; i++){
        // Using bicycle model to measure each particle
        if (fabs(yaw_rate) < 0.00001){
            particles[i].x += velocity * delta_t * cos(particles[i].theta);
            particles[i].y += velocity * delta_t * sin(particles[i].theta);

        }
        else{
            particles[i].x += (velocity/yaw_rate)*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
            particles[i].y += (velocity/yaw_rate)*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t));
            particles[i].theta += yaw_rate*delta_t;
        }

        // Add random Gaussian noise
        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);
    }
}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.


    // Initially, set the nearest_dist as the max value of double
    double nearest_dist = std::numeric_limits<double>::max();
    double d=0;
    int predicted_id=-1;

    for(int i=0; i<observations.size();i++){
        for(int j=0; j<predicted.size(); j++){
            d=dist(predicted[j].x,predicted[j].y,observations[i].x,observations[i].y);
            if ( d < nearest_dist){
                nearest_dist=d;
                predicted_id=predicted[j].id;
            }
        }
        nearest_dist=std::numeric_limits<double>::max();
        observations[i].id=predicted_id;
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   std::vector<LandmarkObs> observations, Map map_landmarks) {
    // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    //   according to the MAP'S coordinate system. You will need to transform between the two systems.
    //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
    //   The following is a good resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    //   and the following is a good resource for the actual equation to implement (look at equation
    //   3.33
    //   http://planning.cs.uiuc.edu/node99.html


    // The steps are:
    // For every particle, we have to do the following things:
    // (1) Transformations: transforms the observations into MAP's coordinate system.
    //     (http://planning.cs.uiuc.edu/node99.html)
    // (2) Associations: find the Landmark ID that nearest to the transformed observations;
    // (3) Calculating the Particle's Final Weight: the particles final weight will be calculated as
    //     the product of each measurement's Multivariate-Gaussian probability.


    weights.clear();

    for (int i = 0; i < num_particles; i++) {

        double standard_x = std_landmark[0];
        double standard_y = std_landmark[1];
        double standard_x_squared_2 = (standard_x * standard_x);
        double standard_y_squared_2 = (standard_y * standard_y);
        long double part_1 = 1 / (2 * M_PI * standard_x * standard_y);


        long double weight_placeholder = 1.0;


        vector<LandmarkObs> tobs; //transformed observations
        vector<LandmarkObs> predicted; // The map_landmarks
        tobs.resize(observations.size());
        predicted.resize(map_landmarks.landmark_list.size());

        // Transformations
        for(int j=0; j<observations.size(); j++){
            tobs[j].x=observations[j].x*cos(particles[i].theta)-observations[j].y*sin(particles[i].theta)+particles[i].x;
            tobs[j].y=observations[j].x*sin(particles[i].theta)+observations[j].y*cos(particles[i].theta)+particles[i].y;
        }

        // conversion from 'vector<Map::single_landmark_s>' to 'vector<LandmarkObs>'
        for(int j=0; j<map_landmarks.landmark_list.size(); j++){
            predicted[j].id=map_landmarks.landmark_list[j].id_i;
            predicted[j].x=map_landmarks.landmark_list[j].x_f;
            predicted[j].y=map_landmarks.landmark_list[j].y_f;
        }



        // Associations
        for (int j = 0; j < tobs.size(); j++) {
            double transformed_observation_x = tobs[j].x;
            double transformed_observation_y = tobs[j].y;

            double nearest_dist = std::numeric_limits<double>::max();
            long predicted_id = -1;

            for (int k = 0; k < predicted.size(); k++) {
                double d = dist(transformed_observation_x, transformed_observation_y,
                                         predicted[k].x, predicted[k].y);

                if (d < nearest_dist) {
                    predicted_id = predicted[k].id - 1;
                    nearest_dist = d;
                }
            }

            // Calculating the Particle's Final Weight
            double x = transformed_observation_x;
            double u_x = map_landmarks.landmark_list[predicted_id].x_f;
            double y = transformed_observation_y;
            double u_y = map_landmarks.landmark_list[predicted_id].y_f;
            double x_ux_squared = (x - u_x) * (x - u_x);
            double y_uy_squared = (y - u_y) * (y - u_y);
            double part_2 = x_ux_squared / (standard_x_squared_2)+
                            y_uy_squared / (standard_y_squared_2);

            weight_placeholder *= part_1 * exp((-1 / 2.) * (part_2));
        }
        particles[i].weight = weight_placeholder;
        weights.push_back(particles[i].weight);
    }






    /*
    weights.clear();

    for(int i=0; i<num_particles; i++){

        vector<LandmarkObs> tobs; //transformed observations
        vector<LandmarkObs> predicted; // The map_landmarks
        tobs.resize(observations.size());
        predicted.resize(map_landmarks.landmark_list.size());

        // Transformations
        for(int j=0; j<observations.size(); j++){
            tobs[j].x=observations[j].x*cos(particles[i].theta)-observations[j].y*sin(particles[i].theta)+particles[i].x;
            tobs[j].y=observations[j].x*sin(particles[i].theta)+observations[j].y*cos(particles[i].theta)+particles[i].y;
        }

        // conversion from 'vector<Map::single_landmark_s>' to 'vector<LandmarkObs>'
        for(int j=0; j<map_landmarks.landmark_list.size(); j++){
            predicted[j].id=map_landmarks.landmark_list[j].id_i;
            predicted[j].x=map_landmarks.landmark_list[j].x_f;
            predicted[j].y=map_landmarks.landmark_list[j].y_f;
        }

        // Associations
        dataAssociation(predicted,tobs);

        // Calculating the Particle's Final Weight

        long double weight_placeholder = 1.0;
        double sigma_x=std_landmark[0];
        double sigma_y=std_landmark[1];
        cout << endl;
        cout << "sigx=" <<sigma_x <<"sigy=" <<sigma_y <<endl;
        long double temp1=1/(2*M_PI*sigma_x*sigma_y);

        for(int j=0; j<tobs.size(); j++){
            long double x=tobs[j].x;
            long double y=tobs[j].y;
            cout << "x= " << x << "y= " << y <<endl;
            long double mu_x=predicted[tobs[j].id].x;
            long double mu_y=predicted[tobs[j].id].y;
            cout << "mu_x= " << mu_x << "mu_y= " << mu_y <<endl;
            long double temp2 =((x-mu_x)*(x-mu_x))/(2*sigma_x*sigma_x);
            long double temp3 =((y-mu_y)*(y-mu_y))/(2*sigma_y*sigma_y);
            cout << "temp2=" <<temp2 <<endl;
            cout << "temp3=" <<temp3 <<endl;
            long double temp4=exp(-(temp2+temp3));
            cout << "temp4=" <<temp4 <<endl;
            weight_placeholder *= temp1 * temp4;
            cout << "put=" << weight_placeholder << endl;
        }

        particles[i].weight = weight_placeholder;
        weights.push_back(particles[i].weight);
    }

    for (int i=0; i<weights.size();i++){
        cout << weights[i] <<endl;
    }
     */



}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    default_random_engine gen;

    vector<Particle> new_particles;

    // all of the current weights
    vector<double> weights;
    for (int i = 0; i < num_particles; i++) {
        weights.push_back(particles[i].weight);
    }

    // generate a resampling wheel for randomly starting index
    uniform_int_distribution<int> uniintdist(0, num_particles-1);
    auto index = uniintdist(gen);

    // The max weight
    double max_weight = *max_element(weights.begin(), weights.end());

    // uniform random distribution: unirealdist(0.0, max_weight)
    uniform_real_distribution<double> unirealdist(0.0, max_weight);

    double beta = 0.0;

    // choose one from the resample wheel
    for (int i = 0; i < num_particles; i++) {
        beta += unirealdist(gen) * 2.0;
        while (beta > weights[index]) {
            beta -= weights[index];
            index = (index + 1) % num_particles;
        }
        new_particles.push_back(particles[index]);
    }

    particles = new_particles;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    //Clear the previous associations
    particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
    vector<int> v = best.associations;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
    vector<double> v = best.sense_x;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
    vector<double> v = best.sense_y;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}