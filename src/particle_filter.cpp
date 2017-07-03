
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

    // Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS)

    //random_device rd;
    //mt19937 gen(rd());
    default_random_engine gen;

    num_particles = 25;
    // Empirical setting: How many particles used in a particle filter. Neither too many nor too few!



    normal_distribution<double> normal_dist_x(x, std[0]);
    normal_distribution<double> normal_dist_y(y, std[1]);
    normal_distribution<double> normal_dist_psi(theta, std[2]);

    particles.resize(num_particles);

    for (int i = 0; i < num_particles; i++) {

        particles[i].id = i;
        particles[i].x = normal_dist_x(gen);  // Add random Gaussian noise to each particle.
        particles[i].y = normal_dist_y(gen);
        particles[i].theta = normal_dist_psi(gen);
        particles[i].weight = 1.0;
        // cout << "Sample " << i + 1 << " " << sample_x << " " << sample_y << " " << sample_psi << endl;
    }

    is_initialized = true;
    //cout << "Initialization complete" << endl;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/

    //Use default random engine
    default_random_engine gen;

    // define normal distributions for sensor noise
    normal_distribution<double> normal_sensor_noise_x(0, std_pos[0]);
    normal_distribution<double> normal_sensor_noise_y(0, std_pos[1]);
    normal_distribution<double> normal_sensor_noise_theta(0, std_pos[2]);

    //Add measurement to each particle:
    for (int i = 0; i < num_particles; i++) {

        // From the bicycle model described as the course material,
        // we can compute the new states, as follows:
        if (fabs(yaw_rate) < 0.00001) { //yaw_rate is zero / almost zero
            particles[i].x += velocity * delta_t * cos(particles[i].theta);
            particles[i].y += velocity * delta_t * sin(particles[i].theta);
        }
        else { //yaw_rate is not zero
            particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
            particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
            particles[i].theta += yaw_rate * delta_t;
        }

        // Add random Gaussian noise:
        particles[i].x += normal_sensor_noise_x(gen);
        particles[i].y += normal_sensor_noise_y(gen);
        particles[i].theta += normal_sensor_noise_theta(gen);
    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.

    for (int i = 0; i < observations.size(); i++) {

        // Get the current observation
        LandmarkObs obs = observations[i];

        // init minimum distance to maximum possible
        double min_dist = numeric_limits<double>::max();

        // init id of landmark from map placeholder to be associated with the observation
        int map_id = -1;

        for (int j = 0; j < predicted.size(); j++) {
            // Get the current prediction
            LandmarkObs pred = predicted[j];

            // get distance between current/predicted landmarks
            double cur_dist = dist(obs.x, obs.y, pred.x, pred.y);

            // find the predicted landmark nearest the current observed landmark
            if (cur_dist < min_dist) {
                min_dist = cur_dist;
                map_id = pred.id;
            }
        }

        // set the observation's id to the nearest predicted landmark's id
        observations[i].id = map_id;
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

    weights.clear();

    for (int i = 0; i < num_particles; i++) {

        double standard_x = std_landmark[0];
        double standard_y = std_landmark[1];
        double standard_x_squared_2 = (standard_x * standard_x);
        double standard_y_squared_2 = (standard_y * standard_y);
        long double part_1 = 1 / (2 * M_PI * standard_x * standard_y);


        long double weight_placeholder = 1.0;

        for (int j = 0; j < observations.size(); j++) {

            double particle_theta = particles[i].theta;
            double cos_theta = cos(particle_theta);
            double sin_theta = sin(particle_theta);

            double transformed_observation_x = particles[i].x + (observations[j].x * cos_theta) -
                                               (observations[j].y * sin_theta);

            double transformed_observation_y = particles[i].y + (observations[j].x * sin_theta) +
                                               (observations[j].y * cos_theta);
            double threshold = 50.0;
            long current_id = -1;

            for (int k = 0; k < map_landmarks.landmark_list.size(); k++) {

                // Find the predicted measurement that is closest to each observed measurement
                double difference = dist(transformed_observation_x, transformed_observation_y,
                                         map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);
                // check if threshold is less than previous
                if (difference < threshold) {
                    // if it is, then we have the current lowest ID
                    // map data is not zero indexed
                    current_id = map_landmarks.landmark_list[k].id_i - 1;
                    // assign the threshold to the current difference so that next loop if the difference is higher it won't update the id.
                    threshold = difference;
                }
            }



            if (current_id >= 0) {

                double x = transformed_observation_x;
                double u_x = map_landmarks.landmark_list[current_id].x_f;
                double y = transformed_observation_y;
                double u_y = map_landmarks.landmark_list[current_id].y_f;
                double x_ux_squared = (x - u_x) * (x - u_x);
                double y_uy_squared = (y - u_y) * (y - u_y);
                double part_2 = x_ux_squared / (standard_x_squared_2)+
                                           y_uy_squared / (standard_y_squared_2);


                weight_placeholder *= part_1 * exp((-1 / 2.) * (part_2));


            }
        }


        particles[i].weight = weight_placeholder;
        weights.push_back(particles[i].weight);


    }

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