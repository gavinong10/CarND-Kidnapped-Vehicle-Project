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

	// Let us start with num_particles first set to 100.
	default_random_engine gen;

	int num_particles = 100;

	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	float default_weight = 1. / num_particles;

	for (int i = 0; i < num_particles; i++) {
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		double theta = dist_theta(gen);

		while (theta >  M_PI) theta -= 2.*M_PI;
    	while (theta < -M_PI) theta += 2.*M_PI;
		p.theta = theta;
		// p.weight = default_weight;

		particles.push_back(p);
		weights.push_back(default_weight);
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
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);

	for(Particle p : particles) {
		p.x += velocity / yaw_rate * (sin(p.theta + yaw_rate) - sin(p.theta)) + dist_x(gen);
		p.y += velocity / yaw_rate * (-cos(p.theta + yaw_rate) + cos(p.theta)) + dist_y(gen);
		p.theta += yaw_rate * delta_t + dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// struct LandmarkObs {
	// 	int id;				// Id of matching landmark in the map.
	// 	double x;			// Local (vehicle coordinates) x position of landmark observation [m]
	// 	double y;			// Local (vehicle coordinates) y position of landmark observation [m]
	// };

	// ???

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

	double sig_x = std_landmark[0];
    double sig_y = std_landmark[1];

	for(Particle &p : particles) {
		double prob;
		double log_p = 0;
		vector<int> associations;
		vector<double> sense_x;
    	vector<double> sense_y;

		for(LandmarkObs obs : observations) {
			// Transform from car coordinate system to map coordinate system:
			double obs_mod = sqrt(obs.x * obs.x + obs.y * obs.y);
			double alpha = atan2(obs.y, obs.x);

			double map_alpha = p.theta + alpha;

			// This is the landmark coordinates in absolute terms of the map
			double lm_map_x = p.x + obs_mod * cos(map_alpha);
			double lm_map_y = p.y + obs_mod * sin(map_alpha);


			// Extract the landmark from the map
			// Iterate through the map and find the matching landmark ID
			// LandmarkObs realLandmarkObs;
			bool found = false;
			Map::single_landmark_s real_landmark;
			for (Map::single_landmark_s candidate_landmark : map_landmarks.landmark_list) {
				if(candidate_landmark.id_i == obs.id) {
					// realLandmarkObs.id = real_landmark.id_i;
					// realLandmarkObs.x = real_landmark.x_f;
					// realLandmarkObs.y = real_landmark.x_f;
					real_landmark = candidate_landmark;
					found = true;
				}
			}

			if(!found) {
				// log_p = 0;
				// break;
				throw "Did not find landmark!";
			}

			associations.push_back(obs.id);
			sense_x.push_back(obs.x);
			sense_y.push_back(obs.y);

			double diff_x = lm_map_x - real_landmark.x_f;
			double diff_y = lm_map_y - real_landmark.y_f;

			log_p += -1 / 2 * (diff_x * diff_x / (sig_x * sig_x) + diff_y * diff_y / (sig_y * sig_y)) - log(2.*M_PI*sig_x*sig_y);

		}
		prob = exp(log_p);
		p.weight = prob;
		p = SetAssociations(p, associations, sense_x, sense_y);
		
	}
	for (int i = 0; i < num_particles; ++i) {
		weights[i] = particles[i].weight;
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> dist(weights.begin(), weights.end());

	std::vector<Particle> new_particles;
    for(int n=0; n<num_particles; ++n) {
        int idx = dist(gen);
		new_particles.push_back(particles[idx]);
		weights[n] = particles[idx].weight;
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
