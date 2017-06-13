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

	num_particles = 100;

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

	normal_distribution<double> dist_x(0, std_pos[0] * delta_t);
	normal_distribution<double> dist_y(0, std_pos[1] * delta_t);
	normal_distribution<double> dist_theta(0, std_pos[2] * delta_t);

	float eps = 0.0001;
	for(Particle &p : particles) {
		if(yaw_rate < eps && yaw_rate > -eps) {
			p.x = p.x + velocity*delta_t*cos(p.theta);
      		p.y = p.y + velocity*delta_t*sin(p.theta);
      		// p.theta = p.theta;
		} else {
			p.x += velocity / yaw_rate * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta)) + dist_x(gen);
			p.y += velocity / yaw_rate * (-cos(p.theta + yaw_rate * delta_t) + cos(p.theta)) + dist_y(gen);
			double theta =  p.theta + yaw_rate * delta_t + dist_theta(gen);
			p.theta = theta;
		}
		// while (p.theta >  M_PI) p.theta -= 2.*M_PI;
		// while (p.theta < -M_PI) p.theta += 2.*M_PI;
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

	for(LandmarkObs& ob : observations) {
		// Find the closest observation from the observation list.
		double bestDist = -1;
		LandmarkObs bestLm;
		for(LandmarkObs lm : predicted) {
			double dist_x = lm.x - ob.x;
			double dist_y = lm.y - ob.y;
			double dist = sqrt(dist_x * dist_x + dist_y * dist_y);
			if(bestDist < 0 || dist < bestDist) {
				bestDist = dist;
				bestLm = lm;
			}
		}
		ob.id = bestLm.id;
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

	double sig_x = std_landmark[0];
    double sig_y = std_landmark[1];
	// 
	// exit(0);

	vector<LandmarkObs> predictedLandmarks;
	for (Map::single_landmark_s candidate_landmark : map_landmarks.landmark_list) {
		predictedLandmarks.push_back(LandmarkObs{
			candidate_landmark.id_i,
			candidate_landmark.x_f,
			candidate_landmark.y_f,
		});
	}


	for(Particle &p : particles) {
		long double prob = 1.0;
		double log_p = 0;
		vector<int> associations;
		vector<double> sense_x;
    	vector<double> sense_y;

		// transform the observations
		std::vector<LandmarkObs> transformed_observations;
		for(LandmarkObs obs : observations) {
			// Transform from car coordinate system to map coordinate system:
			double obs_mod = sqrt(obs.x * obs.x + obs.y * obs.y);
			double alpha = atan2(obs.y, obs.x);

			double map_alpha = p.theta + alpha;

			// This is the landmark coordinates in absolute terms of the map
			double lm_map_x = p.x + obs_mod * cos(map_alpha);
			double lm_map_y = p.y + obs_mod * sin(map_alpha);
			transformed_observations.push_back(LandmarkObs{obs.id, lm_map_x, lm_map_y});
		}
		dataAssociation(predictedLandmarks, transformed_observations);

		
		for(LandmarkObs trans_obs : transformed_observations) {

			// Extract the landmark from the map
			// Iterate through the map and find the matching landmark ID
			// LandmarkObs realLandmarkObs;
			bool found = false;
			Map::single_landmark_s real_landmark;

			for (Map::single_landmark_s candidate_landmark : map_landmarks.landmark_list) {
				if(candidate_landmark.id_i == trans_obs.id) {
					real_landmark = candidate_landmark;
					found = true;
				}
			}

			if(!found) {
				throw "Did not find landmark!";
				continue;
			}
		
			associations.push_back(real_landmark.id_i);
			sense_x.push_back(real_landmark.x_f);
			sense_y.push_back(real_landmark.y_f);

			double diff_x = trans_obs.x - real_landmark.x_f;
			double diff_y = trans_obs.y - real_landmark.y_f;

			long double denomin = 1.0/(2.0 * M_PI * sig_x * sig_y);
			long double x_nom = (diff_x * diff_x) / (sig_x * sig_x);
			long double y_nom = (diff_y * diff_y) / (sig_y * sig_y);
			long double expon = exp(-0.5*(x_nom+y_nom));
			long double meas_likelihood = denomin * expon;

      		prob *= meas_likelihood;
		}
		
		
		p.weight = prob;
		// p = SetAssociations(p, associations, sense_x, sense_y);
		
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
