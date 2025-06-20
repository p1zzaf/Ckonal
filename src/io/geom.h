#ifndef GEOM_H
#define GEOM_H

#include<vector>
#include<string>

class Geometry{
public:
	// ===== atributes =====
	int node_number_;
	std::vector<int> line_code_;// with size of node number
	std::vector<int> node_code_;// with size fo node number
	std::vector<double> node_x_;// with size of node number
	std::vector<double> node_y_;// with size of node number
	std::vector<double> node_elevation_;// with size of node number
	
	// ===== functions =====
	Geometry(const std::string input_geometry_file);
	// ~Geometry();
	int resetElevation(double elevation_new);
	// change x & y
	int exchangexy();
};

#endif


