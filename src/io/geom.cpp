#include<vector>
#include<string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "geom.h"

Geometry::Geometry(const std::string input_geometry_file){
/*
- input geometry file format: line_id node_id x y elevation(up:+)
*/
	std::string line_string_tmp;
	int icount = 0;
	int line_code, node_code;
	double node_x, node_y, node_z;
	
	std::cout << "**Geometry::Geometry** input file: " << input_geometry_file << std::endl;

	std::ifstream fp(input_geometry_file);
	if(!fp.is_open()){
		std::cerr << "Failed to open file: " << input_geometry_file << std::endl;
	}
	line_code_.clear();
	node_code_.clear();
	node_x_.clear();
	node_y_.clear();
	node_elevation_.clear();
	// read geometry per line
	while(std::getline(fp, line_string_tmp)){
		if(line_string_tmp[0] == '#') continue;
		
		std::istringstream iss_tmp(line_string_tmp);
		if(!(iss_tmp >> line_code >> node_code >> node_x >> node_y >> node_z)) break;// 如果读取失败，则退出循环
		line_code_.push_back(line_code);
		node_code_.push_back(node_code);
		
		node_x_.push_back(node_x);
		node_y_.push_back(node_y);
		node_elevation_.push_back(node_z);
		
		icount ++;
	}
	
	node_number_ = icount;
	std::cout << "**Geometry::Geometry** read [" << icount << "] traces from geomtry file." << std::endl;
	fp.close();
}


int Geometry::resetElevation(double elevation_new){
	int itrace = 0;
	int num_trace;
	num_trace = node_number_;
	// std::cout << "**Geometry::resetElevation**" << num_trace << " traces reset elevation " << std::endl;
	for(itrace=0;itrace<num_trace;itrace++){
		// std::cout << itrace << std::endl;
		node_elevation_[itrace] = elevation_new;
	}
	return(0);
}

// change x & y
int Geometry::exchangexy(){
	double var1;
	for(int i=0;i<node_number_;i++){
		var1 = node_x_[i];
		node_x_[i] = node_y_[i];
		node_y_[i] = var1;
	}
	// return
	return(0);
}

