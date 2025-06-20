#include <iostream>
// #include <filesystem>
#include <string>
#include <queue>

#ifndef FILECHECK_H
#define FILECHECK_H

// check whether file exist
bool fileExists(const std::string& path);

// function: list all file with a specific suffix.
std::vector<std::string> listfile(std::string indir, std::string suffix);

// function: check if a file is in queue
bool isFileInDeque(const std::deque<std::string>& file_deque, std::string infile);

#endif
