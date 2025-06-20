// --- START OF FILE filecheck.cpp (已修复) ---

#include "filecheck.h"
#include <filesystem> // 使用 C++17 标准库
#include <deque>
#include <algorithm>
#include <iostream>

// C++17 标准方法检查文件是否存在
bool fileExists(const std::string& path) {
	if (path.empty()) {
		return false;
	}
	// std::filesystem::exists 可以正确处理文件和目录
	return std::filesystem::exists(path);
}

// 使用 C++17 标准库列出目录下指定后缀的文件
std::vector<std::string> listfile(std::string indir, std::string suffix) {
	std::vector<std::string> files_list;

	// 检查目录是否存在
	if (!std::filesystem::exists(indir) || !std::filesystem::is_directory(indir)) {
		return files_list; // 如果目录不存在，返回空列表
	}

	// 遍历目录
	for (const auto& entry : std::filesystem::directory_iterator(indir)) {
		// 检查是否是常规文件并且后缀名匹配
		if (entry.is_regular_file() && entry.path().extension() == suffix) {
			files_list.push_back(entry.path().string()); // 返回完整路径
		}
	}
	return files_list;
}

bool isFileInDeque(const std::deque<std::string>& file_deque, std::string infile) {
	auto it = std::find(file_deque.begin(), file_deque.end(), infile);
	if (it != file_deque.end()) {
		return true;
	}
	else {
		return false;
	}
}