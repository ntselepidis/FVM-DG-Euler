#include <ancse/config.hpp>
#include <fstream>

nlohmann::json get_config(const std::string& fileName) {
    /// The path is relative to the binary.
    /// This assumes you're building in a subfolder and
    /// therefore, the config file is located in the parent.
    std::ifstream file(fileName);
    assert(file.good());

    nlohmann::json config;
    file >> config;
    return config;
}

nlohmann::json get_global_config() {
  static nlohmann::json config = get_config("../config.json");
  return config;
}
