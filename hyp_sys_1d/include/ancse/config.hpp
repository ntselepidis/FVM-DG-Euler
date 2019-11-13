#ifndef HYPSYS1D_CONFIG_HPP
#define HYPSYS1D_CONFIG_HPP

#include <nlohmann/json.hpp>

/// Returns the config stored at `config.json`.
// Note: this is easily misused as essentially a global variable.
nlohmann::json get_global_config();

nlohmann::json get_config(const std::string& fileName);

#endif // HYPSYS1D_CONFIG_HPP
