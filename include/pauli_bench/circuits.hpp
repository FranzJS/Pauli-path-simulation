#pragma once
#include "types.hpp"
namespace pb {
Circuit make_tensor_rx_circuit(int n, int active_qubits, double theta);
std::vector<Circuit> width_suite();
std::vector<Circuit> branching_suite();
std::vector<Circuit> variance_suite();
}
