#include "pauli_bench/reference_registry.hpp"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#ifndef PAULI_REFERENCE_REGISTRY_DEFAULT_PATH
#define PAULI_REFERENCE_REGISTRY_DEFAULT_PATH "references/reference_registry.json"
#endif

namespace pauli_bench {
namespace {

struct JsonValue {
    enum class Kind { Null, Boolean, Number, String, Array, Object };

    Kind kind{Kind::Null};
    std::string scalar;
    std::vector<JsonValue> array;
    std::map<std::string, JsonValue> object;
};

class JsonParser {
  public:
    explicit JsonParser(std::string input) : input_(std::move(input)) {}

    JsonValue parse() {
        JsonValue result = parse_value();
        skip_whitespace();
        if (position_ != input_.size()) {
            fail("unexpected trailing content");
        }
        return result;
    }

  private:
    [[noreturn]] void fail(const std::string& message) const {
        throw std::runtime_error(
            "invalid reference registry JSON at byte " +
            std::to_string(position_) + ": " + message);
    }

    void skip_whitespace() {
        while (position_ < input_.size()) {
            const char c = input_[position_];
            if (c != ' ' && c != '\n' && c != '\r' && c != '\t') {
                break;
            }
            ++position_;
        }
    }

    bool consume(char expected) {
        skip_whitespace();
        if (position_ < input_.size() && input_[position_] == expected) {
            ++position_;
            return true;
        }
        return false;
    }

    void expect(char expected) {
        if (!consume(expected)) {
            fail(std::string("expected '") + expected + "'");
        }
    }

    JsonValue parse_value() {
        skip_whitespace();
        if (position_ == input_.size()) {
            fail("expected a value");
        }
        switch (input_[position_]) {
            case 'n':
                return parse_literal("null", JsonValue::Kind::Null);
            case 't':
                return parse_literal("true", JsonValue::Kind::Boolean);
            case 'f':
                return parse_literal("false", JsonValue::Kind::Boolean);
            case '"': {
                JsonValue value;
                value.kind = JsonValue::Kind::String;
                value.scalar = parse_string();
                return value;
            }
            case '[':
                return parse_array();
            case '{':
                return parse_object();
            default:
                if (input_[position_] == '-' ||
                    (input_[position_] >= '0' && input_[position_] <= '9')) {
                    return parse_number();
                }
                fail("unexpected character");
        }
    }

    JsonValue parse_literal(std::string_view literal, JsonValue::Kind kind) {
        if (input_.compare(position_, literal.size(), literal) != 0) {
            fail("invalid literal");
        }
        position_ += literal.size();
        JsonValue value;
        value.kind = kind;
        value.scalar = std::string(literal);
        return value;
    }

    static int hex_digit(char c) {
        if (c >= '0' && c <= '9') {
            return c - '0';
        }
        if (c >= 'a' && c <= 'f') {
            return 10 + c - 'a';
        }
        if (c >= 'A' && c <= 'F') {
            return 10 + c - 'A';
        }
        return -1;
    }

    unsigned parse_unicode_code_unit() {
        if (position_ + 4 > input_.size()) {
            fail("incomplete unicode escape");
        }
        unsigned codepoint = 0;
        for (int index = 0; index < 4; ++index) {
            const int digit = hex_digit(input_[position_++]);
            if (digit < 0) {
                fail("invalid unicode escape");
            }
            codepoint = (codepoint << 4U) | static_cast<unsigned>(digit);
        }
        return codepoint;
    }

    void append_unicode(std::string& output) {
        unsigned codepoint = parse_unicode_code_unit();
        if (codepoint >= 0xd800 && codepoint <= 0xdbff) {
            if (position_ + 2 > input_.size() || input_[position_] != '\\' ||
                input_[position_ + 1] != 'u') {
                fail("high unicode surrogate has no low surrogate");
            }
            position_ += 2;
            const unsigned low = parse_unicode_code_unit();
            if (low < 0xdc00 || low > 0xdfff) {
                fail("invalid low unicode surrogate");
            }
            codepoint = 0x10000U + ((codepoint - 0xd800U) << 10U) +
                        (low - 0xdc00U);
        } else if (codepoint >= 0xdc00 && codepoint <= 0xdfff) {
            fail("unexpected low unicode surrogate");
        }
        if (codepoint <= 0x7f) {
            output.push_back(static_cast<char>(codepoint));
        } else if (codepoint <= 0x7ff) {
            output.push_back(static_cast<char>(0xc0 | (codepoint >> 6U)));
            output.push_back(static_cast<char>(0x80 | (codepoint & 0x3fU)));
        } else if (codepoint <= 0xffff) {
            output.push_back(static_cast<char>(0xe0 | (codepoint >> 12U)));
            output.push_back(static_cast<char>(0x80 | ((codepoint >> 6U) & 0x3fU)));
            output.push_back(static_cast<char>(0x80 | (codepoint & 0x3fU)));
        } else {
            output.push_back(static_cast<char>(0xf0 | (codepoint >> 18U)));
            output.push_back(static_cast<char>(0x80 | ((codepoint >> 12U) & 0x3fU)));
            output.push_back(static_cast<char>(0x80 | ((codepoint >> 6U) & 0x3fU)));
            output.push_back(static_cast<char>(0x80 | (codepoint & 0x3fU)));
        }
    }

    std::string parse_string() {
        expect('"');
        std::string result;
        while (position_ < input_.size()) {
            const char c = input_[position_++];
            if (c == '"') {
                return result;
            }
            if (static_cast<unsigned char>(c) < 0x20) {
                fail("unescaped control character in string");
            }
            if (c != '\\') {
                result.push_back(c);
                continue;
            }
            if (position_ == input_.size()) {
                fail("incomplete string escape");
            }
            const char escaped = input_[position_++];
            switch (escaped) {
                case '"': result.push_back('"'); break;
                case '\\': result.push_back('\\'); break;
                case '/': result.push_back('/'); break;
                case 'b': result.push_back('\b'); break;
                case 'f': result.push_back('\f'); break;
                case 'n': result.push_back('\n'); break;
                case 'r': result.push_back('\r'); break;
                case 't': result.push_back('\t'); break;
                case 'u': append_unicode(result); break;
                default: fail("invalid string escape");
            }
        }
        fail("unterminated string");
    }

    JsonValue parse_number() {
        const std::size_t start = position_;
        if (input_[position_] == '-') {
            ++position_;
        }
        if (position_ == input_.size()) {
            fail("incomplete number");
        }
        if (input_[position_] == '0') {
            ++position_;
        } else if (input_[position_] >= '1' && input_[position_] <= '9') {
            while (position_ < input_.size() && input_[position_] >= '0' &&
                   input_[position_] <= '9') {
                ++position_;
            }
        } else {
            fail("invalid number");
        }
        if (position_ < input_.size() && input_[position_] == '.') {
            ++position_;
            const std::size_t digits = position_;
            while (position_ < input_.size() && input_[position_] >= '0' &&
                   input_[position_] <= '9') {
                ++position_;
            }
            if (digits == position_) {
                fail("fraction has no digits");
            }
        }
        if (position_ < input_.size() &&
            (input_[position_] == 'e' || input_[position_] == 'E')) {
            ++position_;
            if (position_ < input_.size() &&
                (input_[position_] == '+' || input_[position_] == '-')) {
                ++position_;
            }
            const std::size_t digits = position_;
            while (position_ < input_.size() && input_[position_] >= '0' &&
                   input_[position_] <= '9') {
                ++position_;
            }
            if (digits == position_) {
                fail("exponent has no digits");
            }
        }
        JsonValue value;
        value.kind = JsonValue::Kind::Number;
        value.scalar = input_.substr(start, position_ - start);
        return value;
    }

    JsonValue parse_array() {
        expect('[');
        JsonValue value;
        value.kind = JsonValue::Kind::Array;
        if (consume(']')) {
            return value;
        }
        while (true) {
            value.array.push_back(parse_value());
            if (consume(']')) {
                return value;
            }
            expect(',');
        }
    }

    JsonValue parse_object() {
        expect('{');
        JsonValue value;
        value.kind = JsonValue::Kind::Object;
        if (consume('}')) {
            return value;
        }
        while (true) {
            skip_whitespace();
            if (position_ == input_.size() || input_[position_] != '"') {
                fail("object key must be a string");
            }
            std::string key = parse_string();
            expect(':');
            const auto [iterator, inserted] =
                value.object.emplace(std::move(key), parse_value());
            if (!inserted) {
                fail("duplicate object key '" + iterator->first + "'");
            }
            if (consume('}')) {
                return value;
            }
            expect(',');
        }
    }

    std::string input_;
    std::size_t position_{};
};

const JsonValue& field(const JsonValue& object, std::string_view name) {
    if (object.kind != JsonValue::Kind::Object) {
        throw std::runtime_error("reference registry value must be an object");
    }
    const auto iterator = object.object.find(std::string(name));
    if (iterator == object.object.end()) {
        throw std::runtime_error(
            "reference registry field is missing: " + std::string(name));
    }
    return iterator->second;
}

const JsonValue* optional_field(
    const JsonValue& object, std::string_view name) {
    if (object.kind != JsonValue::Kind::Object) {
        throw std::runtime_error("reference registry value must be an object");
    }
    const auto iterator = object.object.find(std::string(name));
    return iterator == object.object.end() ? nullptr : &iterator->second;
}

std::string string_value(const JsonValue& value, std::string_view name) {
    if (value.kind != JsonValue::Kind::String) {
        throw std::runtime_error(
            "reference registry field must be a string: " + std::string(name));
    }
    return value.scalar;
}

double double_value(const JsonValue& value, std::string_view name) {
    if (value.kind != JsonValue::Kind::Number) {
        throw std::runtime_error(
            "reference registry field must be numeric: " + std::string(name));
    }
    errno = 0;
    char* end = nullptr;
    const double result = std::strtod(value.scalar.c_str(), &end);
    if (errno == ERANGE || end != value.scalar.c_str() + value.scalar.size() ||
        !std::isfinite(result)) {
        throw std::runtime_error(
            "invalid numeric reference registry field: " + std::string(name));
    }
    return result;
}

std::uint64_t uint_value(const JsonValue& value, std::string_view name) {
    if (value.kind != JsonValue::Kind::Number || value.scalar.empty() ||
        value.scalar.front() == '-' || value.scalar.find_first_of(".eE") != std::string::npos) {
        throw std::runtime_error(
            "reference registry field must be a nonnegative integer: " +
            std::string(name));
    }
    std::size_t consumed = 0;
    try {
        const auto result = std::stoull(value.scalar, &consumed);
        if (consumed != value.scalar.size()) {
            throw std::invalid_argument("trailing content");
        }
        return result;
    } catch (const std::exception&) {
        throw std::runtime_error(
            "invalid integer reference registry field: " + std::string(name));
    }
}

int int_value(const JsonValue& value, std::string_view name) {
    const std::uint64_t result = uint_value(value, name);
    if (result > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error(
            "reference registry integer is too large: " + std::string(name));
    }
    return static_cast<int>(result);
}

std::optional<double> optional_double(
    const JsonValue& value, std::string_view name) {
    if (value.kind == JsonValue::Kind::Null) {
        return std::nullopt;
    }
    return double_value(value, name);
}

std::optional<std::uint64_t> optional_uint(
    const JsonValue& value, std::string_view name) {
    if (value.kind == JsonValue::Kind::Null) {
        return std::nullopt;
    }
    return uint_value(value, name);
}

std::uint64_t mask_value(const JsonValue& value, std::string_view name) {
    const std::string text = string_value(value, name);
    std::size_t consumed = 0;
    try {
        const auto result = std::stoull(text, &consumed, 0);
        if (consumed != text.size()) {
            throw std::invalid_argument("trailing content");
        }
        return result;
    } catch (const std::exception&) {
        throw std::runtime_error(
            "invalid Pauli mask in reference registry: " + std::string(name));
    }
}

struct RegistryEntry {
    ReferenceQuery query;
    StoredReference reference;
};

RegistryEntry parse_entry(const JsonValue& value) {
    RegistryEntry entry;
    entry.query.model = string_value(field(value, "model"), "model");
    entry.query.qubits = int_value(field(value, "qubits"), "qubits");
    entry.query.layers = int_value(field(value, "layers"), "layers");
    entry.query.circuit_generation_version = int_value(
        field(value, "circuit_generation_version"), "circuit_generation_version");
    entry.query.circuit_seed = optional_uint(field(value, "circuit_seed"), "circuit_seed");

    const JsonValue& parameters = field(value, "parameters");
    entry.query.t_density = optional_double(field(parameters, "t_density"), "t_density");
    entry.query.depolarizing_probability = optional_double(
        field(parameters, "depolarizing_probability"), "depolarizing_probability");
    entry.query.dt = optional_double(field(parameters, "dt"), "dt");
    entry.query.coupling = optional_double(field(parameters, "coupling"), "coupling");
    entry.query.transverse_field = optional_double(
        field(parameters, "transverse_field"), "transverse_field");
    entry.query.longitudinal_field = optional_double(
        field(parameters, "longitudinal_field"), "longitudinal_field");
    const JsonValue* prefix_depth = optional_field(parameters, "prefix_depth");
    entry.query.prefix_depth = prefix_depth == nullptr
                                   ? std::nullopt
                                   : optional_double(*prefix_depth, "prefix_depth");

    const JsonValue& observable = field(value, "observable");
    entry.query.observable.x = mask_value(field(observable, "x_mask"), "x_mask");
    entry.query.observable.z = mask_value(field(observable, "z_mask"), "z_mask");

    entry.reference.id = string_value(field(value, "id"), "id");
    const JsonValue& reference = field(value, "reference");
    entry.reference.value = double_value(field(reference, "value"), "reference.value");
    entry.reference.method = string_value(field(reference, "method"), "reference.method");
    entry.reference.precision = string_value(
        field(reference, "precision"), "reference.precision");
    entry.reference.uncertainty = optional_double(
        field(reference, "uncertainty"), "reference.uncertainty");
    return entry;
}

std::vector<RegistryEntry> load_registry(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open reference registry: " + path);
    }
    std::ostringstream buffer;
    buffer << input.rdbuf();
    const JsonValue document = JsonParser(buffer.str()).parse();
    const auto schema_version = uint_value(
        field(document, "schema_version"), "schema_version");
    if (schema_version != 1 && schema_version != 2) {
        throw std::runtime_error("unsupported reference registry schema_version");
    }
    const JsonValue& references = field(document, "references");
    if (references.kind != JsonValue::Kind::Array) {
        throw std::runtime_error("reference registry references must be an array");
    }
    std::vector<RegistryEntry> entries;
    entries.reserve(references.array.size());
    for (const JsonValue& value : references.array) {
        entries.push_back(parse_entry(value));
    }
    return entries;
}

bool same_query(const ReferenceQuery& left, const ReferenceQuery& right) {
    return left.model == right.model && left.qubits == right.qubits &&
           left.layers == right.layers &&
           left.circuit_generation_version == right.circuit_generation_version &&
           left.circuit_seed == right.circuit_seed &&
           left.t_density == right.t_density &&
           left.depolarizing_probability == right.depolarizing_probability &&
           left.dt == right.dt && left.coupling == right.coupling &&
           left.transverse_field == right.transverse_field &&
           left.longitudinal_field == right.longitudinal_field &&
           left.prefix_depth == right.prefix_depth &&
           left.observable == right.observable;
}

}  // namespace

std::string reference_registry_path() {
    const char* override_path = std::getenv("PAULI_REFERENCE_REGISTRY");
    if (override_path != nullptr && override_path[0] != '\0') {
        return override_path;
    }
    return PAULI_REFERENCE_REGISTRY_DEFAULT_PATH;
}

std::optional<StoredReference> find_stored_reference(
    const ReferenceQuery& query) {
    std::optional<StoredReference> match;
    for (const RegistryEntry& entry : load_registry(reference_registry_path())) {
        if (!same_query(query, entry.query)) {
            continue;
        }
        if (match.has_value()) {
            throw std::runtime_error(
                "ambiguous reference registry configuration: " + match->id +
                " and " + entry.reference.id);
        }
        match = entry.reference;
    }
    return match;
}

void attach_stored_reference(Circuit& circuit, const ReferenceQuery& query) {
    const auto reference = find_stored_reference(query);
    if (!reference.has_value()) {
        circuit.reference = std::numeric_limits<double>::quiet_NaN();
        circuit.reference_method = "unavailable";
        circuit.reference_id.clear();
        return;
    }
    circuit.reference = reference->value;
    circuit.reference_method = reference->method;
    circuit.reference_id = reference->id;
}

}  // namespace pauli_bench
