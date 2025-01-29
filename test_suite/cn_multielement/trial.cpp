#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>

using namespace std;

// Function to parse a comma-separated list into a vector
vector<string> parseStringList(const string& str) {
    vector<string> result;
    stringstream ss(str);
    string token;
    while (getline(ss, token, ',')) {
        result.push_back(token);
    }
    return result;
}

vector<double> parseDoubleList(const string& str) {
    vector<double> result;
    stringstream ss(str);
    string token;
    while (getline(ss, token, ',')) {
        result.push_back(stod(token));
    }
    return result;
}

// Function to read config file into a key-value map
unordered_map<string, string> readConfig(const string& filename) {
    unordered_map<string, string> config;
    ifstream file(filename);
    
    if (!file) {
        cerr << "Error: Cannot open config file " << filename << endl;
        exit(1);
    }

    string line;
    while (getline(file, line)) {
        // Ignore empty lines or comments
        if (line.empty() || line[0] == '#') continue;

        // Find the '=' character
        size_t pos = line.find('=');
        if (pos == string::npos) continue; // Skip invalid lines

        // Extract key and value
        string key = line.substr(0, pos);
        string value = line.substr(pos + 1);

        // Trim spaces
        key.erase(key.find_last_not_of(" \t") + 1);
        key.erase(0, key.find_first_not_of(" \t"));

        value.erase(value.find_last_not_of(" \t") + 1);
        value.erase(0, value.find_first_not_of(" \t"));

        // Store in the map
        config[key] = value;
    }
    return config;
}

int main(int argc, char* argv[]) {
    // Check if the filename argument is provided
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <config_file>\n";
        return 1;
    }

    string configFile = argv[1];

    // Read the config file
    unordered_map<string, string> config = readConfig(configFile);

    // Retrieve and parse values
    double rcout_2b = stod(config["RCOUT_2B"]);
    double rcout_3b = stod(config["RCOUT_3B"]);
    double rcout_4b = stod(config["RCOUT_4B"]);
    
    vector<string> atom_pairs = parseStringList(config["ATOM_PAIRS"]);
    vector<double> lambda_set = parseDoubleList(config["LAMBDA_SET"]);
    vector<double> rcin_list = parseDoubleList(config["RCIN_LIST"]);

    // Output parsed values
    cout << "Using config file: " << configFile << endl;
    cout << "rcout_2b: " << rcout_2b << endl;
    cout << "rcout_3b: " << rcout_3b << endl;
    cout << "rcout_4b: " << rcout_4b << endl;

    cout << "Atom Pairs: ";
    for (const auto& pair : atom_pairs) cout << pair << " ";
    cout << endl;

    cout << "Lambda Set: ";
    for (const auto& lambda : lambda_set) cout << lambda << " ";
    cout << endl;

    cout << "RCIN List: ";
    for (const auto& rcin : rcin_list) cout << rcin << " ";
    cout << endl;

    return 0;
}
