/*
This script takes in two .xyz(f) files and:
+ Produces distances lists for every 2, 3, and 4-body cluster
+ ...
Expects NON_ORTHO format, but assumes an orthorhombic box
assumes all atoms are the same type
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <sys/time.h>
#include <thread>
#include <omp.h>  // For OpenMP

using namespace std;

// Declare ALPHA as a global variable
double alpha = 0.1;
// Global weights map
map<string, float> weights;

struct xyz {
    double x;
    double y;
    double z;
    string atom_type; // Added for multi-element
};

// Helper: returns the two-character key (sorted) for an atom pair
inline string makeKey(const string &a, const string &b) {
    return (a < b) ? a + b : b + a;
}

// Function to parse a comma-separated list into a vector of strings
vector<string> parseStringList(const string& str) {
    vector<string> result;
    stringstream ss(str);
    string token;
    while(getline(ss, token, ',')) {
        result.push_back(token);
    }
    return result;
}

vector<double> parseDoubleList(const string& str) {
    vector<double> result;
    stringstream ss(str);
    string token;
    while(getline(ss, token, ',')) {
        result.push_back(stod(token));
    }
    return result;
}

// Function to generate unique atom pairs from atom types, sorted alphabetically
vector<string> generateAtomPairs(const vector<string>& atom_types) {
    vector<string> atom_pairs;
    // Generate all unique combinations of atom pairs
    for (size_t i = 0; i < atom_types.size(); ++i) {
        for (size_t j = i; j < atom_types.size(); ++j) {  // Ensure i <= j for uniqueness
            string pair = makeKey(atom_types[i], atom_types[j]);
            atom_pairs.push_back(pair);
        }
    }
    return atom_pairs;
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
    while(getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        size_t pos = line.find('=');
        if (pos == string::npos) continue;
        string key = line.substr(0, pos);
        string value = line.substr(pos + 1);
        // Trim spaces:
        key.erase(key.find_last_not_of(" \t") + 1);
        key.erase(0, key.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t") + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        config[key] = value;
    }
    return config;
}

int split_line(string line, vector<string> &items) {
    // Break a line up into tokens based on space separators.
    // Returns the number of tokens parsed.
    string contents;
    stringstream sstream;
    // Strip comments beginning with ! or ## and terminal new line
    int pos = line.find('!');
    if ( pos != string::npos )
        line.erase(pos, line.length() - pos);
    pos = line.find("##");
    if ( pos != string::npos )
        line.erase(pos, line.length()-pos);
    pos = line.find('\n');
    if ( pos != string::npos )
        line.erase(pos, 1);
    sstream.str(line);
    items.clear();
    while (sstream >> contents)
        items.push_back(contents);
    return items.size();
}

bool get_next_line(istream& str, string & line) {
    getline(str, line);
    return bool(str);
}

// Inline helper: Compute periodic displacement along one coordinate.
inline double periodic_diff(double d, double box_length) {
    return d - box_length * round(d / box_length);
}

// Inlined version of get_dist. Note: returns unweighted distance.
inline double get_dist(const xyz &box, const xyz &a1, const xyz &a2,
                       const map<string, pair<float, double>> &atom_map) {
    string key = makeKey(a1.atom_type, a2.atom_type);
    double rcin = 0.0;
    auto it = atom_map.find(key);
    if (it != atom_map.end()) {
        rcin = it->second.second;
    } else {
        // In production code, you might want to handle this error differently.
        return -1;
    }
    double dx = periodic_diff(a1.x - a2.x, box.x);
    double dy = periodic_diff(a1.y - a2.y, box.y);
    double dz = periodic_diff(a1.z - a2.z, box.z);
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    if (dist < rcin) {
        // Something went wrong; flag an error.
        dist = 1e8;
    }
    return dist;
}

// Inlined version of transform.
inline double transform(const xyz &a1, const xyz &a2,
                        const map<string, pair<float, double>> &atom_map,
                        double rcout, double rij) {
    string key = makeKey(a1.atom_type, a2.atom_type);
    float lambda = 0.0f;
    double rcin = 0.0;
    auto it = atom_map.find(key);
    if (it != atom_map.end()) {
        lambda = it->second.first;
        rcin = it->second.second;
    } else {
        return -1;
    }
    double x_min = exp(-rcin / lambda);
    double x_max = exp(-rcout / lambda);
    double x_avg  = 0.5 * (x_max + x_min);
    double x_diff = -0.5 * (x_max - x_min); // Special for Morse style
    return (exp(-rij / lambda) - x_avg) / x_diff;
}

void writeClusterData(const string& filename, const vector<vector<double>>& data) {
    ofstream ofs(filename);
    if (!ofs) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }
    ostringstream buffer;
    for (const auto& row : data) {
        for (size_t j = 0; j < row.size(); j++) {
            buffer << row[j];
            if (j < row.size() - 1)
                buffer << " ";
        }
        buffer << "\n";
    }
    ofs << buffer.str();
    ofs.close();
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <config_file>\n";
        return 1;
    }
    string configFile = argv[1];
    unordered_map<string, string> config = readConfig(configFile);

    double rcout_2b = stod(config["RCOUT_2B"]);
    double rcout_3b = stod(config["RCOUT_3B"]);
    double rcout_4b = stod(config["RCOUT_4B"]);

    vector<string> atom_types = parseStringList(config["ATOM_TYPES"]);
    vector<string> atom_pairs = generateAtomPairs(atom_types);
    vector<double> rcin_list = parseDoubleList(config["RCIN_LIST"]);
    vector<double> lambda_set;
    vector<double> property_set;
    if (config.find("LAMBDA_SET") != config.end())
        lambda_set = parseDoubleList(config["LAMBDA_SET"]);
    if (config.find("PROPERTY_SET") != config.end())
        property_set = parseDoubleList(config["PROPERTY_SET"]);

    if (property_set.size() != atom_types.size()) {
        cerr << "Error: PROPERTY_SET size does not match ATOM_TYPES size.\n";
        return 1;
    }
    unordered_map<string, double> property_map;
    for (size_t i = 0; i < atom_types.size(); ++i)
        property_map[atom_types[i]] = property_set[i];

    cout << "Using config file: " << configFile << "\n";
    cout << "rcout_2b: " << rcout_2b << "\n";
    cout << "rcout_3b: " << rcout_3b << "\n";
    cout << "rcout_4b: " << rcout_4b << "\n";
    cout << "Atom Pairs: ";
    for (const auto& pair : atom_pairs)
        cout << pair << " ";
    cout << "\n";

    if (!lambda_set.empty()) {
        cout << "Lambda Set: ";
        for (const auto& lambda : lambda_set)
            cout << lambda << " ";
        cout << "\n";
    } else {
        cout << "Lambda set not provided\n";
    }

    cout << "RCIN List: ";
    for (const auto& rcin : rcin_list)
        cout << rcin << " ";
    cout << "\n";

    if (config.find("ALPHA") != config.end()) {
        alpha = stod(config["ALPHA"]);
    } else {
        cerr << "Warning: ALPHA not found in config file. Using default value: " << alpha << "\n";
    }
    
    // Build the atom_map: key -> (lambda, rcin)
    map<string, pair<float, double>> atom_map;
    // Assume the ordering of atom_pairs, rcin_list, and property combination is consistent.
    vector<double> prop_combined;
    for (const auto &pair : atom_pairs) {
        // Assume pair length == 2
        string atom1(1, pair[0]);
        string atom2(1, pair[1]);
        if (property_map.find(atom1) == property_map.end() || property_map.find(atom2) == property_map.end()) {
            cerr << "Error: Atom type not found in PROPERTY_SET: " << pair << "\n";
            return 1;
        }
        double prop_value = sqrt(property_map[atom1] * property_map[atom2]);
        prop_combined.push_back(prop_value);
    }
    for (size_t i = 0; i < atom_pairs.size(); ++i)
        atom_map[atom_pairs[i]] = make_pair(prop_combined[i], rcin_list[i]);

    double min_weight_val = *min_element(prop_combined.begin(), prop_combined.end());
    double max_weight_val = *max_element(prop_combined.begin(), prop_combined.end());
    cout << "Min value for weighting: " << min_weight_val << "\n";
    cout << "Max value for weighting: " << max_weight_val << "\n";

    // For demonstration, we simply set the weights map.
    weights = {{"CC", 1.0f}, {"CN", 10.0f}, {"NN", 100.0f}};

    // Read coordinate file
    string coord_file = "test.xyz";
    ifstream coordstream(coord_file);
    if (!coordstream.good()) {
        cout << "ERROR: Cannot open xyz file " << coord_file << "\n";
        exit(0);
    }
    
    // Read number of atoms
    string line;
    vector<string> line_contents;
    int ncontents;
    get_next_line(coordstream, line);
    ncontents = split_line(line, line_contents);
    int natoms = stoi(line_contents[0]);
    
    // Read box dimensions (assumes NON_ORTHO but orthorhombic box)
    xyz boxdims;
    get_next_line(coordstream, line);
    ncontents = split_line(line, line_contents);
    boxdims.x = stod(line_contents[1]);
    boxdims.y = stod(line_contents[5]);
    boxdims.z = stod(line_contents[9]);
    cout << "Read boxdims: " << boxdims.x << " " << boxdims.y << " " << boxdims.z << "\n";
    
    // Read atom coordinates
    vector<xyz> coords(natoms);
    for (int i = 0; i < natoms; i++) {
        get_next_line(coordstream, line);
        ncontents = split_line(line, line_contents);
        coords[i].atom_type = line_contents[0];
        coords[i].x = stod(line_contents[1]);
        coords[i].y = stod(line_contents[2]);
        coords[i].z = stod(line_contents[3]);
    }
    
    // Containers for clusters (each row holds the distances for one cluster)
    vector<vector<double>> cludists_2b_direct;
    vector<vector<double>> cludists_3b_direct;
    vector<vector<double>> cludists_4b_direct;
    vector<vector<double>> cludists_2b_trans;
    vector<vector<double>> cludists_3b_trans;
    vector<vector<double>> cludists_4b_trans;

    // Reserve some space if possible (optional)
    cludists_2b_direct.reserve(natoms);
    cludists_3b_direct.reserve(natoms);
    cludists_4b_direct.reserve(natoms);
    cludists_2b_trans.reserve(natoms);
    cludists_3b_trans.reserve(natoms);
    cludists_4b_trans.reserve(natoms);

    cout << "Beginning computation\n";
    struct timeval tv;
    gettimeofday(&tv, NULL);
    uint64_t start = tv.tv_sec*(uint64_t)1000000 + tv.tv_usec;

    // Parallelize the outer loop over i using OpenMP.
    // Each thread uses local containers which are merged later.
    vector<vector<vector<double>>> thread_2b_direct;
    vector<vector<vector<double>>> thread_3b_direct;
    vector<vector<vector<double>>> thread_4b_direct;
    vector<vector<vector<double>>> thread_2b_trans;
    vector<vector<vector<double>>> thread_3b_trans;
    vector<vector<vector<double>>> thread_4b_trans;
    
    int num_threads = omp_get_max_threads();
    thread_2b_direct.resize(num_threads);
    thread_3b_direct.resize(num_threads);
    thread_4b_direct.resize(num_threads);
    thread_2b_trans.resize(num_threads);
    thread_3b_trans.resize(num_threads);
    thread_4b_trans.resize(num_threads);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < natoms; i++) {
        int tid = omp_get_thread_num();

        // Local containers for this thread
        vector<double> dist_2b(1);
        vector<double> dist_2b_trans(1);
        vector<double> dist_3b(3);
        vector<double> dist_3b_trans(3);
        vector<double> dist_4b(6);
        vector<double> dist_4b_trans(6);
        
        for (int j = i + 1; j < natoms; j++) {
            // Use sorted key for lookup
            string key_ij = makeKey(coords[i].atom_type, coords[j].atom_type);
            double weight_ij = weights[key_ij];
            double unweighted_dist_ij = get_dist(boxdims, coords[i], coords[j], atom_map);
            double weighted_dist_ij = unweighted_dist_ij * weight_ij;
            dist_2b_trans[0] = transform(coords[i], coords[j], atom_map, rcout_2b, weighted_dist_ij);
            
            if (unweighted_dist_ij < rcout_2b) {
                dist_2b[0] = weighted_dist_ij;
                thread_2b_direct[tid].push_back(dist_2b);
                thread_2b_trans[tid].push_back(dist_2b_trans);
                
                if (unweighted_dist_ij >= rcout_3b)
                    continue;
                
                for (int k = j + 1; k < natoms; k++) {
                    string key_ik = makeKey(coords[i].atom_type, coords[k].atom_type);
                    string key_jk = makeKey(coords[j].atom_type, coords[k].atom_type);
                    double weight_ik = weights[key_ik];
                    double weight_jk = weights[key_jk];
                    
                    double unweighted_dist_ik = get_dist(boxdims, coords[i], coords[k], atom_map);
                    double unweighted_dist_jk = get_dist(boxdims, coords[j], coords[k], atom_map);
                    
                    if (unweighted_dist_ik >= rcout_3b || unweighted_dist_jk >= rcout_3b)
                        continue;
                    
                    dist_3b[0] = weighted_dist_ij;
                    dist_3b[1] = unweighted_dist_ik * weight_ik;
                    dist_3b[2] = unweighted_dist_jk * weight_jk;
                    
                    dist_3b_trans[0] = transform(coords[i], coords[j], atom_map, rcout_3b, dist_3b[0]);
                    dist_3b_trans[1] = transform(coords[i], coords[k], atom_map, rcout_3b, dist_3b[1]);
                    dist_3b_trans[2] = transform(coords[j], coords[k], atom_map, rcout_3b, dist_3b[2]);
                    
                    thread_3b_direct[tid].push_back(dist_3b);
                    thread_3b_trans[tid].push_back(dist_3b_trans);
                    
                    if (unweighted_dist_ij >= rcout_4b || unweighted_dist_ik >= rcout_4b || unweighted_dist_jk >= rcout_4b)
                        continue;
                    
                    for (int l = k + 1; l < natoms; l++) {
                        string key_il = makeKey(coords[i].atom_type, coords[l].atom_type);
                        string key_jl = makeKey(coords[j].atom_type, coords[l].atom_type);
                        string key_kl = makeKey(coords[k].atom_type, coords[l].atom_type);
                        
                        double weight_il = weights[key_il];
                        double weight_jl = weights[key_jl];
                        double weight_kl = weights[key_kl];
                        
                        double unweighted_dist_il = get_dist(boxdims, coords[i], coords[l], atom_map);
                        double unweighted_dist_jl = get_dist(boxdims, coords[j], coords[l], atom_map);
                        double unweighted_dist_kl = get_dist(boxdims, coords[k], coords[l], atom_map);
                        
                        if (unweighted_dist_il >= rcout_4b ||
                            unweighted_dist_jl >= rcout_4b ||
                            unweighted_dist_kl >= rcout_4b)
                            continue;
                        
                        dist_4b[0] = weighted_dist_ij;
                        dist_4b[1] = unweighted_dist_ik * weight_ik;
                        dist_4b[2] = unweighted_dist_jk * weight_jk;
                        dist_4b[3] = unweighted_dist_il * weight_il;
                        dist_4b[4] = unweighted_dist_jl * weight_jl;
                        dist_4b[5] = unweighted_dist_kl * weight_kl;
                        
                        dist_4b_trans[0] = transform(coords[i], coords[j], atom_map, rcout_4b, dist_4b[0]);
                        dist_4b_trans[1] = transform(coords[i], coords[k], atom_map, rcout_4b, dist_4b[1]);
                        dist_4b_trans[2] = transform(coords[j], coords[k], atom_map, rcout_4b, dist_4b[2]);
                        dist_4b_trans[3] = transform(coords[i], coords[l], atom_map, rcout_4b, dist_4b[3]);
                        dist_4b_trans[4] = transform(coords[j], coords[l], atom_map, rcout_4b, dist_4b[4]);
                        dist_4b_trans[5] = transform(coords[k], coords[l], atom_map, rcout_4b, dist_4b[5]);
                        
                        thread_4b_direct[tid].push_back(dist_4b);
                        thread_4b_trans[tid].push_back(dist_4b_trans);
                    } // l loop
                } // k loop
            } // if unweighted_dist_ij < rcout_2b
        } // j loop
    } // i loop

    // Merge thread–local results into the global containers
    for (int t = 0; t < num_threads; t++) {
        cludists_2b_direct.insert(cludists_2b_direct.end(), thread_2b_direct[t].begin(), thread_2b_direct[t].end());
        cludists_3b_direct.insert(cludists_3b_direct.end(), thread_3b_direct[t].begin(), thread_3b_direct[t].end());
        cludists_4b_direct.insert(cludists_4b_direct.end(), thread_4b_direct[t].begin(), thread_4b_direct[t].end());
        cludists_2b_trans.insert(cludists_2b_trans.end(), thread_2b_trans[t].begin(), thread_2b_trans[t].end());
        cludists_3b_trans.insert(cludists_3b_trans.end(), thread_3b_trans[t].begin(), thread_3b_trans[t].end());
        cludists_4b_trans.insert(cludists_4b_trans.end(), thread_4b_trans[t].begin(), thread_4b_trans[t].end());
    }

    gettimeofday(&tv, NULL);
    uint64_t end = tv.tv_sec*(uint64_t)1000000 + tv.tv_usec;
    uint64_t elapsed = end - start;
    cout << "@@@ Elapsed time (usec): " << elapsed << "\n";
    cout << "Processing complete. Preparing output.\n";
    cout << "   ...Calc done, printing results...\n";

    writeClusterData("2b_clu-r.txt", cludists_2b_direct);
    writeClusterData("3b_clu-r.txt", cludists_3b_direct);
    writeClusterData("4b_clu-r.txt", cludists_4b_direct);
    writeClusterData("2b_clu-s.txt", cludists_2b_trans);
    writeClusterData("3b_clu-s.txt", cludists_3b_trans);
    writeClusterData("4b_clu-s.txt", cludists_4b_trans);

    return 0;
}
