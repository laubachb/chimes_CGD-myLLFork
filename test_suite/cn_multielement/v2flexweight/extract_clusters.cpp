/*

This script takes in two .xyz(f) files and:
+ Produces distances lists for every 2, 3, and 4-body cluster
+ ...

Expects NON_ORTHO format, but assumes an orthorhombic box
assumes all atoms are the same type
*/

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<cmath>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <sys/time.h>
#include <thread>

using namespace std;

// Declare ALPHA as a global variable
double alpha = 0.1;
// Global weights map
map<string, float> weights;

struct xyz
{
    double x;
    double y;
    double z;
    string atom_type; // Added for multi-element
};

// Function to parse a comma-separated list into a vector
vector<string> parseStringList(const string& str)
{
    vector<string> result;
    stringstream ss(str);
    string token;
    while (getline(ss, token, ',')) {
        result.push_back(token);
    }
    return result;
}

vector<double> parseDoubleList(const string& str) 
{
    vector<double> result;
    stringstream ss(str);
    string token;
    while (getline(ss, token, ',')) {
        result.push_back(stod(token));
    }
    return result;
}

// Function to generate unique atom pairs from atom types, sorted alphabetically
vector<string> generateAtomPairs(const vector<string>& atom_types) 
{
    vector<string> atom_pairs;
    // Generate all unique combinations of atom pairs
    for (size_t i = 0; i < atom_types.size(); ++i) {
        for (size_t j = i; j < atom_types.size(); ++j) {  // Ensure i <= j for uniqueness
            string pair = atom_types[i] + atom_types[j];
            // Ensure the pair is in alphabetical order (CN instead of NC)
            if (atom_types[i] > atom_types[j]) {
                swap(pair[0], pair[1]);
            }
            atom_pairs.push_back(pair);
        }
    }
    return atom_pairs;
}

// Function to read config file into a key-value map
unordered_map<string, string> readConfig(const string& filename) 
{
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
  
int split_line(string line, vector<string> & items)
{
    // Break a line up into tokens based on space separators.
    // Returns the number of tokens parsed.
    
    string       contents;
    stringstream sstream;

    // Strip comments beginining with ! or ## and terminal new line

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

    while ( sstream >> contents ) 
        items.push_back(contents);

    return items.size();
}

bool get_next_line(istream& str, string & line)
{
    // Read a line and return it, with error checking.
    
        getline(str, line);

        if(!str)
            return false;
    
    return true;
}

// Function to precompute weights for all atom pairs
void compute_weights(const map<string, pair<float, double>>& atom_map, float max_prop, float min_prop, float alpha) 
{
    for (const auto& entry : atom_map) {
        string atom_pair = entry.first;
        float prop = entry.second.first; // Extract lambda value
        std::cout << "Combined property for atom pair " << atom_pair << " = " << prop << endl;
        
        // Compute and store weight
        float weight = alpha + (1 - alpha) * ((prop - min_prop) / (max_prop - min_prop));
        std::cout << "Weight for atom pair " << atom_pair << " = " << weight << endl;
        weights[atom_pair] = weight;
    }
}

double get_dist(xyz box, xyz a1, xyz a2, map<string, pair<float, double>> atom_map)
{
    // Obtain atom type identities
    string atom_type1 = a1.atom_type;
    string atom_type2 = a2.atom_type;
    string atom_pair = atom_type1 + atom_type2;
    
    double rcin = 0.0;
    // Added for multi-element
    // Get rcin value for the concatenated pair
    if (atom_map.find(atom_pair) != atom_map.end()) {
        rcin = atom_map[atom_pair].second;
    } else {
        std::cout << "Atom pair " << atom_pair << " not found in the map." << endl;
        return -1; // Return an error if the atom pair is not found
    }

    double dx = a1.x - a2.x; dx -= box.x * round(dx / box.x);
    double dy = a1.y - a2.y; dy -= box.y * round(dy / box.y);
    double dz = a1.z - a2.z; dz -= box.z * round(dz / box.z);

    double dist = sqrt(dx * dx + dy * dy + dz * dz);

    if (dist < rcin)
    {
        std::cout << "Something went wrong between atoms: " << endl;
        std::cout << "A: " << a1.x << " " << a1.y << " " << a1.z << endl;
        std::cout << "B: " << a2.x << " " << a2.y << " " << a2.z << endl;
        std::cout << dx << " " << dy << " " << dz << endl;
        dist = 100000000;
    }

    return dist;  // Return only the distance
}

// transform(coords[i], coords[j], atom_map, rcout_2b, dist_2b[0])
double transform(xyz a1, xyz a2, map<string, pair<float, double>> atom_map, double rcout, double rij)
{
    // Obtain atom type identities
    string atom_type1 = a1.atom_type;
    string atom_type2 = a2.atom_type;
    string atom_pair = atom_type1 + atom_type2;
    
    float lambda = 0.0f;
    double rcin = 0.0;
    
    // Get the lambda and rcin values for the concatenated pair
    if (atom_map.find(atom_pair) != atom_map.end()) {
        // Access both lambda and rcin from the pair
        lambda = atom_map[atom_pair].first;
        rcin = atom_map[atom_pair].second;

    } else {
        std::cout << "Atom pair " << atom_pair << " not found in the map." << endl;
        return -1; // Return an error if the atom pair is not found
    }
    
    double x_min = exp(-1*rcin/lambda);
    double x_max = exp(-1*rcout/lambda);
    double x_avg   = 0.5 * (x_max + x_min);
    double x_diff  = 0.5 * (x_max - x_min);
    
    x_diff *= -1.0; // Special for Morse style
    
    return (exp(-1*rij/lambda) - x_avg)/x_diff;
}

void writeClusterData(const string& filename, const vector<vector<double>>& data) 
{
    ofstream ofs(filename);
    if (!ofs) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }
    
    ostringstream buffer;
    
    for (const auto& row : data) {
        for (size_t j = 0; j < row.size(); j++) {
            buffer << row[j];
            if (j < row.size() - 1) buffer << " ";  // Add space between elements
        }
        buffer << "\n";
    }

    ofs << buffer.str(); // Single large write operation
    ofs.close();
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
    
    // Read atom types from the config (e.g., "C,N,O")
    vector<string> atom_types = parseStringList(config["ATOM_TYPES"]);
    // Generate atom pairs from atom types
    vector<string> atom_pairs = generateAtomPairs(atom_types);
    vector<double> rcin_list = parseDoubleList(config["RCIN_LIST"]);
    vector<double> lambda_set;
    vector<double> property_set;
    if (config.find("LAMBDA_SET") != config.end()) {
        lambda_set = parseDoubleList(config["LAMBDA_SET"]);
    }
    if (config.find("PROPERTY_SET") != config.end()) {
        property_set = parseDoubleList(config["PROPERTY_SET"]);
    }

    // Compute combined properties
    // Ensure property_set is the same size as atom_types
    if (property_set.size() != atom_types.size()) {
        cerr << "Error: PROPERTY_SET size does not match ATOM_TYPES size.\n";
        return 1;
    }
    // Map atom types to their corresponding property values
    unordered_map<string, double> property_map;
    for (size_t i = 0; i < atom_types.size(); ++i) {
        property_map[atom_types[i]] = property_set[i];
    }
    // Compute prop_combined and x_combined
    vector<double> prop_combined;

    for (const auto& pair : atom_pairs) {
        // Ensure the pair has exactly two characters (assuming valid input like "CN")
        if (pair.length() != 2) {
            cerr << "Error: Invalid atom pair format: " << pair << endl;
            return 1;
        }
        string atom1(1, pair[0]);  // First character
        string atom2(1, pair[1]);  // Second character

        // Get corresponding properties
        if (property_map.find(atom1) == property_map.end() || property_map.find(atom2) == property_map.end()) {
            cerr << "Error: Atom type not found in PROPERTY_SET: " << pair << endl;
            return 1;
        }

        double prop1 = property_map[atom1];
        double prop2 = property_map[atom2];

        // Compute sqrt(prop1 * prop2) and store in prop_combined
        double prop_value = sqrt(prop1 * prop2);
        prop_combined.push_back(prop_value);
    }

    // Output parsed values
    std::cout << "Using config file: " << configFile << endl;
    std::cout << "rcout_2b: " << rcout_2b << endl;
    std::cout << "rcout_3b: " << rcout_3b << endl;
    std::cout << "rcout_4b: " << rcout_4b << endl;

    std::cout << "Atom Pairs: ";
    for (const auto& pair : atom_pairs) std::cout << pair << " ";
    std::cout << endl;

    // Print lambda_set only if it has elements
    if (!lambda_set.empty()) {
        std::cout << "Lambda Set: ";
        for (const auto& lambda : lambda_set) {
            std::cout << lambda << " ";
        }
        std::cout << endl;
    } else {
        std::cout << "Lambda set not provided\n";
        std::cout << "!!!Exiting code!!!\n";
    }

    std::cout << "RCIN List: ";
    for (const auto& rcin : rcin_list) std::cout << rcin << " ";
    std::cout << endl;

    // Read ALPHA if it exists in the config
    if (config.find("ALPHA") != config.end()) {
        alpha = stod(config["ALPHA"]);
    } else {
        cerr << "Warning: ALPHA not found in config file. Using default value: " << alpha << endl;
    }
    
    /////////////////////////////////////////////
    // Hard-coded for now, for a single atom type: rcutin, rcut out, morse lambda
    /////////////////////////////////////////////

    // Create the dictionary (map) with std::pair<float, double> as the value type
    map<string, pair<float, double>> atom_map;
    for (size_t i = 0; i < atom_pairs.size(); ++i) {
        atom_map[atom_pairs[i]] = make_pair(prop_combined[i], rcin_list[i]);
    }

    // Get the maximum lambda value
    double min_rcin = *min_element(rcin_list.begin(), rcin_list.end());
    double min_weight_val = *min_element(prop_combined.begin(), prop_combined.end());
    double max_weight_val = *max_element(prop_combined.begin(), prop_combined.end());
    std::cout << "Min value for weighting: " << min_weight_val << std::endl;
    std::cout << "Max value for weighting: " << max_weight_val << std::endl;

    // Compute the weight list
    //compute_weights(atom_map, max_weight_val, min_weight_val, alpha);
    weights = {{"CC", 1.0f}, {"CN", 10.0f}, {"NN", 100.0f}};

    /////////////////////////////////////////////
    // Read file name
    /////////////////////////////////////////////
    
    string coord_file = "test.xyz";
    
    /////////////////////////////////////////////
    // Open file to read in coordinates
    /////////////////////////////////////////////
    
    ifstream coordstream;
    coordstream.open(coord_file);
    if (!coordstream.good())
    {
        std::cout << "ERROR: Cannot open xyz file " << coord_file << endl;
        exit(0);
    }
    
    
    /////////////////////////////////////////////
    // Read in number of atoms
    /////////////////////////////////////////////
    
    
    string  line;
    vector<string> line_contents;
    int     ncontents;

    get_next_line(coordstream, line);
    ncontents = split_line(line, line_contents);
    int natoms = stoi(line_contents[0]);
    
    
    /////////////////////////////////////////////
    // read the boxdims
    /////////////////////////////////////////////
    
    // Defined struct xyz

    xyz     boxdims;
    
    get_next_line(coordstream, line);
    ncontents = split_line(line, line_contents);
    
    boxdims.x = stod(line_contents[1]);
    boxdims.y = stod(line_contents[5]);
    boxdims.z = stod(line_contents[9]);
    
    std::cout << "Read boxdims: " << boxdims.x << " " << boxdims.y << " " << boxdims.z << endl;
    
    
    /////////////////////////////////////////////
    // Read the atom coordinates
    /////////////////////////////////////////////

    vector<xyz> coords(natoms);
    
    for (int i=0; i<natoms; i++)
    {
        get_next_line(coordstream, line);
        ncontents = split_line(line, line_contents);
        
        coords[i].atom_type = string(line_contents[0]);
        coords[i].x = stod(line_contents[1]);
        coords[i].y = stod(line_contents[2]);
        coords[i].z = stod(line_contents[3]);

    }
    
    
    /////////////////////////////////////////////
    // Extract all 2-, 3-, and 4-body clusters
    /////////////////////////////////////////////
    
    vector<vector<double > > cludists_2b_direct;   // dim 1: cluster index; dim 2, distance between atom pair
    vector<vector<double > > cludists_3b_direct;   // dim 1: cluster index; dim 2, distance between atom pairs
    vector<vector<double > > cludists_4b_direct;   // dim 1: cluster index; dim 2, distance between atom pairs
    vector<vector<double > > cludists_2b_trans;    // dim 1: cluster index; dim 2, distance between atom pair
    vector<vector<double > > cludists_3b_trans;    // dim 1: cluster index; dim 2, distance between atom pairs
    vector<vector<double > > cludists_4b_trans;    // dim 1: cluster index; dim 2, distance between atom pairs
    
    // ... define get dist
    
    vector<double> dist_2b(1);
    vector<double> dist_3b(3);
    vector<double> dist_4b(6);
    vector<double> dist_2b_trans(1);
    vector<double> dist_3b_trans(3);
    vector<double> dist_4b_trans(6);

    std::cout << "Beginning computation\n";

	/* get start timestamp */
 	struct timeval tv;
    	gettimeofday(&tv,NULL);
    	uint64_t start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
 
for (int i = 0; i < natoms; i++)
{
    for (int j = i + 1; j < natoms; j++)
    {
        string atom_pair_ij = coords[i].atom_type + coords[j].atom_type;
        double weight_ij = weights[atom_pair_ij];

        double unweighted_dist_ij = get_dist(boxdims, coords[i], coords[j], atom_map);
        double weighted_dist_ij = unweighted_dist_ij * weight_ij;
        dist_2b_trans[0] = transform(coords[i], coords[j], atom_map, rcout_2b, weighted_dist_ij);
        
        if (unweighted_dist_ij < rcout_2b)
        {
            dist_2b[0] = weighted_dist_ij;
            cludists_2b_direct.push_back(dist_2b);
            cludists_2b_trans.push_back(dist_2b_trans);

            if (unweighted_dist_ij >= rcout_3b)
                continue;

            for (int k = j + 1; k < natoms; k++)
            {
                string atom_pair_ik = coords[i].atom_type + coords[k].atom_type;
                string atom_pair_jk = coords[j].atom_type + coords[k].atom_type;

                double weight_ik = weights[atom_pair_ik];
                double weight_jk = weights[atom_pair_jk];

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
                
                cludists_3b_direct.push_back(dist_3b);
                cludists_3b_trans.push_back(dist_3b_trans);

                if (unweighted_dist_ij >= rcout_4b || unweighted_dist_ik >= rcout_4b || unweighted_dist_jk >= rcout_4b)
                    continue;

                for (int l = k + 1; l < natoms; l++)
                {
                    string atom_pair_il = coords[i].atom_type + coords[l].atom_type;
                    string atom_pair_jl = coords[j].atom_type + coords[l].atom_type;
                    string atom_pair_kl = coords[k].atom_type + coords[l].atom_type;

                    double weight_il = weights[atom_pair_il];
                    double weight_jl = weights[atom_pair_jl];
                    double weight_kl = weights[atom_pair_kl];

                    double unweighted_dist_il = get_dist(boxdims, coords[i], coords[l], atom_map);
                    double unweighted_dist_jl = get_dist(boxdims, coords[j], coords[l], atom_map);
                    double unweighted_dist_kl = get_dist(boxdims, coords[k], coords[l], atom_map);
                    
                    if (unweighted_dist_il >= rcout_4b || unweighted_dist_jl >= rcout_4b || unweighted_dist_kl >= rcout_4b)
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
                    
                    cludists_4b_direct.push_back(dist_4b);
                    cludists_4b_trans.push_back(dist_4b_trans);
                }
            }
        }
    }
}
	/* get elapsed time */
    	gettimeofday(&tv,NULL);
    	uint64_t end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    	uint64_t elapsed = end - start;

	std::cout << "@@@ Elapsed time (usec): " << elapsed << "\n";
	std::cout << "Processing complete.  Preparing output.\n";
    std::cout << "   ...Calc done, printing results..." << endl;
    
    /////////////////////////////////////////////
    // Print UNSORTED cluster distances (rij)
    /////////////////////////////////////////////
    
    writeClusterData("2b_clu-r.txt", cludists_2b_direct);
    writeClusterData("3b_clu-r.txt", cludists_3b_direct);
    writeClusterData("4b_clu-r.txt", cludists_4b_direct);

    writeClusterData("2b_clu-s.txt", cludists_2b_trans);
    writeClusterData("3b_clu-s.txt", cludists_3b_trans);
    writeClusterData("4b_clu-s.txt", cludists_4b_trans);

    return 0;
}