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

using namespace std;

// Declare ALPHA as a global variable
double ALPHA = 0.1;

struct xyz
{
    double x;
    double y;
    double z;
    string atom_type; // Added for multi-element
};

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

// Function to generate unique atom pairs from atom types, sorted alphabetically
vector<string> generateAtomPairs(const vector<string>& atom_types) {
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


double get_dist(xyz box, xyz a1, xyz a2, map<string, pair<float, double>> atom_map, float max_prop)
{
    // Obtain atom type identities
    string atom_type1 = a1.atom_type;
    string atom_type2 = a2.atom_type;
    string atom_pair = atom_type1 + atom_type2;
    
    float prop = 0.0f;
    double rcin = 0.0;
    
    // Get the lambda and rcin values for the concatenated pair
    if (atom_map.find(atom_pair) != atom_map.end()) {
        // Access both lambda and rcin from the pair
        prop = atom_map[atom_pair].first;
        rcin = atom_map[atom_pair].second;

    } else {
        cout << "Atom pair " << atom_pair << " not found in the map." << endl;
        return -1; // Return an error if the atom pair is not found
    }
    
    float weight = prop / max_prop; // Calculate weight
    
    double dx = a1.x - a2.x; dx -= box.x*round(dx/box.x);
    double dy = a1.y - a2.y; dy -= box.y*round(dy/box.y);
    double dz = a1.z - a2.z; dz -= box.z*round(dz/box.z);
    
    
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    
    if (dist < rcin)
    {
        cout << "Something went wrong between atoms: " << endl;
        cout << "A: " << a1.x << " " << a1.y << " " << a1.z << endl;
        cout << "B: " << a2.x << " " << a2.y << " " << a2.z << endl;
        cout << dx << " " << dy << " " << dz << endl;
        dist = 100;
    }
    
    return weight * dist;
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
        cout << "Atom pair " << atom_pair << " not found in the map." << endl;
        return -1; // Return an error if the atom pair is not found
    }
    
    double x_min = exp(-1*rcin/lambda);
    double x_max = exp(-1*rcout/lambda);
    double x_avg   = 0.5 * (x_max + x_min);
    double x_diff  = 0.5 * (x_max - x_min);
    
    x_diff *= -1.0; // Special for Morse style
    
    return (exp(-1*rij/lambda) - x_avg)/x_diff;
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
    cout << "Using config file: " << configFile << endl;
    cout << "rcout_2b: " << rcout_2b << endl;
    cout << "rcout_3b: " << rcout_3b << endl;
    cout << "rcout_4b: " << rcout_4b << endl;

    cout << "Atom Pairs: ";
    for (const auto& pair : atom_pairs) cout << pair << " ";
    cout << endl;

    // Print lambda_set only if it has elements
    if (!lambda_set.empty()) {
        cout << "Lambda Set: ";
        for (const auto& lambda : lambda_set) {
            cout << lambda << " ";
        }
        cout << endl;
    } else {
        cout << "Lambda set not provided\n";
    }

    cout << "RCIN List: ";
    for (const auto& rcin : rcin_list) cout << rcin << " ";
    cout << endl;

    // Read ALPHA if it exists in the config
    if (config.find("ALPHA") != config.end()) {
        ALPHA = stod(config["ALPHA"]);
    } else {
        cerr << "Warning: ALPHA not found in config file. Using default value: " << ALPHA << endl;
    }
    
    /////////////////////////////////////////////
    // Hard-coded for now, for a single atom type: rcutin, rcut out, morse lambda
    /////////////////////////////////////////////

    // If lambda_set is undefined, use prop_combined instead
    if (lambda_set.empty()) {
        lambda_set = prop_combined;
        cout << "Setting lambda set to combined property set\n";
    }

    // Create the dictionary (map) with std::pair<float, double> as the value type
    map<string, pair<float, double>> atom_map;
    for (size_t i = 0; i < atom_pairs.size(); ++i) {
        atom_map[atom_pairs[i]] = make_pair(lambda_set[i], rcin_list[i]);
    }

    // Get the maximum lambda value
    double min_rcin = *min_element(rcin_list.begin(), rcin_list.end());
    double max_lambda = *max_element(lambda_set.begin(), lambda_set.end());
    cout << "Max value for weighting: " << max_lambda << std::endl;

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
        cout << "ERROR: Cannot open xyz file " << coord_file << endl;
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
    
    cout << "Read boxdims: " << boxdims.x << " " << boxdims.y << " " << boxdims.z << endl;
    
    
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
 
    for (int i=0; i<natoms; i++)
    {
        for (int j=i+1; j<natoms; j++)
        {

            //cout << "Calculating distances between atoms with index and coordinates:" << endl;
            //cout << i << " " << coords[i].x << " " << coords[i].y << " " << coords[i].z << endl;
            //cout << j << " " << coords[j].x << " " << coords[j].y << " " << coords[j].z << endl;
            
            dist_2b[0] = get_dist(boxdims, coords[i], coords[j], atom_map, max_lambda);
            dist_2b_trans[0] = transform(coords[i], coords[j], atom_map, rcout_2b, dist_2b[0]);
            
            //cout << dist_2b[0] << endl;
            
            if (dist_2b[0] < rcout_2b)
            {
                //cout << dist_2b[0] << endl;
                
                cludists_2b_direct.push_back(dist_2b);
                cludists_2b_trans.push_back(dist_2b_trans);
                
                if (dist_2b[0] >= rcout_3b)
                    continue;
                
                for (int k=j+1; k<natoms; k++)
                {
                    dist_3b[0] = dist_2b[0]; 
                    dist_3b_trans[0] = dist_2b_trans[0];
                    dist_3b[1] = get_dist(boxdims, coords[i], coords[k], atom_map, max_lambda);
                    dist_3b_trans[1] = transform(coords[i], coords[k], atom_map, rcout_3b, dist_3b[1]);
                    
                    if (dist_3b[1] >= rcout_3b)
                        continue;
                    
                    dist_3b[2] = get_dist(boxdims, coords[j], coords[k], atom_map, max_lambda);
                    dist_3b_trans[2] = transform(coords[j], coords[k], atom_map, rcout_3b, dist_3b[2]);
                    
                    if (dist_3b[2] >= rcout_3b)
                        continue;
                    
                    cludists_3b_direct.push_back(dist_3b);
                    cludists_3b_trans.push_back(dist_3b_trans);
                    
                    if (dist_3b[0] >= rcout_4b)
                        continue; 
                    if (dist_3b[1] >= rcout_4b)
                        continue; 
                    if (dist_3b[2] >= rcout_4b)
                        continue; 
                    
                    
                    for (int l=k+1; l<natoms; l++)
                    {  
                        dist_4b[0] = dist_3b[0]; // ij
                        dist_4b[1] = dist_3b[1]; // ik
                        dist_4b[2] = dist_3b[2]; // jk
                        
                        dist_4b_trans[0] = dist_3b_trans[0]; // ij
                        dist_4b_trans[1] = dist_3b_trans[1]; // ik
                        dist_4b_trans[2] = dist_3b_trans[2]; // jk
                        
                        dist_4b[3] = get_dist(boxdims, coords[i], coords[l], atom_map, max_lambda);
                        dist_4b_trans[3] = transform(coords[i], coords[l], atom_map, rcout_4b, dist_4b[3]);
                    
                        if (dist_4b[3] >= rcout_4b)
                            continue;       

                        dist_4b[4] = get_dist(boxdims, coords[j], coords[l], atom_map, max_lambda);
                        dist_4b_trans[4] = transform(coords[j], coords[l], atom_map, rcout_4b, dist_4b[4]);
                    
                        if (dist_4b[4] >= rcout_4b)
                            continue;      
                        
                        dist_4b[5] = get_dist(boxdims, coords[k], coords[l], atom_map, max_lambda);
                        dist_4b_trans[5] = transform(coords[k], coords[l], atom_map, rcout_4b, dist_4b[5]);
                    
                        if (dist_4b[5] >= rcout_4b)
                            continue;  
                        
                        cludists_4b_direct.push_back(dist_4b); 
                        cludists_4b_trans.push_back(dist_4b_trans); 
                    }
                                        
                }
            }
        }
    }
    
    cout << "   ...Calc done, printing results..." << endl;
    
    /////////////////////////////////////////////
    // Print UNSORTED cluster distances (rij)
    /////////////////////////////////////////////
    
    ofstream ofstream_2b_r;
    ofstream ofstream_3b_r;
    ofstream ofstream_4b_r;
    
    ofstream_2b_r.open("2b_clu-r.txt");
    ofstream_3b_r.open("3b_clu-r.txt");
    ofstream_4b_r.open("4b_clu-r.txt");
    
    for (int i=0; i<cludists_2b_direct.size(); i++)
        ofstream_2b_r << cludists_2b_direct[i][0] << endl;
    
    for (int i=0; i<cludists_3b_direct.size(); i++)
    {
        for (int j=0; j<3; j++)
            ofstream_3b_r << cludists_3b_direct[i][j] <<  " ";
        ofstream_3b_r << endl;
    }
    
    for (int i=0; i<cludists_4b_direct.size(); i++)
    {
        for (int j=0; j<6; j++)
            ofstream_4b_r << cludists_4b_direct[i][j] <<  " "; 
        ofstream_4b_r << endl;  
    }
    
    ofstream_2b_r.close();
    ofstream_3b_r.close();
    ofstream_4b_r.close();
    
    
    /////////////////////////////////////////////
    // Print UNSORTED cluster TRANSFORMED distances (sij)
    /////////////////////////////////////////////    
    
    
    ofstream ofstream_2b_s;
    ofstream ofstream_3b_s;
    ofstream ofstream_4b_s;
    
    ofstream_2b_s.open("2b_clu-s.txt");
    ofstream_3b_s.open("3b_clu-s.txt");
    ofstream_4b_s.open("4b_clu-s.txt");
    
    for (int i=0; i<cludists_2b_trans.size(); i++)
        ofstream_2b_s << cludists_2b_trans[i][0] << endl;
    
    for (int i=0; i<cludists_3b_trans.size(); i++)
    {
        for (int j=0; j<3; j++)
            ofstream_3b_s << cludists_3b_trans[i][j] <<  " ";
        ofstream_3b_s << endl;
    }
    
    for (int i=0; i<cludists_4b_trans.size(); i++)
    {
        for (int j=0; j<6; j++)
            ofstream_4b_s << cludists_4b_trans[i][j] <<  " "; 
        ofstream_4b_s << endl;  
    }

    ofstream_2b_s.close();
    ofstream_3b_s.close();
    ofstream_4b_s.close();  
}
