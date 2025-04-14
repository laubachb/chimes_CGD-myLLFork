/*
   CUDA version for computing 2-, 3-, and 4-body clusters.
   (The host code still reads the configuration file and coordinate file.)
   
   NOTE:
   - Atom types are converted from strings to integer codes.
   - The “atom_map” parameters (rcin, lambda, weight) are stored in arrays of size ntypes*ntypes.
   - For simplicity, the kernels write clusters into preallocated output arrays
     (whose maximum sizes are nchoose2, nchoose3, and nchoose4).
   - Each kernel uses an atomic counter to record how many clusters passed the cutoff.
   - This code must be compiled with nvcc (e.g., nvcc -O3 your_cuda_code.cu -o your_program)
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <sys/time.h>
// CUDA runtime
#include <cuda_runtime.h>

using namespace std;
const long long BUFFER_SIZE_3B = 1000000; // 1 million clusters per buffer
const long long BUFFER_SIZE_4B = 1000000;
//---------------------------------------------------------------------
// Data structure for coordinates. (Note: atom_type is now an int code.)
struct xyz {
    double x;
    double y;
    double z;
    int atom_type;
};

int split_line(string line, vector<string> & items){
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

bool get_next_line(istream& str, string & line){
    // Read a line and return it, with error checking.
    
        getline(str, line);

        if(!str)
            return false;
    
    return true;
}

//---------------------------------------------------------------------
// Read the configuration file into an unordered_map<string,string>
unordered_map<string, string> readConfig(const string &filename) {
    unordered_map<string, string> config;
    ifstream file(filename);
    if (!file) {
        cerr << "Error: Cannot open config file: " << filename << endl;
        exit(1);
    }
    string line;
    while(getline(file, line)) {
        // Skip empty lines and comments
        if(line.empty() || line[0] == '#')
            continue;
        size_t pos = line.find('=');
        if (pos == string::npos)
            continue;
        string key = line.substr(0, pos);
        string value = line.substr(pos+1);
        // Trim spaces
        key.erase(key.find_last_not_of(" \t")+1);
        key.erase(0, key.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t")+1);
        value.erase(0, value.find_first_not_of(" \t"));
        config[key] = value;
    }
    return config;
}

//---------------------------------------------------------------------
// Parse a comma-separated list into a vector of strings.
vector<string> parseStringList(const string &str) {
    vector<string> tokens;
    stringstream ss(str);
    string token;
    while(getline(ss, token, ',')) {
        // Trim spaces
        token.erase(token.find_last_not_of(" \t")+1);
        token.erase(0, token.find_first_not_of(" \t"));
        tokens.push_back(token);
    }
    return tokens;
}

//---------------------------------------------------------------------
// Parse a comma-separated list into a vector of doubles.
vector<double> parseDoubleList(const string &str) {
    vector<double> tokens;
    stringstream ss(str);
    string token;
    while(getline(ss, token, ',')) {
        token.erase(token.find_last_not_of(" \t")+1);
        token.erase(0, token.find_first_not_of(" \t"));
        tokens.push_back(stod(token));
    }
    return tokens;
}

__device__ inline double periodic_diff(double d, double box_length) {
    return d - box_length * round(d / box_length);
}

// ---------------------------------------------------------------------
// Kernel for 2-body clusters
// For each unique pair (i,j) with i < j we compute the (unweighted) distance,
// multiply by the appropriate weight, and if the unweighted distance is less than rcout_2b,
// we store the weighted distance and its transformed value.
struct PairData {
    int i, j;          // Atom indices (i < j)
    double dist;       // Distance after periodic boundary conditions
    int type_a, type_b; // Sorted atom types
    double rcin;       // Inner cutoff for this pair type
    float lambda;      // Lambda parameter for transformation
    float weight;      // Weight parameter (if used)
};

__device__ bool check_pair(const PairData& pair, double rcout) {
    double dist = pair.dist;
    if (dist < pair.rcin) dist = 1e8;
    return (dist * pair.weight) >= rcout;
}

__device__ double compute_trans(const PairData& pair, double rcout) {
    double x_min = exp(-pair.rcin / pair.lambda);
    double x_max = exp(-rcout / pair.lambda);
    double x_avg = 0.5 * (x_min + x_max);
    double x_diff = -0.5 * (x_max - x_min);
    return (exp(-(pair.dist * pair.weight) / pair.lambda) - x_avg) / x_diff;
}

__global__ void kernel_pairs(
    const xyz* d_coords, int natoms, xyz box,
    const double* d_rcin, const float* d_lambda, const float* d_weight,
    PairData* d_pairs, int ntypes
) {
    long long totalPairs = ((long long)natoms * (natoms - 1)) / 2;
    long long idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= totalPairs) return;

    // Unrank idx into (i, j) where i < j
    int i, j;
    long long rem = idx;
    for (i = 0; i < natoms - 1; i++) {
        long long count_i = natoms - 1 - i;
        if (rem < count_i) {
            j = i + 1 + rem;
            break;
        }
        rem -= count_i;
    }

    // Compute pair data
    int type_i = d_coords[i].atom_type;
    int type_j = d_coords[j].atom_type;
    int a = min(type_i, type_j);
    int b = max(type_i, type_j);

    float weight = d_weight[a * ntypes + b];
    double rcin = d_rcin[a * ntypes + b];
    float lambda = d_lambda[a * ntypes + b];

    double dx = d_coords[i].x - d_coords[j].x;
    double dy = d_coords[i].y - d_coords[j].y;
    double dz = d_coords[i].z - d_coords[j].z;
    dx = periodic_diff(dx, box.x);
    dy = periodic_diff(dy, box.y);
    dz = periodic_diff(dz, box.z);
    double dist = sqrt(dx*dx + dy*dy + dz*dz);

    // Store pair data
    d_pairs[idx] = {i, j, dist, a, b, rcin, lambda, weight};
}

__device__ long long get_pair_index(int i, int j, int natoms) {
    // Manual swap implementation for device code
    if (i > j) {
        int temp = i;
        i = j;
        j = temp;
    }
    return ((long long)i * (2 * natoms - i - 1)) / 2 + (j - i - 1);
}

__global__ void kernel_2b(const xyz* d_coords, int natoms, xyz box, double rcout_2b, 
                            const double* d_rcin, const float* d_lambda, const float* d_weight,
                            double* d_2b_direct, double* d_2b_trans, int* d_count, int ntypes) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if(i >= natoms || j >= natoms) return;
    if(j <= i) return;  // process only i < j

    int type_i = d_coords[i].atom_type;
    int type_j = d_coords[j].atom_type;
    // get sorted indices
    int a = (type_i < type_j) ? type_i : type_j;
    int b = (type_i < type_j) ? type_j : type_i;
    // float weight = d_weight[a * ntypes + b];
    float weight = 1.0;
    double rcin = d_rcin[a * ntypes + b];
    
    double dx = d_coords[i].x - d_coords[j].x;
    double dy = d_coords[i].y - d_coords[j].y;
    double dz = d_coords[i].z - d_coords[j].z;
    // printf("Coords: %f, %f, %f\n", d_coords[i].x, d_coords[i].y, d_coords[i].z);
    dx = periodic_diff(dx, box.x);
    dy = periodic_diff(dy, box.y);
    dz = periodic_diff(dz, box.z);
    // printf("Distances: %f, %f, %f\n", dx, dy, dz);
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    // printf("distance between atoms (rcin = %f): %f, rcout_2b: %f\n", rcin, dist, rcout_2b);
    // if (dist < rcout_2b) {
    //     printf("Cluster formed!\n");
    // }
    if(dist < rcin) { 
        dist = 1e8;  // flag error condition 
    }
    
    double weighted_dist = dist * 1;
    
    // Compute the transformation using the Morse–style formula:
    double lambda_val = d_lambda[a * ntypes + b];
    double x_min_val = exp(-rcin / lambda_val);
    double x_max_val = exp(-rcout_2b / lambda_val);
    double x_avg = 0.5 * (x_min_val + x_max_val);
    double x_diff = -0.5 * (x_max_val - x_min_val);
    double t_val = (exp(-weighted_dist / lambda_val) - x_avg) / x_diff;
    
    // Only if the unweighted distance is within cutoff do we store:
    if(dist < rcout_2b) {
        int index = atomicAdd(d_count, 1);
        d_2b_direct[index] = weighted_dist;
        d_2b_trans[index] = t_val;
        printf("%i \n", index);
    }
}

// ---------------------------------------------------------------------
// Kernel for 3-body clusters (triplets)
// We flatten the triplet (i,j,k) with i < j < k into a single linear index.
// The unranking routine converts a thread’s linear index to (i,j,k).
__global__ void kernel_3b(
    const PairData* d_pairs, int natoms, double rcout_3b,
    double* d_3b_direct, double* d_3b_trans, int* d_count, int ntypes,
    long long start_idx, long long end_idx
) {
    long long idx = blockIdx.x * blockDim.x + threadIdx.x;
    long long global_idx = start_idx + idx;
    if (global_idx >= end_idx) return;

    // Unrank global_idx into unique (i,j,k)
    int i, j, k;
    long long rem = global_idx;
    for(i = 0; i < natoms - 2; i++){
        long long count_i = ((long long)(natoms - i - 1)*(natoms - i - 2))/2;
        if(rem < count_i) break;
        rem -= count_i;
    }
    for(j = i+1; j < natoms - 1; j++){
        long long count_j = natoms - j - 1;
        if(rem < count_j) break;
        rem -= count_j;
    }
    k = j + 1 + rem;

    // Get pair indices
    long long idx_ij = get_pair_index(i, j, natoms);
    long long idx_ik = get_pair_index(i, k, natoms);
    long long idx_jk = get_pair_index(j, k, natoms);

    PairData pair_ij = d_pairs[idx_ij];
    PairData pair_ik = d_pairs[idx_ik];
    PairData pair_jk = d_pairs[idx_jk];

    // Check pair cutoffs
    if (check_pair(pair_ij, rcout_3b)) return;
    if (check_pair(pair_ik, rcout_3b)) return;
    if (check_pair(pair_jk, rcout_3b)) return;

    // Compute transformed values
    double trans_ij = compute_trans(pair_ij, rcout_3b);
    double trans_ik = compute_trans(pair_ik, rcout_3b);
    double trans_jk = compute_trans(pair_jk, rcout_3b);

    // Store results with buffer bounds checking
    int index = atomicAdd(d_count, 1);
    if(index < BUFFER_SIZE_3B) {
        d_3b_direct[index*3 + 0] = pair_ij.dist * pair_ij.weight;
        d_3b_direct[index*3 + 1] = pair_ik.dist * pair_ik.weight;
        d_3b_direct[index*3 + 2] = pair_jk.dist * pair_jk.weight;
        d_3b_trans[index*3 + 0] = trans_ij;
        d_3b_trans[index*3 + 1] = trans_ik;
        d_3b_trans[index*3 + 2] = trans_jk;
    }
}

// ---------------------------------------------------------------------
// Kernel for 4-body clusters (quadruplets)
// Here we “unrank” a linear index into a unique (i,j,k,l) combination (with i < j < k < l).
__global__ void kernel_4b(
    const PairData* d_pairs, int natoms, double rcout_4b,
    double* d_4b_direct, double* d_4b_trans, int* d_count, int ntypes,
    long long start_idx, long long end_idx
) {
    long long idx = blockIdx.x * blockDim.x + threadIdx.x;
    long long global_idx = start_idx + idx;
    if (global_idx >= end_idx) return;

    // Unrank global_idx into unique (i,j,k,l)
    int i, j, k, l;
    long long rem = global_idx;
    for(i = 0; i < natoms - 3; i++){
        long long count_i = ((long long)(natoms - i - 1)*(natoms - i - 2)*(natoms - i - 3))/6;
        if(rem < count_i) break;
        rem -= count_i;
    }
    for(j = i+1; j < natoms - 2; j++){
        long long count_j = ((long long)(natoms - j - 1)*(natoms - j - 2))/2;
        if(rem < count_j) break;
        rem -= count_j;
    }
    for(k = j+1; k < natoms - 1; k++){
        long long count_k = natoms - k - 1;
        if(rem < count_k) break;
        rem -= count_k;
    }
    l = k + 1 + rem;

    // Get all six pair indices
    long long idx_ij = get_pair_index(i, j, natoms);
    long long idx_ik = get_pair_index(i, k, natoms);
    long long idx_il = get_pair_index(i, l, natoms);
    long long idx_jk = get_pair_index(j, k, natoms);
    long long idx_jl = get_pair_index(j, l, natoms);
    long long idx_kl = get_pair_index(k, l, natoms);

    PairData pairs[6] = {
        d_pairs[idx_ij], d_pairs[idx_ik], d_pairs[idx_il],
        d_pairs[idx_jk], d_pairs[idx_jl], d_pairs[idx_kl]
    };

    // Check all pairs against cutoff
    for(int p=0; p<6; p++) {
        if(check_pair(pairs[p], rcout_4b)) return;
    }

    // Compute transformed values
    double trans[6];
    for(int p=0; p<6; p++) {
        trans[p] = compute_trans(pairs[p], rcout_4b);
    }

    // Store results with buffer bounds checking
    int index = atomicAdd(d_count, 1);
    if(index < BUFFER_SIZE_4B) {
        for(int p=0; p<6; p++) {
            d_4b_direct[index*6 + p] = pairs[p].dist * pairs[p].weight;
            d_4b_trans[index*6 + p] = trans[p];
        }
    }
}

void write_results_to_file(const std::vector<double>& h_2b_direct, 
                            const std::vector<double>& h_2b_trans,
                            const std::vector<double>& h_3b_direct, 
                            const std::vector<double>& h_3b_trans,
                            const std::vector<double>& h_4b_direct, 
                            const std::vector<double>& h_4b_trans) {
    // Open output files for 2-body, 3-body, and 4-body results
    std::ofstream out_2b_direct("000.2b_clu-r.txt");
    std::ofstream out_2b_trans("000.2b_clu-s.txt");
    std::ofstream out_3b_direct("000.3b_clu-r.txt");
    std::ofstream out_3b_trans("000.3b_clu-s.txt");
    std::ofstream out_4b_direct("000.4b_clu-r.txt");
    std::ofstream out_4b_trans("000.4b_clu-s.txt");

    // Write the 2-body results
    for (size_t i = 0; i < h_2b_direct.size(); ++i) {
        out_2b_direct << h_2b_direct[i] << std::endl;
    }
    for (size_t i = 0; i < h_2b_trans.size(); ++i) {
        out_2b_trans << h_2b_trans[i] << std::endl;
    }

    // Write the 3-body results
    for (size_t i = 0; i < h_3b_direct.size(); i += 3) {
        out_3b_direct << h_3b_direct[i] << " "
                      << h_3b_direct[i+1] << " "
                      << h_3b_direct[i+2] << std::endl;
    }
    for (size_t i = 0; i < h_3b_trans.size(); i += 3) {
        out_3b_trans << h_3b_trans[i] << " "
                     << h_3b_trans[i+1] << " "
                     << h_3b_trans[i+2] << std::endl;
    }

    // Write the 4-body results
    for (size_t i = 0; i < h_4b_direct.size(); i += 6) {
        out_4b_direct << h_4b_direct[i]   << " "
                      << h_4b_direct[i+1] << " "
                      << h_4b_direct[i+2] << " "
                      << h_4b_direct[i+3] << " "
                      << h_4b_direct[i+4] << " "
                      << h_4b_direct[i+5] << std::endl;
    }
    for (size_t i = 0; i < h_4b_trans.size(); i += 6) {
        out_4b_trans << h_4b_trans[i]   << " "
                     << h_4b_trans[i+1] << " "
                     << h_4b_trans[i+2] << " "
                     << h_4b_trans[i+3] << " "
                     << h_4b_trans[i+4] << " "
                     << h_4b_trans[i+5] << std::endl;
    }

    // Close the output files
    out_2b_direct.close();
    out_2b_trans.close();
    out_3b_direct.close();
    out_3b_trans.close();
    out_4b_direct.close();
    out_4b_trans.close();

    std::cout << "Results written to files." << std::endl;
}

// ---------------------------------------------------------------------
// Host Code: main()
// (Here we assume that the configuration file reading and coordinate file reading
//  have been done on the host. We then convert the atom types (string) to integers,
//  fill the h_coords array, set up box dimensions, and fill the atom_map arrays.)
//---------------------------------------------------------------------
// Main
int main(int argc, char* argv[]) {
    cout << "Beginning computation\n";
    struct timeval tv;
    gettimeofday(&tv, NULL);
    uint64_t start = tv.tv_sec*(uint64_t)1000000 + tv.tv_usec;

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <config_file>" << endl;
        return 1;
    }
    
    // === Read configuration file ===
    unordered_map<string, string> config = readConfig(argv[1]);

    // Extract string values (e.g. the trajectory file path)
    string traj_path = config["TRAJ_PATH"];
    // Remove surrounding quotes if present.
    if (!traj_path.empty() && traj_path.front() == '"' && traj_path.back() == '"') {
        traj_path = traj_path.substr(1, traj_path.size() - 2);
    }
    
    // Extract cutoff values
    double rcout_2b = stod(config["RCOUT_2B"]);
    double rcout_3b = stod(config["RCOUT_3B"]);
    double rcout_4b = stod(config["RCOUT_4B"]);
    int n_frames = stoi(config["N_FRAMES"]);
    cout << "RCOUT_2B: " << rcout_2b << ", RCOUT_3B: " << rcout_3b
     << ", RCOUT_4B: " << rcout_4b << endl;
    
    // Extract system-specific parameters
    vector<string> atom_types = parseStringList(config["ATOM_TYPES"]);
    int ntypes = atom_types.size();
    vector<double> property_set = parseDoubleList(config["PROPERTY_SET"]);
    vector<double> rcin_list = parseDoubleList(config["RCIN_LIST"]);
    vector<double> lambda_set = parseDoubleList(config["LAMBDA_SET"]);
    double alpha = stod(config["ALPHA"]);
    
    // Debug output:
    cout << "Trajectory Path: " << traj_path << endl;
    cout << "RCOUT_2B: " << rcout_2b << ", RCOUT_3B: " << rcout_3b
         << ", RCOUT_4B: " << rcout_4b << endl;
    cout << "ATOM_TYPES: ";
    for (const auto &at : atom_types) cout << at << " ";
    cout << "\nPROPERTY_SET: ";
    for (double p : property_set) cout << p << " ";
    cout << "\nRCIN_LIST: ";
    for (double r : rcin_list) cout << r << " ";
    cout << "\nLAMBDA_SET: ";
    for (double l : lambda_set) cout << l << " ";
    cout << "\nALPHA: " << alpha << endl;

    // Read in data file
    ifstream coordstream;
    coordstream.open(traj_path);
    if (!coordstream.good())
    {
        std::cout << "ERROR: Cannot open xyz file " << traj_path << endl;
        exit(0);
    }
    
    // === Set up the system ===
    // For demonstration, we assume a fixed number of atoms.
    string  line;
    vector<string> line_contents;
    int     ncontents;

    get_next_line(coordstream, line);
    ncontents = split_line(line, line_contents);

    int natoms = stoi(line_contents[0]);
    
    // Set up box dimensions (example values)
    xyz     h_box;
    
    get_next_line(coordstream, line);
    ncontents = split_line(line, line_contents);
    
    h_box.x = stod(line_contents[1]);
    h_box.y = stod(line_contents[5]);
    h_box.z = stod(line_contents[9]);
    
    // Allocate and (in a real code) fill the host coordinate array.

    xyz* h_coords = new xyz[natoms];
    for (int i = 0; i < natoms; i++) {
        get_next_line(coordstream, line);
        ncontents = split_line(line, line_contents);
        h_coords[i].x = stod(line_contents[1]);  // Replace with real coordinate data
        h_coords[i].y = stod(line_contents[2]);
        h_coords[i].z = stod(line_contents[3]);
        // For example, assign atom types: first half type 0, second half type 1.
        h_coords[i].atom_type = 0;
    }

    // === Build host arrays for the atom map parameters ===
    // We need arrays of size ntypes*ntypes. For two atom types (e.g. "C" and "N")
    // the unique pairs (with i<=j) are:
    // (0,0): CC, (0,1): CN, (1,1): NN.
    // We assume rcin_list and lambda_set are provided for these three entries.
    int pairCount = ntypes * ntypes;
    vector<double> h_rcin(pairCount, 0.0);
    vector<float> h_lambda(pairCount, 0.0f);
    vector<float> h_weight(pairCount, 0.0f);
    
    // Fill rcin and lambda.
    // We assume the order in the config lists is the order of unique pairs when looping i=0..ntypes-1 and j=i..ntypes-1.
    int idx = 0;
    for (int i = 0; i < ntypes; i++) {
        for (int j = i; j < ntypes; j++) {
            double rcin_val = rcin_list[idx];
            double lambda_val = lambda_set[idx];
            h_rcin[i*ntypes + j] = rcin_val;
            h_rcin[j*ntypes + i] = rcin_val;  // symmetric assignment
            h_lambda[i*ntypes + j] = static_cast<float>(lambda_val);
            h_lambda[j*ntypes + i] = static_cast<float>(lambda_val);
            idx++;
        }
    }
    
    // Compute weights based on the property_set.
    // For each unique pair, compute the combined property: sqrt(prop_i * prop_j)
    vector<double> prop_combined;
    for (int i = 0; i < ntypes; i++) {
        for (int j = i; j < ntypes; j++) {
            double prop_val = sqrt(property_set[i] * property_set[j]);
            prop_combined.push_back(prop_val);
        }
    }
    double min_prop = *min_element(prop_combined.begin(), prop_combined.end());
    double max_prop = *max_element(prop_combined.begin(), prop_combined.end());
    idx = 0;
    for (int i = 0; i < ntypes; i++) {
        for (int j = i; j < ntypes; j++) {
            double prop_val = sqrt(property_set[i] * property_set[j]);
            float weight = static_cast<float>(1);
            h_weight[i*ntypes + j] = weight;
            h_weight[j*ntypes + i] = weight;
            idx++;
        }
    }
    
    // === Copy data to device ===
    xyz* d_coords;
    cudaMalloc((void**)&d_coords, natoms * sizeof(xyz));
    cudaMemcpy(d_coords, h_coords, natoms * sizeof(xyz), cudaMemcpyHostToDevice);
    
    double* d_rcin;
    cudaMalloc((void**)&d_rcin, pairCount * sizeof(double));
    cudaMemcpy(d_rcin, h_rcin.data(), pairCount * sizeof(double), cudaMemcpyHostToDevice);
    
    float* d_lambda;
    cudaMalloc((void**)&d_lambda, pairCount * sizeof(float));
    cudaMemcpy(d_lambda, h_lambda.data(), pairCount * sizeof(float), cudaMemcpyHostToDevice);
    
    float* d_weight;
    cudaMalloc((void**)&d_weight, pairCount * sizeof(float));
    cudaMemcpy(d_weight, h_weight.data(), pairCount * sizeof(float), cudaMemcpyHostToDevice);
    
    // Preallocate maximum output space for clusters.
    long long max_2b = ((long long)natoms*(natoms-1))/2;
    long long max_3b = ((long long)natoms*(natoms-1)*(natoms-2))/6;
    long long max_4b = ((long long)natoms*(natoms-1)*(natoms-2)*(natoms-3))/24;
    printf("Max 2b: %lld, Max 3b: %lld, Max 4b: %lld\n", max_2b, max_3b, max_4b);
    // int max_2b = 10000;
    // int max_3b = 10000;
    // int max_4b = 10000;
    
    double* d_2b_direct; cudaMalloc((void**)&d_2b_direct, max_2b * sizeof(double));
    double* d_2b_trans;  cudaMalloc((void**)&d_2b_trans,  max_2b * sizeof(double));
    int* d_count_2b; cudaMalloc((void**)&d_count_2b, sizeof(int));
    cudaMemset(d_count_2b, 0, sizeof(int));
    
   
    
    // === Launch the kernels ===
    // (Assuming your CUDA kernels kernel_2b, kernel_3b, and kernel_4b are defined elsewhere.)
    // For 2-body clusters, use a 2D grid:
    // ==== AFTER COPYING DATA TO DEVICE ====

    // --- Precompute pairs ---
    long long totalPairs = ((long long)natoms * (natoms - 1)) / 2;
    PairData* d_pairs;
    cudaMalloc(&d_pairs, totalPairs * sizeof(PairData));

    int threadsPairs = 256;
    long long blocksPairs = (totalPairs + threadsPairs - 1) / threadsPairs;
    kernel_pairs<<<blocksPairs, threadsPairs>>>(d_coords, natoms, h_box, d_rcin, d_lambda, d_weight, d_pairs, ntypes);

    cudaDeviceSynchronize();

    // --- Launch 2-body kernel (unchanged) ---
    dim3 block2(16, 16);
    dim3 grid2((natoms + block2.x - 1)/block2.x, (natoms + block2.y - 1)/block2.y);
    kernel_2b<<<grid2, block2>>>(d_coords, natoms, h_box, rcout_2b, d_rcin, d_lambda, d_weight,
                                 d_2b_direct, d_2b_trans, d_count_2b, ntypes);
    const long long BUFFER_SIZE_3B = 1000000;  // 1 million clusters per chunk
    const long long totalTriplets = ((long long)natoms*(natoms-1)*(natoms-2))/6;
    
    // Allocate device buffers
    double* d_3b_direct;
    double* d_3b_trans;
    int* d_count_3b;
    cudaMalloc(&d_3b_direct, BUFFER_SIZE_3B * 3 * sizeof(double));
    cudaMalloc(&d_3b_trans, BUFFER_SIZE_3B * 3 * sizeof(double));
    cudaMalloc(&d_count_3b, sizeof(int));
    
    // Open output files
    std::ofstream out_3b_direct("000.3b_clu-r.txt");
    std::ofstream out_3b_trans("000.3b_clu-s.txt");
    
    // Process in chunks
    for(long long start=0; start<totalTriplets; start+=BUFFER_SIZE_3B) {
        long long end = std::min(start + BUFFER_SIZE_3B, totalTriplets);
        long long chunk_size = end - start;
        
        // Reset counter and launch kernel
        cudaMemset(d_count_3b, 0, sizeof(int));
        int threads = 256;
        int blocks = (chunk_size + threads - 1) / threads;
        
        kernel_3b<<<blocks, threads>>>(
            d_pairs, natoms, rcout_3b,
            d_3b_direct, d_3b_trans, d_count_3b, ntypes,
            start, end
        );
        
        // Copy results
        int count;
        cudaMemcpy(&count, d_count_3b, sizeof(int), cudaMemcpyDeviceToHost);
        count = std::min(count, static_cast<int>(BUFFER_SIZE_3B));
    
        if(count > 0) {
            std::vector<double> h_direct(count * 3);
            std::vector<double> h_trans(count * 3);
            cudaMemcpy(h_direct.data(), d_3b_direct, count*3*sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(h_trans.data(), d_3b_trans, count*3*sizeof(double), cudaMemcpyDeviceToHost);
    
            // Write to files
            for(int i=0; i<count; i++) {
                out_3b_direct << h_direct[i*3] << " " 
                            << h_direct[i*3+1] << " " 
                            << h_direct[i*3+2] << "\n";
                out_3b_trans << h_trans[i*3] << " " 
                            << h_trans[i*3+1] << " " 
                            << h_trans[i*3+2] << "\n";
            }
        }
    }
    
    // Cleanup
    cudaFree(d_3b_direct);
    cudaFree(d_3b_trans);
    cudaFree(d_count_3b);
    
    // --------------------------------------------------
    // 4-Body Cluster Processing
    // --------------------------------------------------
    const long long BUFFER_SIZE_4B = 1000000;  // 1 million clusters per chunk
    const long long totalQuadruplets = ((long long)natoms*(natoms-1)*(natoms-2)*(natoms-3))/24;
    // const long long totalQuadruplets = 5000000;
    
    // Allocate device buffers
    double* d_4b_direct;
    double* d_4b_trans;
    int* d_count_4b;
    cudaMalloc(&d_4b_direct, BUFFER_SIZE_4B * 6 * sizeof(double));
    cudaMalloc(&d_4b_trans, BUFFER_SIZE_4B * 6 * sizeof(double));
    cudaMalloc(&d_count_4b, sizeof(int));
    
    // Open output files
    std::ofstream out_4b_direct("000.4b_clu-r.txt");
    std::ofstream out_4b_trans("000.4b_clu-s.txt");
    
    // Process in chunks
    for(long long start=0; start<totalQuadruplets; start+=BUFFER_SIZE_4B) {
        long long end = std::min(start + BUFFER_SIZE_4B, totalQuadruplets);
        long long chunk_size = end - start;
        
        // Reset counter and launch kernel
        cudaMemset(d_count_4b, 0, sizeof(int));
        int threads = 256;
        int blocks = (chunk_size + threads - 1) / threads;
        
        kernel_4b<<<blocks, threads>>>(
            d_pairs, natoms, rcout_4b,
            d_4b_direct, d_4b_trans, d_count_4b, ntypes,
            start, end
        );
        
        // Copy results
        int count;
        cudaMemcpy(&count, d_count_4b, sizeof(int), cudaMemcpyDeviceToHost);
        count = std::min(count, static_cast<int>(BUFFER_SIZE_4B));
    
        if(count > 0) {
            std::vector<double> h_direct(count * 6);
            std::vector<double> h_trans(count * 6);
            cudaMemcpy(h_direct.data(), d_4b_direct, count*6*sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(h_trans.data(), d_4b_trans, count*6*sizeof(double), cudaMemcpyDeviceToHost);
    
            // Write to files
            for(int i=0; i<count; i++) {
                out_4b_direct << h_direct[i*6] << " " 
                            << h_direct[i*6+1] << " " 
                            << h_direct[i*6+2] << " " 
                            << h_direct[i*6+3] << " " 
                            << h_direct[i*6+4] << " " 
                            << h_direct[i*6+5] << "\n";
                out_4b_trans << h_trans[i*6] << " " 
                            << h_trans[i*6+1] << " " 
                            << h_trans[i*6+2] << " " 
                            << h_trans[i*6+3] << " " 
                            << h_trans[i*6+4] << " " 
                            << h_trans[i*6+5] << "\n";
            }
        }
    }
    
    // Cleanup
    cudaFree(d_4b_direct);
    cudaFree(d_4b_trans);
    cudaFree(d_count_4b);
    cudaError_t err;


    

    // (You would now write the output to disk as needed.)
    int h_count_2b;
cudaMemcpy(&h_count_2b, d_count_2b, sizeof(int), cudaMemcpyDeviceToHost);

vector<double> h_2b_direct(h_count_2b);
vector<double> h_2b_trans(h_count_2b);
cudaMemcpy(h_2b_direct.data(), d_2b_direct, h_count_2b * sizeof(double), cudaMemcpyDeviceToHost);
cudaMemcpy(h_2b_trans.data(), d_2b_trans, h_count_2b * sizeof(double), cudaMemcpyDeviceToHost);

// === Write 2-body results ===
std::ofstream out_2b_direct("000.2b_clu-r.txt");
std::ofstream out_2b_trans("000.2b_clu-s.txt");
for(int i=0; i<h_count_2b; i++) {
    out_2b_direct << h_2b_direct[i] << endl;
    out_2b_trans << h_2b_trans[i] << endl;
}

    
    // === Free device and host memory ===
    cudaFree(d_pairs);  // Add this line
    cudaFree(d_coords);
    cudaFree(d_coords);
    cudaFree(d_rcin);
    cudaFree(d_lambda);
    cudaFree(d_weight);
    cudaFree(d_2b_direct);
    cudaFree(d_2b_trans);
    cudaFree(d_count_2b);
    cudaFree(d_3b_direct);
    cudaFree(d_3b_trans);
    cudaFree(d_count_3b);
    cudaFree(d_4b_direct);
    cudaFree(d_4b_trans);
    cudaFree(d_count_4b);
    
    delete[] h_coords;
    gettimeofday(&tv, NULL);
    uint64_t end = tv.tv_sec*(uint64_t)1000000 + tv.tv_usec;
    uint64_t elapsed = end - start;
    cout << "@@@ Elapsed time (usec): " << elapsed << "\n";
    cout << "Processing complete. Preparing output.\n";
    return 0;
}