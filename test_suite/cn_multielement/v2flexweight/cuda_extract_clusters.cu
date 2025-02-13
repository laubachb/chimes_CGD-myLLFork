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

// CUDA runtime
#include <cuda_runtime.h>

using namespace std;

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
    float weight = d_weight[a * ntypes + b];
    double rcin = d_rcin[a * ntypes + b];
    
    double dx = d_coords[i].x - d_coords[j].x;
    double dy = d_coords[i].y - d_coords[j].y;
    double dz = d_coords[i].z - d_coords[j].z;
    printf("Coords: %f, %f, %f\n", d_coords[i].x, d_coords[i].y, d_coords[i].z);
    dx = periodic_diff(dx, box.x);
    dy = periodic_diff(dy, box.y);
    dz = periodic_diff(dz, box.z);
    printf("Distances: %f, %f, %f\n", dx, dy, dz);
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    printf("distance between atoms (rcin = %f): %f, rcout_2b: %f\n", rcin, dist, rcout_2b);
    // if (dist < rcout_2b) {
    //     printf("Cluster formed!\n");
    // }
    if(dist < rcin) { 
        dist = 1e8;  // flag error condition 
    }

    double weighted_dist = dist * weight;
    
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
    }
}

// ---------------------------------------------------------------------
// Kernel for 3-body clusters (triplets)
// We flatten the triplet (i,j,k) with i < j < k into a single linear index.
// The unranking routine converts a thread’s linear index to (i,j,k).
__global__ void kernel_3b(const xyz* d_coords, int natoms, xyz box, double rcout_3b, 
                           const double* d_rcin, const float* d_lambda, const float* d_weight,
                           double* d_3b_direct, double* d_3b_trans, int* d_count, int ntypes) {
    // total number of triplets = n*(n-1)*(n-2)/6
    long long totalTriplets = ((long long)natoms*(natoms-1)*(natoms-2))/6;
    long long idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= totalTriplets) return;
    
    // Unrank idx into unique (i,j,k)
    int i, j, k;
    long long rem = idx;
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
    
    // For the triplet (i,j,k), first check the (i,j) pair
    int type_i = d_coords[i].atom_type;
    int type_j = d_coords[j].atom_type;
    int type_k = d_coords[k].atom_type;
    
    // (i,j)
    int a = (type_i < type_j) ? type_i : type_j;
    int b = (type_i < type_j) ? type_j : type_i;
    float weight_ij = d_weight[a * ntypes + b];
    double rcin_ij = d_rcin[a * ntypes + b];
    double dx = d_coords[i].x - d_coords[j].x;
    double dy = d_coords[i].y - d_coords[j].y;
    double dz = d_coords[i].z - d_coords[j].z;
    dx = periodic_diff(dx, box.x);
    dy = periodic_diff(dy, box.y);
    dz = periodic_diff(dz, box.z);
    double dist_ij = sqrt(dx*dx + dy*dy + dz*dz);
    if(dist_ij < rcin_ij) dist_ij = 1e8;
    double weighted_ij = dist_ij * weight_ij;
    if(dist_ij >= rcout_3b) return;  // skip if too long

    // (i,k)
    a = (type_i < type_k) ? type_i : type_k;
    b = (type_i < type_k) ? type_k : type_i;
    float weight_ik = d_weight[a * ntypes + b];
    double rcin_ik = d_rcin[a * ntypes + b];
    dx = d_coords[i].x - d_coords[k].x;
    dy = d_coords[i].y - d_coords[k].y;
    dz = d_coords[i].z - d_coords[k].z;
    dx = periodic_diff(dx, box.x);
    dy = periodic_diff(dy, box.y);
    dz = periodic_diff(dz, box.z);
    double dist_ik = sqrt(dx*dx + dy*dy + dz*dz);
    if(dist_ik < rcin_ik) dist_ik = 1e8;
    double weighted_ik = dist_ik * weight_ik;
    if(dist_ik >= rcout_3b) return;
    
    // (j,k)
    a = (type_j < type_k) ? type_j : type_k;
    b = (type_j < type_k) ? type_k : type_j;
    float weight_jk = d_weight[a * ntypes + b];
    double rcin_jk = d_rcin[a * ntypes + b];
    dx = d_coords[j].x - d_coords[k].x;
    dy = d_coords[j].y - d_coords[k].y;
    dz = d_coords[j].z - d_coords[k].z;
    dx = periodic_diff(dx, box.x);
    dy = periodic_diff(dy, box.y);
    dz = periodic_diff(dz, box.z);
    double dist_jk = sqrt(dx*dx + dy*dy + dz*dz);
    if(dist_jk < rcin_jk) dist_jk = 1e8;
    double weighted_jk = dist_jk * weight_jk;
    if(dist_jk >= rcout_3b) return;
    
    // Compute transformed values for each pair:
    double lambda_ij = d_lambda[((type_i < type_j) ? type_i : type_j)*ntypes + ((type_i < type_j) ? type_j : type_i)];
    double x_min_val = exp(-rcin_ij / lambda_ij);
    double x_max_val = exp(-rcout_3b / lambda_ij);
    double x_avg = 0.5*(x_min_val + x_max_val);
    double x_diff = -0.5*(x_max_val - x_min_val);
    double trans_ij = (exp(-weighted_ij / lambda_ij) - x_avg) / x_diff;
    
    double lambda_ik = d_lambda[((type_i < type_k) ? type_i : type_k)*ntypes + ((type_i < type_k) ? type_k : type_i)];
    x_min_val = exp(-rcin_ik / lambda_ik);
    x_max_val = exp(-rcout_3b / lambda_ik);
    x_avg = 0.5*(x_min_val + x_max_val);
    x_diff = -0.5*(x_max_val - x_min_val);
    double trans_ik = (exp(-weighted_ik / lambda_ik) - x_avg) / x_diff;
    
    double lambda_jk = d_lambda[((type_j < type_k) ? type_j : type_k)*ntypes + ((type_j < type_k) ? type_k : type_j)];
    x_min_val = exp(-rcin_jk / lambda_jk);
    x_max_val = exp(-rcout_3b / lambda_jk);
    x_avg = 0.5*(x_min_val + x_max_val);
    x_diff = -0.5*(x_max_val - x_min_val);
    double trans_jk = (exp(-weighted_jk / lambda_jk) - x_avg) / x_diff;
    
    // Store the triplet (each cluster has three distances and three transformed values)
    int index = atomicAdd(d_count, 1);
    d_3b_direct[index*3 + 0] = weighted_ij;
    d_3b_direct[index*3 + 1] = weighted_ik;
    d_3b_direct[index*3 + 2] = weighted_jk;
    d_3b_trans[index*3 + 0] = trans_ij;
    d_3b_trans[index*3 + 1] = trans_ik;
    d_3b_trans[index*3 + 2] = trans_jk;
}

// ---------------------------------------------------------------------
// Kernel for 4-body clusters (quadruplets)
// Here we “unrank” a linear index into a unique (i,j,k,l) combination (with i < j < k < l).
__global__ void kernel_4b(const xyz* d_coords, int natoms, xyz box, double rcout_4b, 
                           const double* d_rcin, const float* d_lambda, const float* d_weight,
                           double* d_4b_direct, double* d_4b_trans, int* d_count, int ntypes) {
    long long totalQuad = ((long long)natoms*(natoms-1)*(natoms-2)*(natoms-3))/24;
    long long idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= totalQuad) return;
    int i, j, k, l;
    long long rem = idx;
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
    
    int type_i = d_coords[i].atom_type;
    int type_j = d_coords[j].atom_type;
    int type_k = d_coords[k].atom_type;
    int type_l = d_coords[l].atom_type;
    
    // Compute all six pair distances and check against rcout_4b.
    // (i,j)
    int a = (type_i < type_j) ? type_i : type_j;
    int b = (type_i < type_j) ? type_j : type_i;
    float weight_ij = d_weight[a * ntypes + b];
    double rcin_ij = d_rcin[a * ntypes + b];
    double dx = d_coords[i].x - d_coords[j].x;
    double dy = d_coords[i].y - d_coords[j].y;
    double dz = d_coords[i].z - d_coords[j].z;
    dx = periodic_diff(dx, box.x);
    dy = periodic_diff(dy, box.y);
    dz = periodic_diff(dz, box.z);
    double dist_ij = sqrt(dx*dx + dy*dy + dz*dz);
    if(dist_ij < rcin_ij) dist_ij = 1e8;
    double weighted_ij = dist_ij * weight_ij;
    if(dist_ij >= rcout_4b) return;
    
    // (i,k)
    a = (type_i < type_k) ? type_i : type_k;
    b = (type_i < type_k) ? type_k : type_i;
    float weight_ik = d_weight[a * ntypes + b];
    double rcin_ik = d_rcin[a * ntypes + b];
    dx = d_coords[i].x - d_coords[k].x;
    dy = d_coords[i].y - d_coords[k].y;
    dz = d_coords[i].z - d_coords[k].z;
    dx = periodic_diff(dx, box.x);
    dy = periodic_diff(dy, box.y);
    dz = periodic_diff(dz, box.z);
    double dist_ik = sqrt(dx*dx + dy*dy + dz*dz);
    if(dist_ik < rcin_ik) dist_ik = 1e8;
    double weighted_ik = dist_ik * weight_ik;
    if(dist_ik >= rcout_4b) return;
    
    // (j,k)
    a = (type_j < type_k) ? type_j : type_k;
    b = (type_j < type_k) ? type_k : type_j;
    float weight_jk = d_weight[a * ntypes + b];
    double rcin_jk = d_rcin[a * ntypes + b];
    dx = d_coords[j].x - d_coords[k].x;
    dy = d_coords[j].y - d_coords[k].y;
    dz = d_coords[j].z - d_coords[k].z;
    dx = periodic_diff(dx, box.x);
    dy = periodic_diff(dy, box.y);
    dz = periodic_diff(dz, box.z);
    double dist_jk = sqrt(dx*dx + dy*dy + dz*dz);
    if(dist_jk < rcin_jk) dist_jk = 1e8;
    double weighted_jk = dist_jk * weight_jk;
    if(dist_jk >= rcout_4b) return;
    
    // (i,l)
    a = (type_i < type_l) ? type_i : type_l;
    b = (type_i < type_l) ? type_l : type_i;
    float weight_il = d_weight[a * ntypes + b];
    double rcin_il = d_rcin[a * ntypes + b];
    dx = d_coords[i].x - d_coords[l].x;
    dy = d_coords[i].y - d_coords[l].y;
    dz = d_coords[i].z - d_coords[l].z;
    dx = periodic_diff(dx, box.x);
    dy = periodic_diff(dy, box.y);
    dz = periodic_diff(dz, box.z);
    double dist_il = sqrt(dx*dx + dy*dy + dz*dz);
    if(dist_il < rcin_il) dist_il = 1e8;
    double weighted_il = dist_il * weight_il;
    if(dist_il >= rcout_4b) return;
    
    // (j,l)
    a = (type_j < type_l) ? type_j : type_l;
    b = (type_j < type_l) ? type_l : type_j;
    float weight_jl = d_weight[a * ntypes + b];
    double rcin_jl = d_rcin[a * ntypes + b];
    dx = d_coords[j].x - d_coords[l].x;
    dy = d_coords[j].y - d_coords[l].y;
    dz = d_coords[j].z - d_coords[l].z;
    dx = periodic_diff(dx, box.x);
    dy = periodic_diff(dy, box.y);
    dz = periodic_diff(dz, box.z);
    double dist_jl = sqrt(dx*dx + dy*dy + dz*dz);
    if(dist_jl < rcin_jl) dist_jl = 1e8;
    double weighted_jl = dist_jl * weight_jl;
    if(dist_jl >= rcout_4b) return;
    
    // (k,l)
    a = (type_k < type_l) ? type_k : type_l;
    b = (type_k < type_l) ? type_l : type_k;
    float weight_kl = d_weight[a * ntypes + b];
    double rcin_kl = d_rcin[a * ntypes + b];
    dx = d_coords[k].x - d_coords[l].x;
    dy = d_coords[k].y - d_coords[l].y;
    dz = d_coords[k].z - d_coords[l].z;
    dx = periodic_diff(dx, box.x);
    dy = periodic_diff(dy, box.y);
    dz = periodic_diff(dz, box.z);
    double dist_kl = sqrt(dx*dx + dy*dy + dz*dz);
    if(dist_kl < rcin_kl) dist_kl = 1e8;
    double weighted_kl = dist_kl * weight_kl;
    if(dist_kl >= rcout_4b) return;
    
    // If all conditions are met, store the six weighted distances and their transformed values.
    int index = atomicAdd(d_count, 1);
    d_4b_direct[index*6 + 0] = weighted_ij;
    d_4b_direct[index*6 + 1] = weighted_ik;
    d_4b_direct[index*6 + 2] = weighted_jk;
    d_4b_direct[index*6 + 3] = weighted_il;
    d_4b_direct[index*6 + 4] = weighted_jl;
    d_4b_direct[index*6 + 5] = weighted_kl;
    
    // For each pair compute the transformed value (similar to above)
    double lambda_ij = d_lambda[((type_i < type_j) ? type_i : type_j)*ntypes + ((type_i < type_j) ? type_j : type_i)];
    double x_min_val = exp(-rcin_ij / lambda_ij);
    double x_max_val = exp(-rcout_4b / lambda_ij);
    double x_avg = 0.5*(x_min_val + x_max_val);
    double x_diff = -0.5*(x_max_val - x_min_val);
    double trans_ij = (exp(-weighted_ij / lambda_ij) - x_avg) / x_diff;
    
    double lambda_ik = d_lambda[((type_i < type_k) ? type_i : type_k)*ntypes + ((type_i < type_k) ? type_k : type_i)];
    x_min_val = exp(-rcin_ik / lambda_ik);
    x_max_val = exp(-rcout_4b / lambda_ik);
    x_avg = 0.5*(x_min_val + x_max_val);
    x_diff = -0.5*(x_max_val - x_min_val);
    double trans_ik = (exp(-weighted_ik / lambda_ik) - x_avg) / x_diff;
    
    double lambda_jk = d_lambda[((type_j < type_k) ? type_j : type_k)*ntypes + ((type_j < type_k) ? type_k : type_j)];
    x_min_val = exp(-rcin_jk / lambda_jk);
    x_max_val = exp(-rcout_4b / lambda_jk);
    x_avg = 0.5*(x_min_val + x_max_val);
    x_diff = -0.5*(x_max_val - x_min_val);
    double trans_jk = (exp(-weighted_jk / lambda_jk) - x_avg) / x_diff;
    
    double lambda_il = d_lambda[((type_i < type_l) ? type_i : type_l)*ntypes + ((type_i < type_l) ? type_l : type_i)];
    x_min_val = exp(-rcin_il / lambda_il);
    x_max_val = exp(-rcout_4b / lambda_il);
    x_avg = 0.5*(x_min_val + x_max_val);
    x_diff = -0.5*(x_max_val - x_min_val);
    double trans_il = (exp(-weighted_il / lambda_il) - x_avg) / x_diff;
    
    double lambda_jl = d_lambda[((type_j < type_l) ? type_j : type_l)*ntypes + ((type_j < type_l) ? type_l : type_j)];
    x_min_val = exp(-rcin_jl / lambda_jl);
    x_max_val = exp(-rcout_4b / lambda_jl);
    x_avg = 0.5*(x_min_val + x_max_val);
    x_diff = -0.5*(x_max_val - x_min_val);
    double trans_jl = (exp(-weighted_jl / lambda_jl) - x_avg) / x_diff;
    
    double lambda_kl = d_lambda[((type_k < type_l) ? type_k : type_l)*ntypes + ((type_k < type_l) ? type_l : type_k)];
    x_min_val = exp(-rcin_kl / lambda_kl);
    x_max_val = exp(-rcout_4b / lambda_kl);
    x_avg = 0.5*(x_min_val + x_max_val);
    x_diff = -0.5*(x_max_val - x_min_val);
    double trans_kl = (exp(-weighted_kl / lambda_kl) - x_avg) / x_diff;
    
    int base = index * 6;
    d_4b_trans[base + 0] = trans_ij;
    d_4b_trans[base + 1] = trans_ik;
    d_4b_trans[base + 2] = trans_jk;
    d_4b_trans[base + 3] = trans_il;
    d_4b_trans[base + 4] = trans_jl;
    d_4b_trans[base + 5] = trans_kl;
}

void write_results_to_file(const std::vector<double>& h_2b_direct, 
                            const std::vector<double>& h_2b_trans,
                            const std::vector<double>& h_3b_direct, 
                            const std::vector<double>& h_3b_trans,
                            const std::vector<double>& h_4b_direct, 
                            const std::vector<double>& h_4b_trans) {
    // Open output files for 2-body, 3-body, and 4-body results
    std::ofstream out_2b_direct("2b_direct.txt");
    std::ofstream out_2b_trans("2b_trans.txt");
    std::ofstream out_3b_direct("3b_direct.txt");
    std::ofstream out_3b_trans("3b_trans.txt");
    std::ofstream out_4b_direct("4b_direct.txt");
    std::ofstream out_4b_trans("4b_trans.txt");

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
            float weight = static_cast<float>(alpha + (1 - alpha) * ((prop_val - min_prop) / (max_prop - min_prop)));
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
    // int max_2b = (natoms*(natoms-1))/2;
    // int max_3b = (natoms*(natoms-1)*(natoms-2))/6;
    // int max_4b = (natoms*(natoms-1)*(natoms-2)*(natoms-3))/24;
    int max_2b = 10000;
    int max_3b = 10000;
    int max_4b = 10000;
    
    double* d_2b_direct; cudaMalloc((void**)&d_2b_direct, max_2b * sizeof(double));
    double* d_2b_trans;  cudaMalloc((void**)&d_2b_trans,  max_2b * sizeof(double));
    int* d_count_2b; cudaMalloc((void**)&d_count_2b, sizeof(int));
    cudaMemset(d_count_2b, 0, sizeof(int));
    
    double* d_3b_direct; cudaMalloc((void**)&d_3b_direct, max_3b * 3 * sizeof(double));
    double* d_3b_trans;  cudaMalloc((void**)&d_3b_trans,  max_3b * 3 * sizeof(double));
    int* d_count_3b; cudaMalloc((void**)&d_count_3b, sizeof(int));
    cudaMemset(d_count_3b, 0, sizeof(int));
    
    double* d_4b_direct; cudaMalloc((void**)&d_4b_direct, max_4b * 6 * sizeof(double));
    double* d_4b_trans;  cudaMalloc((void**)&d_4b_trans,  max_4b * 6 * sizeof(double));
    int* d_count_4b; cudaMalloc((void**)&d_count_4b, sizeof(int));
    cudaMemset(d_count_4b, 0, sizeof(int));
    
    // === Launch the kernels ===
    // (Assuming your CUDA kernels kernel_2b, kernel_3b, and kernel_4b are defined elsewhere.)
    // For 2-body clusters, use a 2D grid:
    dim3 block2(16, 16);
    dim3 grid2((natoms + block2.x - 1)/block2.x, (natoms + block2.y - 1)/block2.y);
    kernel_2b<<<grid2, block2>>>(d_coords, natoms, h_box, rcout_2b, d_rcin, d_lambda, d_weight,
                                 d_2b_direct, d_2b_trans, d_count_2b, ntypes);
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        cerr << "CUDA kernel launch failed for 2b: " << cudaGetErrorString(error) << endl;
        return 1;
    }

    cudaDeviceSynchronize();
    
    // For 3-body clusters, launch a 1D grid:
    int threads3 = 256;
    int blocks3 = (max_3b + threads3 - 1) / threads3;
    kernel_3b<<<blocks3, threads3>>>(d_coords, natoms, h_box, rcout_3b, d_rcin, d_lambda, d_weight,
                                     d_3b_direct, d_3b_trans, d_count_3b, ntypes);
    cudaError_t error2 = cudaGetLastError();
    if (error != cudaSuccess) {
        cerr << "CUDA kernel launch failed for 3b: " << cudaGetErrorString(error2) << endl;
        return 1;
    }
    cudaDeviceSynchronize();
    
    // For 4-body clusters, launch a 1D grid:
    int threads4 = 256;
    int blocks4 = (max_4b + threads4 - 1) / threads4;
    kernel_4b<<<blocks4, threads4>>>(d_coords, natoms, h_box, rcout_4b, d_rcin, d_lambda, d_weight,
                                     d_4b_direct, d_4b_trans, d_count_4b, ntypes);
    cudaError_t error3 = cudaGetLastError();
    if (error != cudaSuccess) {
        cerr << "CUDA kernel launch failed for 4b: " << cudaGetErrorString(error3) << endl;
        return 1;
    }
    cudaDeviceSynchronize();
    
    // === Copy results back to host ===
    int h_count_2b, h_count_3b, h_count_4b;
    cudaMemcpy(&h_count_2b, d_count_2b, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_count_3b, d_count_3b, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_count_4b, d_count_4b, sizeof(int), cudaMemcpyDeviceToHost);
    
    vector<double> h_2b_direct(h_count_2b);
    vector<double> h_2b_trans(h_count_2b);
    cudaMemcpy(h_2b_direct.data(), d_2b_direct, h_count_2b * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_2b_trans.data(), d_2b_trans, h_count_2b * sizeof(double), cudaMemcpyDeviceToHost);
    
    vector<double> h_3b_direct(h_count_3b * 3);
    vector<double> h_3b_trans(h_count_3b * 3);
    cudaMemcpy(h_3b_direct.data(), d_3b_direct, h_count_3b * 3 * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_3b_trans.data(), d_3b_trans, h_count_3b * 3 * sizeof(double), cudaMemcpyDeviceToHost);
    
    vector<double> h_4b_direct(h_count_4b * 6);
    vector<double> h_4b_trans(h_count_4b * 6);
    cudaMemcpy(h_4b_direct.data(), d_4b_direct, h_count_4b * 6 * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_4b_trans.data(), d_4b_trans, h_count_4b * 6 * sizeof(double), cudaMemcpyDeviceToHost);

    cudaError_t err;

    err = cudaMemcpy(&h_count_2b, d_count_2b, sizeof(int), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        cerr << "CUDA memcpy failed for 2b count: " << cudaGetErrorString(err) << endl;
        return 1;
    }

    err = cudaMemcpy(h_2b_direct.data(), d_2b_direct, h_count_2b * sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        cerr << "CUDA memcpy failed for 2b direct: " << cudaGetErrorString(err) << endl;
        return 1;
    }

    
    cout << "h_count_2b: " << h_count_2b << endl;
    cout << "h_count_3b: " << h_count_3b << endl;
    cout << "h_count_4b: " << h_count_4b << endl;

    // (You would now write the output to disk as needed.)
    write_results_to_file(h_2b_direct, h_2b_trans, h_3b_direct, h_3b_trans, h_4b_direct, h_4b_trans);

    
    // === Free device and host memory ===
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
    
    return 0;
}