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

using namespace std;


struct xyz
{
    double x;
    double y;
    double z;
    string atom_type; // Added for multi-element
};


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


double get_dist(xyz box, xyz a1, xyz a2, map<string, pair<float, double>> atom_map, float max_lambda)
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

        // cout << "Lambda for " << atom_pair << ": " << lambda << endl;
        // cout << "Rcin for " << atom_pair << ": " << rcin << endl;

    } else {
        cout << "Atom pair " << atom_pair << " not found in the map." << endl;
        return -1; // Return an error if the atom pair is not found
    }
    
    float weight = lambda / max_lambda; // Calculate weight
    
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

        // cout << "Lambda for " << atom_pair << ": " << lambda << endl;
        // cout << "Rcin for " << atom_pair << ": " << rcin << endl;

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

int main()
{
    
    /////////////////////////////////////////////
    // Hard-coded for now, for a single atom type: rcutin, rcut out, morse lambda
    /////////////////////////////////////////////
    
    double rcout_2b = 5.0;
    double rcout_3b = 5.0;
    double rcout_4b = 4.5;

    vector<string> atom_pairs = {"CC", "NN", "CN"};
    vector<float> lambda_set = {1.4, 1.09, 1.34};
    vector<double> rcin_list = {0.98, 0.86, 0.90};

    // Create the dictionary (map) with std::pair<float, double> as the value type
    map<string, pair<float, double>> atom_map;
    for (size_t i = 0; i < atom_pairs.size(); ++i) {
        atom_map[atom_pairs[i]] = make_pair(lambda_set[i], rcin_list[i]);
    }

    // Get the maximum lambda value
    double min_rcin = *min_element(rcin_list.begin(), rcin_list.end());
    double max_lambda = *max_element(lambda_set.begin(), lambda_set.end());
    cout << "Max lambda: " << max_lambda << std::endl;

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
//     for (int i=0; i<cludists_2b.size(); i++)
//         ofstream_2b_s << transform(rcin, rcout_2b, lambda, cludists_2b[i][0]) << endl;

    
//     for (int i=0; i<cludists_3b.size(); i++)
//     {
//         for (int j=0; j<3; j++)
//             ofstream_3b_s << transform(rcin, rcout_3b, lambda, cludists_3b[i][j]) << " ";
//         ofstream_3b_s << endl;
//     }
//     for (int i=0; i<cludists_4b.size(); i++)
//     {
//         for (int j=0; j<6; j++)
//             ofstream_4b_s << transform(rcin, rcout_4b, lambda, cludists_4b[i][j]) <<  " ";  
//         ofstream_4b_s << endl;  
//     }    
    ofstream_2b_s.close();
    ofstream_3b_s.close();
    ofstream_4b_s.close();
     
    
    
}
