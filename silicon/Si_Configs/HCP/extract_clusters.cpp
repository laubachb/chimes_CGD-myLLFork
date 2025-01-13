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

using namespace std;


struct xyz
{
    double x;
    double y;
    double z;
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


double get_dist(xyz box, xyz a1, xyz a2)
{
    double dx = a1.x - a2.x; dx -= box.x*round(dx/box.x);
    double dy = a1.y - a2.y; dy -= box.y*round(dy/box.y);
    double dz = a1.z - a2.z; dz -= box.z*round(dz/box.z);
    
    
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    
    if (dist < 0.98)
    {
        cout << "Something went wrong between atoms: " << endl;
        cout << "A: " << a1.x << " " << a1.y << " " << a1.z << endl;
        cout << "B: " << a2.x << " " << a2.y << " " << a2.z << endl;
        cout << dx << " " << dy << " " << dz << endl;
        exit(0);
    }
        
        
    return sqrt(dx*dx + dy*dy + dz*dz);
}


double transform(double rcin, double rcout,double lambda, double rij)
{
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
    
    double rcin     = 1.55;
    double rcout_2b = 6.0;
    double rcout_3b = 6.0;
    double rcout_4b = 5.5;
    double lambda   = 2.33;

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
        
        coords[i].x = stod(line_contents[1]);
        coords[i].y = stod(line_contents[2]);
        coords[i].z = stod(line_contents[3]);

    }
    
    
    /////////////////////////////////////////////
    // Extract all 2-, 3-, and 4-body clusters
    /////////////////////////////////////////////
    
    vector<vector<double > > cludists_2b;   // dim 1: cluster index; dim 2, distance between atom pair
    vector<vector<double > > cludists_3b;   // dim 1: cluster index; dim 2, distance between atom pairs
    vector<vector<double > > cludists_4b;   // dim 1: cluster index; dim 2, distance between atom pairs
    
    // ... define get dist
    
    vector<double> dist_2b(1);
    vector<double> dist_3b(3);
    vector<double> dist_4b(6);
 
    for (int i=0; i<natoms; i++)
    {
        for (int j=i+1; j<natoms; j++)
        {

            //cout << "Calculating distances between atoms with index and coordinates:" << endl;
            //cout << i << " " << coords[i].x << " " << coords[i].y << " " << coords[i].z << endl;
            //cout << j << " " << coords[j].x << " " << coords[j].y << " " << coords[j].z << endl;
            
            dist_2b[0] = get_dist(boxdims, coords[i], coords[j]);
            
            //cout << dist_2b[0] << endl;
            
            if (dist_2b[0] < rcout_2b)
            {
                //cout << dist_2b[0] << endl;
                
                cludists_2b.push_back(dist_2b);
                
                if (dist_2b[0] >= rcout_3b)
                    continue;
                
                for (int k=j+1; k<natoms; k++)
                {
                    dist_3b[0] = dist_2b[0]; 
                    dist_3b[1] = get_dist(boxdims, coords[i], coords[k]);
                    
                    if (dist_3b[1] >= rcout_3b)
                        continue;
                    
                    dist_3b[2] = get_dist(boxdims, coords[j], coords[k]);
                    
                    if (dist_3b[2] >= rcout_3b)
                        continue;
                    
                    cludists_3b.push_back(dist_3b);
                    
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
                        
                        dist_4b[3] = get_dist(boxdims, coords[i], coords[l]);
                    
                        if (dist_4b[3] >= rcout_4b)
                            continue;       

                        dist_4b[4] = get_dist(boxdims, coords[j], coords[l]);
                    
                        if (dist_4b[4] >= rcout_4b)
                            continue;      
                        
                        dist_4b[5] = get_dist(boxdims, coords[k], coords[l]);
                    
                        if (dist_4b[5] >= rcout_4b)
                            continue;  
                        
                        cludists_4b.push_back(dist_4b);                                                                     
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
    
    for (int i=0; i<cludists_2b.size(); i++)
        ofstream_2b_r << cludists_2b[i][0] << endl;
    
    for (int i=0; i<cludists_3b.size(); i++)
    {
        for (int j=0; j<3; j++)
            ofstream_3b_r << cludists_3b[i][j] <<  " ";
        ofstream_3b_r << endl;
    }
    
    for (int i=0; i<cludists_4b.size(); i++)
    {
        for (int j=0; j<6; j++)
            ofstream_4b_r << cludists_4b[i][j] <<  " "; 
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
    
    for (int i=0; i<cludists_2b.size(); i++)
        ofstream_2b_s << transform(rcin, rcout_2b, lambda, cludists_2b[i][0]) << endl;

    
    for (int i=0; i<cludists_3b.size(); i++)
    {
        for (int j=0; j<3; j++)
            ofstream_3b_s << transform(rcin, rcout_3b, lambda, cludists_3b[i][j]) << " ";
        ofstream_3b_s << endl;
    }
    for (int i=0; i<cludists_4b.size(); i++)
    {
        for (int j=0; j<6; j++)
            ofstream_4b_s << transform(rcin, rcout_4b, lambda, cludists_4b[i][j]) <<  " ";  
        ofstream_4b_s << endl;  
    }    
    ofstream_2b_s.close();
    ofstream_3b_s.close();
    ofstream_4b_s.close();
     
    
    
}
