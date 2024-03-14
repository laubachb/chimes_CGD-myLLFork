#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<cmath>
#include<algorithm> // For sort

using namespace std;

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

int main(int argc, char *argv[])
{
    
    string f1 = argv[1];
    string f2 = argv[2];
    
    ifstream f1stream;
    ifstream f2stream;
    
    f1stream.open(f1); 
    f2stream.open(f2); 

    string                  line1, line2;
    vector<string>          line1_contents, line2_contents;
    int                     items1, items2;
    
    double                  similiarity = 0;


    while ( get_next_line(f1stream, line1) && get_next_line(f2stream, line2))
    {
        items1 = split_line(line1, line1_contents);
        items2 = split_line(line2, line2_contents);
        
        if (line1_contents[0] != line2_contents[0])
        {
            cout << "ERROR: histogram x-value mismatch!" << endl;
            exit(0);
        }
        
        similiarity += abs(stod(line1_contents[1])-stod(line2_contents[1]));
        
    }
    
    cout << similiarity << endl;
}

    