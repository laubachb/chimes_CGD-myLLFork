#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstdlib>

using namespace std;

int nprocs;
int my_rank;

//----------------------------------------------------------------
// Utility: splits a line into tokens and removes comments.
int split_line(string line, vector<string> &items)
{
    string contents;
    stringstream sstream;

    // Remove any comments starting with '!' or '##'
    int pos = line.find('!');
    if (pos != string::npos)
        line.erase(pos, line.length() - pos);

    pos = line.find("##");
    if (pos != string::npos)
        line.erase(pos, line.length() - pos);

    pos = line.find('\n');
    if (pos != string::npos)
        line.erase(pos, 1);

    sstream.str(line);
    items.clear();
    while (sstream >> contents)
        items.push_back(contents);

    return items.size();
}

//----------------------------------------------------------------
// Utility: safely reads the next line from a stream.
bool get_next_line(istream &str, string &line)
{
    getline(str, line);
    return bool(str);
}

//----------------------------------------------------------------
// New function to read clusters along with composition data.
// The first npairs_per_cluster tokens are numerical distances,
// and any remaining tokens are concatenated (with spaces)
// to form a composition string.
void read_flat_clusters_with_composition(string clufile, int npairs_per_cluster,
                                         vector<double> &clusters,
                                         vector<string> &compositions)
{
    ifstream clustream(clufile);
    if (!clustream.is_open()) {
        cerr << "Error: Cannot open file " << clufile << endl;
        exit(1);
    }
    
    string line;
    vector<string> line_contents;
    int n_contents;
    vector<double> one_cluster(npairs_per_cluster);

    while (get_next_line(clustream, line))
    {
        n_contents = split_line(line, line_contents);
        if (n_contents < npairs_per_cluster + 1) {
            cout << "ERROR: Expected at least " << npairs_per_cluster + 1
                 << " tokens in " << clufile << ", but got " << n_contents << endl;
            exit(0);
        }
        
        // Read the first npairs_per_cluster tokens as doubles.
        for (int i = 0; i < npairs_per_cluster; i++)
            one_cluster[i] = stod(line_contents[i]);
        
        // Optionally sort the numerical part (as in your original code)
        sort(one_cluster.begin(), one_cluster.end());
        
        // Insert the numerical data into the clusters vector.
        clusters.insert(clusters.end(), one_cluster.begin(), one_cluster.end());
        
        // Join any remaining tokens as the composition (separated by space)
        string comp = "";
        for (int i = npairs_per_cluster; i < n_contents; i++) {
            if (i > npairs_per_cluster) comp += " ";
            comp += line_contents[i];
        }
        compositions.push_back(comp);
    }
    clustream.close();
}

//----------------------------------------------------------------
// Return a histogram bin index based on a distance.
int get_bin(double binw, double maxval, double dist)
{
    int bin = floor(dist/binw);
    if (dist == maxval)
        return bin - 1;
    else
        return bin;
}

//----------------------------------------------------------------
// Divide tasks among MPI processes.
void divide_task(int &my_rank_start, int &my_rank_end, int tasks)
{
    int procs_used;
    if (tasks <= 0) {
        my_rank_start = 1;
        my_rank_end = 0;
        return;
    }
    if (nprocs <= tasks)
        procs_used = nprocs;
    else
        procs_used = tasks;
    
    my_rank_start = ceil((double) my_rank * tasks / procs_used);
    
    if (my_rank > tasks) {
        my_rank_start = tasks + 1;
        my_rank_end = tasks - 1;
    } else if (my_rank == procs_used - 1) {
        my_rank_end = tasks - 1;
    } else {
        my_rank_end = ceil((double)(my_rank+1) * tasks / procs_used) - 1;
        if (my_rank_end > tasks - 1)
            my_rank_end = tasks - 1;
    }
}

//----------------------------------------------------------------
// Helper function: computes a weight for a single cluster based on its composition.
// The more 'N's, the closer to 2.0; the more 'C's, the closer to 0.1.
// Assumes the composition string contains tokens like "C" or "N" (ignores other characters).
double get_cluster_weight(const string &comp)
{
    int countN = 0, countC = 0, total = 0;
    // We'll iterate over each character; if tokens are space separated,
    // this will count letters 'N' and 'C' in the string.
    for (char ch : comp) {
        if (ch == 'N') { countN++; total++; }
        else if (ch == 'C') { countC++; total++; }
    }
    // If none found, default weight (you might choose to handle this differently)
    if (total == 0) return 1.0;
    
    double fracN = double(countN) / total;
    // Linear interpolation:
    // If fracN==1.0 -> weight = 2.0; if fracN==0.0 -> weight = 0.1.
    double weight = 0.1 + fracN * (2.0 - 0.1);
    return weight;
}

//----------------------------------------------------------------
// Computes a combined weight for two clusters based on their compositions.
// Here we average the individual cluster weights.
double composition_weight(const string &comp1, const string &comp2)
{
    double w1 = get_cluster_weight(comp1);
    double w2 = get_cluster_weight(comp2);
    return (w1 + w2) / 2.0;
}

//----------------------------------------------------------------
// Modified histogram generation function.
// It now takes two extra vectors (comp1 and comp2) that hold the composition
// for each cluster. When computing the distance between clusters,
// the Euclidean distance computed from the numerical part is multiplied by a weight
// derived from the cluster compositions.
void gen_flat_hists(vector<double> &clu1, vector<double> &clu2,
                    vector<string> &comp1, vector<string> &comp2,
                    int n_cluster_pairs, int nbin, double binw, double maxd,
                    string histfile, double rcout, bool same = false)
{
    int bin;
    double dist;
    vector<long long int> my_hist(nbin, 0);
    vector<long long int> hist(nbin, 0);
    long long int my_nsamples = 0;
    long long int nsamples = 0;
    int my_rank_start, my_rank_end;
    int looptwo_start;
    int total_tasks;
    int status;
    
    int nclusters1 = clu1.size() / n_cluster_pairs;
    int nclusters2 = clu2.size() / n_cluster_pairs;
    
    divide_task(my_rank_start, my_rank_end, nclusters1);
    total_tasks = my_rank_end - my_rank_start;
    
    if (my_rank == 0)
        cout << "Dividing " << nclusters1 << " clusters across " << nprocs << " processors" << endl;
    
    if (total_tasks > 0)
    {
        for (int i = my_rank_start; i <= my_rank_end; i++)
        {
            status = double(i - my_rank_start) / (total_tasks) * 100.0;
            if (my_rank == 0) {
                if ((total_tasks / 10) == 0)
                    cout << histfile << " Completion percent: " << status << " " << i
                         << " of " << total_tasks << " assigned" << endl;
                else if (i % (total_tasks / 10) == 0)
                    cout << histfile << " Completion percent: " << status << " " << i
                         << " of " << total_tasks << " assigned" << endl;
            }
            
            if (same)
                looptwo_start = i + 1;
            else
                looptwo_start = 0;
            
            for (int j = looptwo_start; j < nclusters2; j++)
            {
                dist = 0;
                for (int k = 0; k < n_cluster_pairs; k++)
                    dist += pow(clu1[i * n_cluster_pairs + k] - clu2[j * n_cluster_pairs + k], 2.0);
                
                double d = sqrt(dist);
                // Weight the distance based on the compositions from the two clusters.
                double weight = composition_weight(comp1[i], comp2[j]);
                double weighted_dist = d * weight;
                
                bin = get_bin(binw, maxd, weighted_dist);
                if (bin >= nbin || bin < 0)
                    continue;
                
                my_hist[bin] += 1;
                my_nsamples += 1;
            }
        }
    }
    
    if (my_rank == 0)
    {
        cout << "Loop done, printing results:" << endl;
        cout << "Counted nsamples: " << my_nsamples << endl;
    }
    
    MPI_Reduce(my_hist.data(), hist.data(), nbin, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&my_nsamples, &nsamples, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (my_rank == 0)
    {
        ofstream cluhist(histfile);
        for (int i = 0; i < hist.size(); i++)
            cluhist << i * binw + 0.5 * binw << " " << double(hist[i]) / nsamples << endl;
        cluhist.close();
        cout << "Printed." << endl;
    }
}

//----------------------------------------------------------------
// Main function.
int main(int argc, char *argv[])
{
    my_rank = 0;
    nprocs  = 1;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == 0)
        cout << "Code compiled in MPI mode.";
    if (my_rank == 0)
        cout << " Will run on " << nprocs << " processor(s)." << endl;

    string f1_idx = argv[1];  // e.g., "0050"
    string f2_idx = argv[2];  // e.g., "0075"
    string style  = argv[3];  // e.g., "s" (for using transformed distances)

    double rcout_2b = 5.0;
    double rcout_3b = 5.0;
    double rcout_4b = 4.5;

    int nbin_2b = 50;
    int nbin_3b = 50;
    int nbin_4b = 50;

    /////////////////////////////////////////////
    // Read in 2B clusters with composition.
    /////////////////////////////////////////////
    string f1_2b = f1_idx + ".2b_clu-" + style + ".txt";
    string f2_2b = f2_idx + ".2b_clu-" + style + ".txt";
    vector<double> f1_2b_flat_clusters, f2_2b_flat_clusters;
    vector<string> f1_2b_compositions, f2_2b_compositions;
    int npairs_2b = 1;
    read_flat_clusters_with_composition(f1_2b, npairs_2b, f1_2b_flat_clusters, f1_2b_compositions);
    read_flat_clusters_with_composition(f2_2b, npairs_2b, f2_2b_flat_clusters, f2_2b_compositions);

    /////////////////////////////////////////////
    // Read in 3B clusters with composition.
    /////////////////////////////////////////////
    string f1_3b = f1_idx + ".3b_clu-" + style + ".txt";
    string f2_3b = f2_idx + ".3b_clu-" + style + ".txt";
    vector<double> f1_3b_flat_clusters, f2_3b_flat_clusters;
    vector<string> f1_3b_compositions, f2_3b_compositions;
    int npairs_3b = 3;
    read_flat_clusters_with_composition(f1_3b, npairs_3b, f1_3b_flat_clusters, f1_3b_compositions);
    read_flat_clusters_with_composition(f2_3b, npairs_3b, f2_3b_flat_clusters, f2_3b_compositions);

    /////////////////////////////////////////////
    // Read in 4B clusters with composition.
    /////////////////////////////////////////////
    string f1_4b = f1_idx + ".4b_clu-" + style + ".txt";
    string f2_4b = f2_idx + ".4b_clu-" + style + ".txt";
    vector<double> f1_4b_flat_clusters, f2_4b_flat_clusters;
    vector<string> f1_4b_compositions, f2_4b_compositions;
    int npairs_4b = 6;
    read_flat_clusters_with_composition(f1_4b, npairs_4b, f1_4b_flat_clusters, f1_4b_compositions);
    read_flat_clusters_with_composition(f2_4b, npairs_4b, f2_4b_flat_clusters, f2_4b_compositions);

    /////////////////////////////////////////////
    // Determine maximum possible distance between clusters.
    /////////////////////////////////////////////
    double maxd_2b = rcout_2b;                      if (style == "s") maxd_2b = 2.0;
    double maxd_3b = sqrt(3.0 * pow(rcout_3b, 2.0)); if (style == "s") maxd_3b = sqrt(3.0 * pow(2.0, 2.0));
    double maxd_4b = sqrt(6.0 * pow(rcout_4b, 2.0)); if (style == "s") maxd_4b = sqrt(6.0 * pow(2.0, 2.0));
    
    if (my_rank == 0)
        cout << "Setting maximum histogram values: " << maxd_2b << " " << maxd_3b << " " << maxd_4b << endl;

    double binw_2b = maxd_2b / nbin_2b;
    double binw_3b = maxd_3b / nbin_3b;
    double binw_4b = maxd_4b / nbin_4b;
    
    bool same = (f1_idx == f2_idx);

    // Generate histograms using the weighted distances.
    gen_flat_hists(f1_2b_flat_clusters, f2_2b_flat_clusters, f1_2b_compositions, f2_2b_compositions,
                   npairs_2b, nbin_2b, binw_2b, maxd_2b, f1_idx + "-" + f2_idx + ".2b_clu-" + style + ".hist", rcout_2b, same);
    gen_flat_hists(f1_3b_flat_clusters, f2_3b_flat_clusters, f1_3b_compositions, f2_3b_compositions,
                   npairs_3b, nbin_3b, binw_3b, maxd_3b, f1_idx + "-" + f2_idx + ".3b_clu-" + style + ".hist", rcout_3b, same);
    gen_flat_hists(f1_4b_flat_clusters, f2_4b_flat_clusters, f1_4b_compositions, f2_4b_compositions,
                   npairs_4b, nbin_4b, binw_4b, maxd_4b, f1_idx + "-" + f2_idx + ".4b_clu-" + style + ".hist", rcout_4b, same);

    MPI_Finalize();
    return 0;
}
