#include <vector>
using namespace std;

class KMeans
{
private:
    int dimension_data;
    int number_clusters;
public:
    vector<vector<double>> centroid;

    KMeans(int dimension_data, int number_clusters);
    ~KMeans();

    // forgy method
    void Initialize(int number_data, const vector<vector<double>>& data);

    int Classify(const vector<double>& data);

    double Cluster(int number_data, const vector<vector<double>>& data);
};
