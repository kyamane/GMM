#include <string>
#include <vector>

using namespace std;

class Gaussian_Mixture_Model
{
private:
    string type_covariance;

    int dimension_data;
    int number_gaussian_components;
public:
    vector<double> weight;  // [number_gaussian_components]

    vector<vector<double>> mean;  // [number_gaussian_components][dimension_data]

    vector<vector<double>> diagonal_covariance;  // [number_gaussian_components][dimension_data]
    vector<vector<vector<double>>> covariance;  // [number_gaussian_components][dimension_data][dimension_data]

    Gaussian_Mixture_Model(string type_covariance, int dimension_data, int number_gaussian_components);
    ~Gaussian_Mixture_Model();

    void Initialize(int number_data, const vector<vector<double>>& data);
    void Load_Parameter(string path);
    void Save_Parameter(string path);

    int Classify(const vector<double>& data);

    double Calculate_Likelihood(const vector<double>& data);
    double Calculate_Likelihood(const vector<double>& data, const vector<double>& gaussian_distribution);
    double Expectaion_Maximization(int number_data, const vector<vector<double>>& data);
    double Gaussian_Distribution(const vector<double>& data, int component_index);
    double Gaussian_Distribution(const vector<double>& data, const vector<double>& mean, const vector<double>& diagonal_covariance);
    double Gaussian_Distribution(const vector<double>& data, const vector<double>& mean, const vector<vector<double>>& covariance);
};
