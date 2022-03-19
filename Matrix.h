#include <string>
#include <vector>
using namespace std;

class Matrix
{
private:
    int number_rows;
    int number_columns;

    vector<vector<int>> index_M;
    vector<vector<int>> index_N;
public:
    Matrix();
    ~Matrix();

    void Inverse(string type_matrix, int number_rows, const vector<vector<float>>& M, vector<vector<float>>& N);
    void Inverse(string type_matrix, int number_rows, const vector<vector<double>>& M, vector<vector<double>>& N);
    void Multiplication(int M_rows, int M_columns, int N_columns, const vector<vector<float>>& M, const vector<vector<float>>& N, vector<vector<float>>& O);
    void Multiplication(int M_rows, int M_columns, int N_columns, const vector<vector<double>>& M, const vector<vector<double>>& N, vector<vector<double>>& O);
    void Transpose(int number_rows, int number_columns, const vector<vector<float>>& M, vector<vector<float>>& N);
    void Transpose(int number_rows, int number_columns, const vector<vector<double>>& M, vector<vector<double>>& N);

    int LU_Decomposition(int number_rows, const vector<vector<float>>& M, vector<vector<float>>& L, vector<vector<float>>& U);
    int LU_Decomposition(int number_rows, const vector<vector<double>>& M, vector<vector<double>>& L, vector<vector<double>>& U);

    float Determinant(string type_matrix, int number_rows, const vector<vector<float>>& M);

    double Determinant(string type_matrix, int number_rows, const vector<vector<double>>& M);
};
