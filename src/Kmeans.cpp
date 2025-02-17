#include <cmath>

#include "gmm/Kmeans.h"

KMeans::KMeans(int dimension_data, int number_clusters)
{
    this->dimension_data = dimension_data;
    this->number_clusters = number_clusters;

    centroid.resize(number_clusters);

    for (int i = 0; i < number_clusters; i++){
        centroid[i].resize(dimension_data, 0.0);
    }
}

KMeans::~KMeans()
{
}

void KMeans::Initialize(int number_data, const vector<vector<double>>& data)
{
    for (int i = 0; i < number_clusters; i++){
        int number_sample = number_data / number_clusters;

        for (int j = 0; j < dimension_data; j++){
            double sum = 0;

            for (int k = i * number_sample; k < (i + 1) * number_sample; k++){
                sum += data[k][j];
            }
            centroid[i][j] = sum / number_sample;
        }
    }
}

int KMeans::Classify(const vector<double>& data)
{
    int argmin;

    double min = -1;

    for (int j = 0; j < number_clusters; j++){
        double distance = 0;

        for (int k = 0; k < dimension_data; k++){
            distance += (data[k] - centroid[j][k]) * (data[k] - centroid[j][k]);
        }
        distance = sqrt(distance);

        if (min == -1 || min > distance){
            argmin = j;
            min = distance;
        }
    }
    return argmin;
}

double KMeans::Cluster(int number_data, const vector<vector<double>>& data)
{
    double movements_centroids = 0;

    vector<int> label(number_data, -1);

    vector<double> mean(dimension_data, 0.0);

    for (int i = 0; i < number_data; i++){
        label[i] = Classify(data[i]);
    }

    for (int j = 0; j < number_clusters; j++){
        int number_sample = 0;

        double movements = 0;

        for (int k = 0; k < dimension_data; k++){
            mean[k] = 0;
        }
        for (int i = 0; i < number_data; i++){
            if (label[i] == j){
                for (int k = 0; k < dimension_data; k++){
                    mean[k] += data[i][k];
                }
                number_sample++;
            }
        }
        for (int k = 0; k < dimension_data; k++){
            if (number_sample) mean[k] /= number_sample;

            movements += (centroid[j][k] - mean[k]) * (centroid[j][k] - mean[k]);
            centroid[j][k] = mean[k];
        }
        movements_centroids += sqrt(movements);
    }

    return movements_centroids;
}
