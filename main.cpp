#include <iostream>
#include <vector>
using namespace std;


vector<double> substituicoesRetroativas (int n, vector<vector<double>> A, vector<double> b) {

  vector<double> x (n);

  x[n-1] = b[n-1]/A[n-1][n-1];
  for (int i = n-2; i >= 0; i--) {
    double soma = 0;
    for (int j = i+1; j <= n-1; j++) {
      soma = soma + A[i][j]*x[j];
    }
    x[i] = (b[i] - soma)/A[i][i];
  }

  return x;

}

vector<double> eliminacaoDeGauss (int n, vector<vector<double>> A, vector<double> b) {

  for (int k = 0; k <= n-2; k++) {
    for (int i = k+1; i <= n-1; i++) {
      double m = - A[i][k]/A[k][k];
      A[i][k] = 0;
      for (int j = k+1; j <= n-1; j++) {
        A[i][j] = A[i][j] + m*A[k][j];
      }
      b[i] = b[i] + m*b[k];
    }
  }

  return substituicoesRetroativas(n, A, b);

}

vector<vector<double>> getKi(double Li) {

  vector<vector<double>> Ki(2, vector<double>(2));

  Ki[0][0] = 1.0/Li;
  Ki[0][1] = -1.0/Li;
  Ki[1][0] = -1.0/Li;
  Ki[1][1] = 1.0/Li;
  
  Ki[0][0] += Li/3.0;
  Ki[0][1] += Li/6.0;
  Ki[1][0] += Li/6.0;
  Ki[1][1] += Li/3.0;

  return Ki;

} 

int main() {

  int N = 8;
  double y_at0 = 0;
  double y_at1 = 1;

  vector<vector<double>> K(N-1, vector<double>(N-1));
  vector<double> B(N-1);
  vector<double> U(N-1);

  K[0][0] += getKi(1.0/N)[1][1];
  for (int i = 1; i < N-1; i++) {
    vector<vector<double>> Ki = getKi(1.0/N);
    K[i-1][i-1] += Ki[0][0];
    K[i-1][i] += Ki[0][1];
    K[i][i-1] += Ki[1][0];
    K[i][i] += Ki[1][1];
  }
  K[N-2][N-2] += getKi(1.0/N)[0][0];
  B[0] = -y_at0*(-1.0/(1.0/N)+(1.0/N)/6.0);
  B[N-2] = -y_at1*(-1.0/(1.0/N)+(1.0/N)/6.0);

  U = eliminacaoDeGauss(N-1, K, B);
  cout << "Soluções:\n";

  for (int i = 0; i < N-1; i++) {
    cout << "y_" << i+2 << ": " << U[i] << "\n";
  }

}