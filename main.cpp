#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

std::vector<double> mvm(std::vector<std::vector<double>> A, std::vector<double> b, int size)
{
  std::vector<double> result(size, 0);
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      result[i] += A[i][j] * b[j];

  return result;
}

double ip(std::vector<double> a, std::vector<double> b, int size)
{
  double result = 0;
  for (int i = 0; i < size; i++)
    result += a[i] * b[i];

  return result;
}

std::vector<std::vector<double>> matrixpow (std::vector<std::vector<double>> A, int n, int size)
{
  for (int i = 0; i < size; i ++)
  {
    for (int j = 0; j < size; j++)
    {
      A[i][j] = pow(A[i][j], n);
    }
  }

  return A;
}

std::vector<double> inputvec(std::string filename, int size)
{
  std::vector<double> input(size, 0);
  std::ifstream ifs(filename);
  std::string str;

  if (ifs.fail())
    std::cout<< "Failded to open vec file." << std::endl;

  int count = 0;
  while (getline(ifs, str))
  {
    input[count] = std::stod(str);
    count++;
    if (size <= count)
      break;
  }

  return input;
}

std::vector<std::vector<double>> inputmatrix(std::string filename, int size)
{
  std::vector<std::vector<double>> input(size, std::vector<double>(size, 0));
  std::ifstream ifs(filename);
  std::string str;

  if (ifs.fail())
    std::cout << "Failded to open matrix file." << std::endl;

  int count1 = 0;
  int count2 = 0;
  while (getline(ifs, str))
  {
    input[count1][count2] = std::stod(str);
    count2++;
    if (size <= count2)
    {
      count2 = 0;
      count1++;
      if (size <= count1)
        break;
    }
  }

  return input;
}

double norm(std::vector<double> vec, int size)
{
  double s = 0.0;

  for (int i = 0; i < size; i++)
  {
    s += vec[i] * vec[i];
  }

  return sqrt(s);
}

int main() {
  int n = 0;
  int size = 64;
  int k = 4;
  double eps = pow(2, -50);

  std::vector<std::vector<double>> A = inputmatrix("../input.txt", size);
  std::vector<double> b = inputvec("../inputb.txt", size);

  std::vector<std::vector<double>> x(size, std::vector<double>(size, 0));

  std::vector<double> Ax = mvm(A, b, size);

  std::vector<std::vector<double>> r(size, std::vector<double>(size, 0));
  // r = b-Ax
  for (int i = 0; i < size; i++)
    r[0][i] = b[i] - Ax[i];

  std::vector<std::vector<double>> p = r;
  // p = r

  std::vector<double> gamma(size, 0);
  std::vector<std::vector<double>> delta(size, std::vector<double>(k+1, 0));
  std::vector<std::vector<double>> eta(size, std::vector<double>(k+2, 0));
  std::vector<std::vector<double>> zeta(size, std::vector<double>(k+3, 0));
  std::vector<double> alpha(size*k, 0);
  std::vector<double> beta(size*k, 0);
  std::vector<double> residual(size, 0);

  while(true)
  {
    std::cout << n << std::endl;

//    residual[step] = norm(r[n], size) / norm(b, size);
    gamma[n] = ip(r[n], r[n], size);
    for (int i = 1; i <= k; i++)
    {
      std::vector<double> Ar = mvm(matrixpow(A, i, size), r[n], size);
      delta[n][i] = ip(r[n], Ar, size);
    }

    for (int i = 1; i <= k+1; i++)
    {
      std::vector<double> Ap = mvm(matrixpow(A, i, size), p[n], size);
      eta[n][i] = ip(r[n], Ap, size);
    }

    for (int i = 1; i <= k+2; i++)
    {
      std::vector<double> Ap = mvm(matrixpow(A, i, size), p[n], size);
      zeta[n][i] = ip(p[n], Ap, size);
    }

    // for i = n to n + k do
    for (int i = n; i < n+k; i++)
    {
      alpha[i] = gamma[i] / zeta[i][1];
//      std::cout << "alpha:" << alpha[i] << std::endl;
//      if (i == 3)
//      {
//        std::cout << "return" << std::endl;
//        return 0;
//      }
      beta[i] = (alpha[i]*zeta[i][2] / zeta[i][1]) - 1;
//      if (i == 3)
//      {
//        for (int l = 0; l < 4; l++)
//        {
//          std::cout << "beta[i]=(" << alpha[l] << "*" << zeta[l][2] << "/" << zeta[l][1] << ")-1" << std::endl;
//        }
//        return 0;
//      }
      gamma[i+1] = beta[i] * gamma[i];
//      if (i == 3) {
//        for (int l = 0; l < 4; l++)
//          std::cout << beta[l] << " * " << gamma[l] << std::endl;
//        return 0;
//      }

      for (int j = 1; j < (2*k) - 2*(i-n); j++)
      {
        delta[i+1][j] = delta[i][j] - (2*alpha[i]*eta[i][j+1]) + (alpha[i]*alpha[i]*zeta[i][j+2]);
        eta[i+1][j] = delta[i+1][j] + (beta[i]*eta[i][j]) - (alpha[i]*beta[i]*zeta[i][j+1]);
        zeta[i+1][j] = eta[i+1][j] + (beta[i]*eta[i][j]) + (beta[i]*beta[i]*zeta[i][j]) - (alpha[i]*beta[i]*zeta[i][j+1]);
      }

      for (int l = 0; l < size; l++)
        x[i+1][l] = x[i][l] + (alpha[i]*p[i][l]);

//      for (int l = 0; l < size; l++)
//        std::cout << "x[" << i+1 << "][" << l << "] = " << x[i+1][l] << std::endl;

      std::vector<double> Ap = mvm(A, p[i], size);
      for (int l = 0; l < size; l++)
        r[i+1][l] = r[i][l] - (alpha[i]*Ap[l]);

//      for (int l = 0; l < size; l++)
//        std::cout << "r[" << i+1 << "][" << l << "] = " << r[i+1][l] << std::endl;

      for (int l = 0; l < size; l++)
        p[i+1][l] = r[i+1][l] + (beta[i]*p[i][l]);

//      for (int l = 0; l < size; l++)
//        std::cout << "p[" << i+1 << "][" << l << "] = " << p[i+1][l] << std::endl;
    }

    n = n+k+1;

    if (size < n)
      break;
  }
}
